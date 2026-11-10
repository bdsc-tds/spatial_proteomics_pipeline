"""Local refinement of the insitupy alignment: one residual shift per tile, after the global affine.

A single affine cannot follow parts of the section that moved differently from the rest (e.g. detached tissue
fragments). On the Xenium DAPI and the affine-registered proteomics DAPI (dapi_registered.ome.tif of
align_to_xenium_insitupy, same Xenium pixel grid):
  1. both are downsampled by an integer step to about --work_max_px on the long side, 8-bit (2-98th percentile)
  2. the image is cut into --grid x --grid tiles; each tile is matched in a window extended by half a tile on
     every side (skipped if blank)
  3. per window: SIFT keypoints, Lowe ratio test, 4-parameter similarity (rotation, scale, shift) with RANSAC
  4. the shift of the tile centre under that fit is the tile's control point; control points shifted by more
     than --max_shift_um are dropped (wrong fits: the thin-plate spline would pass through them exactly)
The control points are applied to cell coordinates with multimodal_readwrite.apply_local_refinement (thin-plate
spline through them). Per-tile loop adapted from run_alignment.py (estimate_residual_field) by a colleague.
Writes:
  - control_points.csv: one row per tile with a fit: centre (cx, cy) and shift (dx, dy) in working-grid px,
    scale_x / scale_y (working-grid px per Xenium px), shift_um, n_inliers, kept (False: > --max_shift_um)
  - a QC figure: control points over both DAPIs, and the magnitude of the resulting field
"""

import argparse
import os
import sys

import cv2
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
import zarr

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(script_dir, "../utils/"))
from plotting import savefig, set_mpl_rcparams
from multimodal_readwrite import apply_local_refinement, xenium_pixel_size

# SIFT / matching settings of run_alignment.py
OVERLAP = 0.5      # window = tile + OVERLAP tile on every side
NFEATURES = 10000  # max SIFT keypoints per window and image
RATIO = 0.75       # Lowe ratio test
MIN_MATCHES = 8    # min ratio-test matches and RANSAC inliers per window
MIN_STD = 8        # windows with a lower 8-bit std in either image are blank


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registered", required=True,
                        help="dapi_registered.ome.tif (channel 0 Xenium DAPI, 1 proteomics DAPI after the affine)")
    parser.add_argument("--xenium_dir", required=True, help="Xenium output folder (experiment.xenium, pixel size)")
    parser.add_argument("--grid", type=int, default=12, help="tiles per side")
    parser.add_argument("--work_max_px", type=int, default=6000, help="long side of the downsampled images (px)")
    parser.add_argument("--max_shift_um", type=float, default=100, help="control points shifted more are dropped")
    parser.add_argument("--sample_id", default="", help="used in the figure title")
    parser.add_argument("--out_control_points", required=True)
    parser.add_argument("--out_figure", required=True)
    return parser.parse_args()


def to_uint8(a):
    """2-98 percentile contrast stretch to uint8 (SIFT needs 8-bit input)."""
    p1, p2 = np.percentile(a, (2, 98))
    a = np.clip((a.astype("float32") - p1) / max(1e-6, p2 - p1), 0, 1)
    return (a * 255).astype("uint8")


def tile_shifts(template, image, grid):
    """Per tile: centre (cx, cy), shift (dx, dy) of the centre under the tile's similarity fit (template ->
    image), and the number of RANSAC inliers, in working-grid px. Tiles without a fit are left out."""
    H, W = template.shape
    t_u8, i_u8 = to_uint8(template), to_uint8(image)
    sift = cv2.SIFT_create(nfeatures=NFEATURES)
    matcher = cv2.BFMatcher()
    step_y, step_x = H / grid, W / grid
    half_y, half_x = step_y * (0.5 + OVERLAP), step_x * (0.5 + OVERLAP)
    rows = []
    for iy in range(grid):
        for ix in range(grid):
            # window around the tile centre
            cy, cx = (iy + 0.5) * step_y, (ix + 0.5) * step_x
            y0, y1 = int(max(0, cy - half_y)), int(min(H, cy + half_y))
            x0, x1 = int(max(0, cx - half_x)), int(min(W, cx + half_x))
            win_t, win_i = t_u8[y0:y1, x0:x1], i_u8[y0:y1, x0:x1]
            if win_t.std() < MIN_STD or win_i.std() < MIN_STD:
                continue
            # SIFT keypoints in both windows, kept if the best match is clearly better than the second
            kp_t, desc_t = sift.detectAndCompute(win_t, None)
            kp_i, desc_i = sift.detectAndCompute(win_i, None)
            if desc_t is None or desc_i is None or len(kp_t) < MIN_MATCHES or len(kp_i) < MIN_MATCHES:
                continue
            good = [m for m, n in matcher.knnMatch(desc_t, desc_i, k=2) if m.distance < RATIO * n.distance]
            if len(good) < MIN_MATCHES:
                continue
            # similarity fit with RANSAC (OpenCV default threshold 3 px), refined on the inliers
            src = np.float32([kp_t[m.queryIdx].pt for m in good])
            dst = np.float32([kp_i[m.trainIdx].pt for m in good])
            M, inliers = cv2.estimateAffinePartial2D(src, dst, method=cv2.RANSAC)
            if M is None or inliers is None or int(inliers.sum()) < MIN_MATCHES:
                continue
            # shift of the tile centre under the fit = control point
            centre = np.float32([cx - x0, cy - y0])
            moved = M[:, :2] @ centre + M[:, 2]
            rows.append(dict(cx=cx, cy=cy, dx=float(moved[0] - centre[0]), dy=float(moved[1] - centre[1]),
                             n_inliers=int(inliers.sum())))
    return pd.DataFrame(rows, columns=["cx", "cy", "dx", "dy", "n_inliers"])


def main():
    args = parse_args()
    set_mpl_rcparams()
    pixel_size = xenium_pixel_size(args.xenium_dir)

    # both DAPIs on the working grid: every step-th pixel of the registered image
    with tifffile.TiffFile(args.registered) as tif:
        z = zarr.open(tif.aszarr(series=0, level=0), mode="r")
        _, full_h, full_w = z.shape
        step = max(1, max(full_h, full_w) // args.work_max_px)
        template = np.asarray(z[0, ::step, ::step])
        image = np.asarray(z[1, ::step, ::step])
    work_h, work_w = template.shape
    work_um = step * pixel_size
    print(f"registered {full_h} x {full_w} px -> working grid {work_h} x {work_w} px (step {step}, {work_um:.2f} µm/px)")

    # control points, with the max-shift guard
    cps = tile_shifts(template, image, args.grid)
    cps["scale_x"], cps["scale_y"] = work_w / full_w, work_h / full_h
    cps["shift_um"] = np.hypot(cps["dx"], cps["dy"]) * work_um
    cps["kept"] = cps["shift_um"] <= args.max_shift_um
    cps.to_csv(args.out_control_points, index=False)
    print(f"{len(cps)} of {args.grid ** 2} tiles fitted, {(~cps['kept']).sum()} dropped (> {args.max_shift_um} µm)")
    for r in cps[~cps["kept"]].itertuples():
        print(f"  dropped: centre ({r.cx / r.scale_x:.0f}, {r.cy / r.scale_y:.0f}) Xenium px, {r.shift_um:.1f} µm, "
              f"{r.n_inliers} inliers")
    kept = cps[cps["kept"]]
    if len(kept):
        print(f"kept shifts (µm): median {kept['shift_um'].median():.2f}, max {kept['shift_um'].max():.2f}")
    if len(kept) < 4:
        print("fewer than 4 control points kept: the refinement leaves the cells where the affine put them")

    # field magnitude on a coarse grid (Xenium px), from the same function the readers use
    gx, gy = np.meshgrid(np.linspace(0, full_w - 1, 150), np.linspace(0, full_h - 1, 300))
    grid_px = np.column_stack([gx.ravel(), gy.ravel()])
    shift = np.linalg.norm(apply_local_refinement(grid_px, args.out_control_points) - grid_px, axis=1)
    shift_um = (shift * pixel_size).reshape(gx.shape)

    # QC figure: control points over both DAPIs (left), field magnitude (right)
    ext = [0, full_w * pixel_size, full_h * pixel_size, 0]
    stretch = lambda a: to_uint8(a).astype(np.float32) / 255
    fig, axes = plt.subplots(1, 2, figsize=(12, 12 * work_h / work_w / 2 + 1.5), layout="constrained")
    axes[0].imshow(np.stack([stretch(template), stretch(image), stretch(template)], axis=-1), extent=ext)
    for sel, color, label in [(cps["kept"], "cyan", "kept"), (~cps["kept"], "red", f"> {args.max_shift_um:g} µm, dropped")]:
        c = cps[sel]
        x, y = c["cx"] / c["scale_x"] * pixel_size, c["cy"] / c["scale_y"] * pixel_size
        axes[0].plot(x, y, "o", ms=3, color=color, label=f"{label} ({len(c)})")
        if color == "cyan":
            # shifts x20 (µm)
            axes[0].quiver(x, y, c["dx"] * work_um, c["dy"] * work_um, color=color, angles="xy",
                           scale_units="xy", scale=0.05, width=0.003)
    axes[0].legend(loc="lower left", fontsize=7)
    axes[0].set_title("Xenium DAPI (magenta), proteomics DAPI after the affine (green)\ncontrol points, shifts x20", fontsize=9)
    im = axes[1].imshow(shift_um, extent=ext, cmap="magma", vmin=0)
    fig.colorbar(im, ax=axes[1], shrink=0.6, label="local shift (µm)")
    axes[1].set_title(f"shift applied by the refinement ({len(kept)} control points)", fontsize=9)
    for ax in axes:
        ax.set_xlabel("x (µm)")
        ax.set_ylabel("y (µm)")
    fig.suptitle(f"Local refinement {args.sample_id}: {args.grid} x {args.grid} tiles, {work_um:.2f} µm/px", fontsize=10)
    savefig(None, args.out_figure, dpi=150)
    print(f"written: {args.out_control_points}, {args.out_figure}")


if __name__ == "__main__":
    main()
