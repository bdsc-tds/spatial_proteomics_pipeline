"""Per-cell quantification of every channel from a cell label mask (any segmentation method).

Writes a feature table with one row per cell: cell_id, centroid (px and µm), area (px and µm²) and
the mean intensity of every image channel (columns named after the markers), after removing each
channel's pixel value offset from the image metadata (+100 on Xenium protein channels, 0 for DAPI,
the morphology stains and every other technology), values below the offset kept negative. The image is read one
channel plane at a time and the mask is processed in row blocks, so memory stays bounded for large
images.

With --qc_mask (Xenium only: aux_outputs/morphology_focus_qc_masks), pixels flagged in the channel's saturation QC
mask are left out of the mean. A cell without any unflagged pixel gets 0.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import tifffile
import zarr

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(script_dir, "../utils/"))
from readwrite import MultiplexImage


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--image", required=True, help="image the mask was computed on (same Y, X)")
    parser.add_argument("--mask", required=True, help="cell label mask (tiled TIFF, 0 = background)")
    parser.add_argument("--pixel_size_fallback", type=float, required=True, help="µm/px used when the image metadata has none")
    parser.add_argument("--out_csv", required=True)
    parser.add_argument("--block_rows", type=int, default=2048, help="mask rows processed at once")
    parser.add_argument("--qc_mask", default="",
                        help="Xenium saturation QC masks (same channels as --image): flagged pixels are left out of "
                             "the mean; empty = every pixel counts")
    parser.add_argument("--threads", type=int, default=1)
    return parser.parse_args()


def blockwise_bincount(mask, n_labels, block_rows, weights=None):
    """Sum of weights (or pixel count) per label, accumulated over row blocks of the mask."""
    total = np.zeros(n_labels, dtype=np.float64)
    for y in range(0, mask.shape[0], block_rows):
        labels = np.asarray(mask[y : y + block_rows]).ravel()
        w = None if weights is None else weights(y, y + labels.size // mask.shape[1]).ravel()
        total += np.bincount(labels, weights=w, minlength=n_labels)
    return total


def blockwise_valid_sums(mask, n_labels, block_rows, plane, valid):
    """Sum of the valid plane values and their pixel count per label, in one pass over the mask.

    valid(y0, y1) returns the boolean validity of rows y0:y1.
    """
    total = np.zeros(n_labels, dtype=np.float64)
    count = np.zeros(n_labels, dtype=np.float64)
    for y in range(0, mask.shape[0], block_rows):
        labels = np.asarray(mask[y : y + block_rows]).ravel()
        y1 = y + labels.size // mask.shape[1]
        # only the valid pixels of the block contribute
        keep = valid(y, y1).ravel()
        total += np.bincount(labels[keep], weights=plane[y:y1].ravel()[keep], minlength=n_labels)
        count += np.bincount(labels[keep], minlength=n_labels)
    return total, count


def main():
    args = parse_args()

    # label mask, read lazily (only the tiles of one row block are decompressed at a time)
    mask = zarr.open(tifffile.imread(args.mask, aszarr=True), mode="r")
    n_labels = max(int(np.asarray(mask[y : y + args.block_rows]).max())
                   for y in range(0, mask.shape[0], args.block_rows)) + 1
    print(f"mask: {mask.shape}, {n_labels - 1} labels")

    # pixel size: image metadata first, config fallback otherwise
    with MultiplexImage(args.image) as img:
        pixel_size = img.pixel_size_um or args.pixel_size_fallback
    print(f"pixel size: {pixel_size:.5f} µm/px")

    # area and centroid of every label from pixel counts and coordinate sums
    area = blockwise_bincount(mask, n_labels, args.block_rows)
    ys = blockwise_bincount(mask, n_labels, args.block_rows,
                            lambda y0, y1: np.broadcast_to(np.arange(y0, y1)[:, None], (y1 - y0, mask.shape[1])))
    xs = blockwise_bincount(mask, n_labels, args.block_rows,
                            lambda y0, y1: np.broadcast_to(np.arange(mask.shape[1])[None, :], (y1 - y0, mask.shape[1])))
    cells = np.flatnonzero(area[1:] > 0) + 1
    with np.errstate(invalid="ignore", divide="ignore"):
        table = pd.DataFrame({
            "cell_id": cells,
            "x_px": xs[cells] / area[cells],
            "y_px": ys[cells] / area[cells],
            "area_px": area[cells].astype(int),
        })
    table["x_um"] = table["x_px"] * pixel_size
    table["y_um"] = table["y_px"] * pixel_size
    table["area_um2"] = table["area_px"] * pixel_size**2

    # mean intensity of every channel, one plane at a time, offset removed (Xenium protein: 100)
    with MultiplexImage(args.image) as img:
        if img.shape[1:] != mask.shape:
            sys.exit(f"image {img.shape[1:]} and mask {mask.shape} differ in size")
        # Xenium saturation QC masks, one plane per image channel (matched by name)
        qc = MultiplexImage(args.qc_mask) if args.qc_mask else None
        if qc is not None and qc.shape != img.shape:
            sys.exit(f"QC mask {qc.shape} and image {img.shape} differ in size")
        intensities = {}
        for c, name in enumerate(img.channel_names):
            # raw plane, offset removed from the mean (mean(raw - offset) = mean(raw) - offset, values below the
            # offset count as negative) so the plane stays in its integer dtype
            plane = img.read_channel(c, args.threads)
            offset = img.channel_offsets[c]
            if qc is None:
                sums = blockwise_bincount(mask, n_labels, args.block_rows, lambda y0, y1: plane[y0:y1])
                intensities[name] = sums[cells] / area[cells] - offset
                print(f"channel {c + 1}/{img.n_channels}: {name} (offset {offset} removed)")
                continue

            # sum and count of the pixels not flagged in the channel's QC mask
            qc_plane = qc.read_channel(name, args.threads)
            sums, n_valid = blockwise_valid_sums(mask, n_labels, args.block_rows, plane,
                                                 lambda y0, y1: qc_plane[y0:y1] == 0)

            # mean over unflagged pixels; a cell without any gets 0, the background level of XOA's
            # background-subtracted signal (NaN would void the whole marker in TACIT's z-score)
            with np.errstate(invalid="ignore", divide="ignore"):
                intensities[name] = np.where(n_valid[cells] > 0, sums[cells] / n_valid[cells] - offset, 0.0)
            n_flagged = area[cells] - n_valid[cells]
            print(f"channel {c + 1}/{img.n_channels}: {name} (offset {offset} removed; {int(n_flagged.sum())} QC-flagged "
                  f"pixels in {int((n_flagged > 0).sum())} cells excluded, {int((n_valid[cells] == 0).sum())} cells "
                  f"fully flagged set to 0)")
            del qc_plane
        if qc is not None:
            qc.close()
    table = pd.concat([table, pd.DataFrame(intensities)], axis=1)

    table.to_csv(args.out_csv, index=False)
    print(f"written: {args.out_csv} ({len(table)} cells x {img.n_channels} channels)")


if __name__ == "__main__":
    main()
