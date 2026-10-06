"""Register the proteomics DAPI image onto the Xenium nuclei (DAPI) image with insitupy.

The Xenium DAPI (morphology focus channel 0) is read with read_xenium_sample, as in analysis_repo.
Uses insitupy's ImageRegistration directly: SIFT features are matched on downscaled copies of both
images and the affine matrix is fitted on the full-resolution feature coordinates. Only that matrix
is saved (3x3, maps proteomics image pixels (x, y) to Xenium morphology image pixels), together with
the feature-matching QC figure; no registered image is written.

If the images only match after flipping the proteomics image, insitupy keeps that flip separate from
its matrix; here it is folded into the saved matrix so that it applies to the unflipped image.
If too few features match (samples not from the same section, or too different), the step fails;
lowering --min_good_matches_per_area allows alignment with fewer matches.
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(script_dir, "../utils/"))
from plotting import savefig, set_mpl_rcparams
from readwrite import MultiplexImage
from xenium_readwrite import read_xenium_sample, read_xenium_specs


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xenium_dir", required=True, help="Xenium output folder (holding experiment.xenium)")
    parser.add_argument("--dapi", required=True, help="proteomics DAPI image (dapi.ome.tif from prepare_image)")
    parser.add_argument("--out_affine", required=True)
    parser.add_argument("--out_matches", required=True, help="QC figure of the matched features")
    parser.add_argument("--min_good_matches_per_area", type=int, default=5, help="minimum good matches per mm² of Xenium image")
    parser.add_argument("--max_width", default="4000",
                        help="both images are shrunk to about max_width² px before feature matching; "
                             "'codex': the proteomics image size, so only the Xenium image is shrunk")
    return parser.parse_args()


def flip_matrix(flip_axis, shape):
    """3x3 matrix of the flip insitupy applied to the (Y, X) image before matching."""
    height, width = shape
    if flip_axis == 0:  # vertical flip: y -> height - 1 - y
        return np.array([[1, 0, 0], [0, -1, height - 1], [0, 0, 1]], dtype=float)
    if flip_axis == 1:  # horizontal flip: x -> width - 1 - x
        return np.array([[-1, 0, width - 1], [0, 1, 0], [0, 0, 1]], dtype=float)
    return np.eye(3)


def main():
    args = parse_args()
    set_mpl_rcparams()
    from insitupy.tools.registration import ImageRegistration
    from spatialdata import get_pyramid_levels
    from xarray import DataTree

    # Xenium DAPI: channel 0 of the morphology focus image, full resolution
    sdata = read_xenium_sample(args.xenium_dir, morphology_focus=True)
    focus = sdata.images["morphology_focus"]
    if isinstance(focus, DataTree):
        focus = get_pyramid_levels(focus, n=0)
    xenium_dapi = focus.isel(c=0)
    template = xenium_dapi.values

    # Xenium pixel size from experiment.xenium
    xenium_pixel_size = read_xenium_specs(os.path.join(args.xenium_dir, "experiment.xenium"))["pixel_size"]
    print(f"Xenium {xenium_dapi.c.item()} image: {template.shape}, {xenium_pixel_size} µm/px")

    # proteomics DAPI plane
    with MultiplexImage(args.dapi) as img:
        dapi = img.read_channel(0)
    print(f"proteomics DAPI: {dapi.shape}, {img.pixel_size_um} µm/px")

    # minimum number of good matches, scaled with the Xenium image area (as insitupy.register_images)
    area_mm2 = template.shape[0] * template.shape[1] * xenium_pixel_size**2 / 1000**2
    min_good_matches = int(args.min_good_matches_per_area * area_mm2)
    print(f"minimum good matches: {min_good_matches} ({args.min_good_matches_per_area}/mm² x {area_mm2:.1f} mm²)")

    # max_width: insitupy scales each image to max_width² px area, so sqrt of the proteomics area keeps it unscaled (to within 1 px)
    max_width = round(np.sqrt(dapi.shape[0] * dapi.shape[1])) if args.max_width == "codex" else int(args.max_width)

    # feature matching and affine fit (full-resolution coordinates)
    print(f"feature matching on copies shrunk to ~{max_width}² px")
    reg = ImageRegistration(
        image=dapi,
        template=template,
        axes_image="YX",
        axes_template="YX",
        max_width=max_width,
        convert_to_grayscale=False,
        perspective_transform=False,
        min_good_matches=min_good_matches,
    )
    reg.load_and_scale_images()
    print(f"scaled for matching: proteomics {reg.image_scaled.shape}, Xenium {reg.template_scaled.shape}")
    reg.extract_features(test_flipping=True)
    reg.calculate_transformation_matrix()

    # 3x3 matrix for the unflipped proteomics image
    flip_axis = getattr(reg, "flip_axis", None)
    matrix = np.vstack([reg.T, [0, 0, 1]]) @ flip_matrix(flip_axis, dapi.shape)
    if flip_axis is not None:
        print(f"image matched after a {'vertical' if flip_axis == 0 else 'horizontal'} flip; flip folded into the matrix")
    print(f"affine matrix (proteomics px -> Xenium px):\n{matrix}")
    np.savetxt(args.out_affine, matrix, delimiter=",")

    # QC figure of the matched features, shrunk to at most 4000 px wide (full size with --max_width codex)
    vis = reg.matchedVis
    step = max(1, int(np.ceil(vis.shape[1] / 4000)))
    plt.figure(figsize=(12, 6))
    plt.imshow(vis[::step, ::step])
    plt.axis("off")
    savefig(None, args.out_matches, dpi=200)
    print(f"written: {args.out_affine}, {args.out_matches}")


if __name__ == "__main__":
    main()
