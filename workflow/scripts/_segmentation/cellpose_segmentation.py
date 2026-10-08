"""Segment cells with Cellpose 3 (through SPACEc) on the prepared segmentation input.

The input is the two-channel image written by prepare_image: "nuclear" and "membrane", each already
the max projection of the panel's nuclear / membrane channels. Cellpose runs twice, both tiled and
stitched by SPACEc for large images: the nuclear model on the nuclear signal alone, and the
whole-cell model on membrane + nuclear (SPACEc's own nucleus/cytoplasm mode would run the whole-cell
model on the nuclei as well). Cellpose works in pixels: the expected diameters are given in µm and
converted with the image's pixel size. A fixed diameter also stops Cellpose from estimating it again
on every tile.
Only the three label masks are kept, all labelled with the cell id (the nucleus and cytoplasm of a
cell carry that cell's label); per-cell features are computed by quantify_cells on all channels.
SPACEc's own feature tables are written to a temporary folder and discarded.
"""

import argparse
import gc
import os
import sys
import tempfile
import warnings

import numpy as np

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(script_dir, "../utils/"))
from label_masks import save_label_masks
from readwrite import SEGMENTATION_MEMBRANE, SEGMENTATION_NUCLEAR, MultiplexImage

warnings.filterwarnings("ignore")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--image", required=True, help=f"segmentation input (OME-TIFF, channels {SEGMENTATION_NUCLEAR!r} and {SEGMENTATION_MEMBRANE!r})")
    parser.add_argument("--model_dir", required=True, help="folder holding (or receiving) the downloaded Cellpose models")
    parser.add_argument("--model", required=True, help="whole-cell Cellpose model, run on membrane + nuclear (e.g. cyto3)")
    parser.add_argument("--nuclear_model", required=True, help="nuclear Cellpose model, run on the nuclear signal (e.g. nuclei)")
    parser.add_argument("--diameter_um", required=True, type=float, help="expected cell diameter (µm)")
    parser.add_argument("--nuclear_diameter_um", required=True, type=float, help="expected nucleus diameter (µm)")
    parser.add_argument("--out_cells", required=True)
    parser.add_argument("--out_nuclei", required=True)
    parser.add_argument("--out_cytoplasm", required=True)
    return parser.parse_args()


def patch_spacec(spacec_segmentation, use_gpu):
    """Fix SPACEc's Cellpose import, make it fail loudly on a failed tile and take its GPU decision from torch."""
    from cellpose import models as cellpose_models

    # SPACEc 0.1.2 imports `from cellpose import models` but calls `cellpose_models`, so every
    # Cellpose run fails with a NameError
    spacec_segmentation.cellpose_models = cellpose_models

    perform = spacec_segmentation._perform_segmentation

    # SPACEc skips a tile whose segmentation fails, leaving a hole in the stitched mask: raise
    # instead (cell_segmentation catches it and returns None)
    def perform_or_raise(*args, **kwargs):
        mask = perform(*args, **kwargs)
        if mask is None:
            raise RuntimeError("Cellpose failed on a tile, see the messages above")
        return mask

    spacec_segmentation._perform_segmentation = perform_or_raise
    # SPACEc asks TensorFlow whether there is a GPU, but Cellpose runs on torch
    spacec_segmentation.setup_gpu = lambda *args, **kwargs: use_gpu


def run_cellpose(sp, image, channel_file, out_dir, model, diameter_px, membrane_channels, pixel_size):
    """One SPACEc Cellpose run (nuclear signal alone, or membrane + nuclear); returns its label mask."""
    result = sp.tl.cell_segmentation(
        file_name=image,
        channel_file=channel_file,
        output_dir=out_dir,
        output_fname=f"cellpose_{model}",
        seg_method="cellpose",
        nuclei_channel=SEGMENTATION_NUCLEAR,
        membrane_channel_list=membrane_channels,
        input_format="Multichannel",
        model=model,
        diameter=diameter_px,
        image_mpp=pixel_size,
        resize_factor=1,
        size_cutoff=0,
    )

    # SPACEc returns None (it does not raise) when segmentation fails
    if result is None or result.get("masks") is None:
        sys.exit(f"Cellpose segmentation with model {model!r} failed, see the SPACEc messages above")
    mask = result["masks"].astype(np.int32, copy=False)
    # drop SPACEc's copies of the image and masks
    result.clear()
    gc.collect()
    return mask


def main():
    args = parse_args()

    # Cellpose reads its model folder from the environment when it is imported (through spacec)
    os.makedirs(args.model_dir, exist_ok=True)
    os.environ["CELLPOSE_LOCAL_MODELS_PATH"] = args.model_dir

    import skimage.io
    import spacec as sp
    import spacec.tools._segmentation as spacec_segmentation
    import tifffile
    import torch

    # SPACEc loads the image with skimage.io.imread, which moves 3- or 4-channel stacks to
    # channels-last (read as RGB/RGBA) and then fails its (C, Y, X) check; tifffile keeps (C, Y, X)
    skimage.io.imread = tifffile.imread

    # report whether torch sees the GPU (Cellpose silently falls back to CPU otherwise)
    use_gpu = torch.cuda.is_available()
    print(f"torch {torch.__version__}, CUDA available: {use_gpu}")
    patch_spacec(spacec_segmentation, use_gpu)

    # channel names and pixel size of the segmentation input
    with MultiplexImage(args.image) as img:
        channel_names, pixel_size = img.channel_names, img.pixel_size_um
    print(f"channels: {channel_names}, pixel size: {pixel_size} µm/px")
    missing = [c for c in (SEGMENTATION_NUCLEAR, SEGMENTATION_MEMBRANE) if c not in channel_names]
    if missing:
        sys.exit(f"channels {missing} not in {args.image}")
    if pixel_size is None:
        sys.exit(f"{args.image} carries no pixel size; Cellpose needs it to convert the diameters to pixels")

    # expected diameters in pixels
    diameter_px = args.diameter_um / pixel_size
    nuclear_diameter_px = args.nuclear_diameter_um / pixel_size
    print(f"diameters: cells {args.diameter_um} µm = {diameter_px:.1f} px, nuclei {args.nuclear_diameter_um} µm = {nuclear_diameter_px:.1f} px")

    # SPACEc writes feature tables next to its output: keep them in a temporary folder
    out_dir = os.path.dirname(os.path.abspath(args.out_cells))
    with tempfile.TemporaryDirectory(dir=out_dir, prefix="spacec_tmp_") as tmp:
        channel_file = os.path.join(tmp, "channels.txt")
        with open(channel_file, "w") as f:
            f.write("\n".join(channel_names))

        # nuclei: nuclear model on the nuclear signal alone; only the nuclear pixels are kept
        nuclei = run_cellpose(sp, args.image, channel_file, tmp, args.nuclear_model, nuclear_diameter_px, None, pixel_size)
        in_nucleus = nuclei > 0
        del nuclei
        gc.collect()

        # whole cells: whole-cell model on membrane + nuclear
        cells = run_cellpose(sp, args.image, channel_file, tmp, args.model, diameter_px, [SEGMENTATION_MEMBRANE], pixel_size)

    # write the three masks, all labelled with the cell id
    save_label_masks(cells, in_nucleus, args.out_cells, args.out_nuclei, args.out_cytoplasm)


if __name__ == "__main__":
    main()
