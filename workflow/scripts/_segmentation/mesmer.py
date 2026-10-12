"""Segment cells with Mesmer (deepcell, through SPACEc) on the prepared segmentation input.

The input is the two-channel image written by prepare_segmentation_input: "nuclear" and "membrane", each already
the max projection of the panel's nuclear / membrane channels. SPACEc segments nuclei and whole cells
separately (tiled for large images) and derives the cytoplasm mask. Mesmer is run at the image's real
pixel size (image_mpp).
Only the three label masks are kept, all labelled with the cell id (the nucleus and cytoplasm of a
cell carry that cell's label); per-cell features are computed by quantify_cells on all channels. SPACEc's own feature tables are written to a temporary folder and discarded.
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
    parser.add_argument("--model_dir", required=True, help="folder holding (or receiving) the downloaded Mesmer model")
    parser.add_argument("--out_cells", required=True)
    parser.add_argument("--out_nuclei", required=True)
    parser.add_argument("--out_cytoplasm", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    import skimage.io
    import spacec as sp
    import tensorflow as tf
    import tifffile

    # SPACEc loads the image with skimage.io.imread, which moves 3- or 4-channel stacks to
    # channels-last (read as RGB/RGBA) and then fails its (C, Y, X) check; tifffile keeps (C, Y, X)
    skimage.io.imread = tifffile.imread

    # report whether TensorFlow sees the GPU (Mesmer silently falls back to CPU otherwise)
    print(f"TensorFlow {tf.__version__}, GPUs visible: {tf.config.list_physical_devices('GPU')}")

    # channel names and pixel size of the segmentation input
    with MultiplexImage(args.image) as img:
        channel_names, pixel_size = img.channel_names, img.pixel_size_um
    print(f"channels: {channel_names}, pixel size: {pixel_size} µm/px")
    missing = [c for c in (SEGMENTATION_NUCLEAR, SEGMENTATION_MEMBRANE) if c not in channel_names]
    if missing:
        sys.exit(f"channels {missing} not in {args.image}")
    if pixel_size is None:
        sys.exit(f"{args.image} carries no pixel size; Mesmer needs it (image_mpp)")

    # SPACEc writes feature tables next to its output: keep them in a temporary folder
    out_dir = os.path.dirname(os.path.abspath(args.out_cells))
    with tempfile.TemporaryDirectory(dir=out_dir, prefix="spacec_tmp_") as tmp:
        channel_file = os.path.join(tmp, "channels.txt")
        with open(channel_file, "w") as f:
            f.write("\n".join(channel_names))

        # nuclear and whole-cell segmentation, cytoplasm = whole cell minus nucleus
        result = sp.tl.cell_segmentation(
            file_name=args.image,
            channel_file=channel_file,
            output_dir=tmp,
            output_fname="spacec",
            seg_method="mesmer",
            nuclei_channel=SEGMENTATION_NUCLEAR,
            membrane_channel_list=[SEGMENTATION_MEMBRANE],
            differentiate_nucleus_cytoplasm=True,
            input_format="Multichannel",
            image_mpp=pixel_size,
            model_path=args.model_dir,
            resize_factor=1,
            size_cutoff=0,
        )

    # SPACEc returns None (it does not raise) when segmentation fails
    masks = {"masks": args.out_cells, "masks_nuclei": args.out_nuclei, "masks_cytoplasm": args.out_cytoplasm}
    if result is None or any(result.get(k) is None for k in masks):
        sys.exit("Mesmer segmentation failed, see the SPACEc messages above")

    # SPACEc labels nuclei (separate Mesmer run) and cytoplasm (connected components) independently
    # of the cells: keep the cell mask and the nuclear pixels, the nuclei and cytoplasm masks are
    # rebuilt with the cell ids by save_label_masks
    cells = result["masks"].astype(np.int32, copy=False)
    in_nucleus = result["masks_nuclei"] > 0
    # drop SPACEc's images and masks before building the derived masks (a full-resolution Xenium
    # mask is ~9 GB as int32)
    result.clear()
    gc.collect()

    # write the three masks, all labelled with the cell id
    save_label_masks(cells, in_nucleus, args.out_cells, args.out_nuclei, args.out_cytoplasm)


if __name__ == "__main__":
    main()
