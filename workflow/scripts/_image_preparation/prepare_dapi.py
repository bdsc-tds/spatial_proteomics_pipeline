"""Write one sample's channel table and DAPI image (input of the alignment to Xenium).

Reads the raw (or background-subtracted) multiplexed image and writes:
  - channels.csv: channel_index, marker_name of every channel
  - dapi.ome.tif: the DAPI channel alone, raw values (used for alignment to Xenium DAPI)
Kept separate from prepare_segmentation_input so that changing the segmentation markers does not
re-run the alignment. Fails if the DAPI channel is missing from the image.
"""

import argparse
import os
import sys

import pandas as pd

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(script_dir, "../utils/"))
from readwrite import MultiplexImage, resolve_pixel_size, write_ome_tiff


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--image", required=True, help="qptiff or OME-TIFF image")
    parser.add_argument("--dapi_channel", required=True, help="channel written to dapi.ome.tif")
    parser.add_argument("--pixel_size_fallback", type=float, required=True, help="µm/px used when the image metadata has none")
    parser.add_argument("--out_channels", required=True)
    parser.add_argument("--out_dapi", required=True)
    parser.add_argument("--threads", type=int, default=1)
    return parser.parse_args()


def main():
    args = parse_args()

    with MultiplexImage(args.image) as img:
        print(f"image: {args.image}\nshape (C, Y, X): {img.shape}, dtype: {img.dtype}")

        # the DAPI channel must exist before anything is written
        if args.dapi_channel not in img.channel_names:
            sys.exit(f"channel {args.dapi_channel!r} not in image; available: {img.channel_names}")

        pixel_size = resolve_pixel_size(img, args.pixel_size_fallback)
        print(f"pixel size: {pixel_size:.5f} µm/px")

        # channel table
        pd.DataFrame({"channel_index": range(img.n_channels), "marker_name": img.channel_names}).to_csv(
            args.out_channels, index=False)

        # DAPI channel alone (raw values), for the alignment to Xenium
        write_ome_tiff(args.out_dapi, [img.read_channel(args.dapi_channel, args.threads)],
                       shape=img.shape[1:], dtype=img.dtype, channel_names=[args.dapi_channel],
                       pixel_size_um=pixel_size)
        print(f"written: {args.out_dapi}")


if __name__ == "__main__":
    main()
