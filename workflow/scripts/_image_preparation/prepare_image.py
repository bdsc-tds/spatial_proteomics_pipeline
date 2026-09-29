"""Prepare one sample's image for the technology-agnostic downstream steps.

Reads the raw (or background-subtracted) multiplexed image one channel at a time and writes:
  - channels.csv:              channel_index, marker_name of every channel
  - dapi.ome.tif:              the nuclear channel (used for alignment to Xenium DAPI)
  - segmentation_input.ome.tif: nuclear + membrane channels, in that order (Mesmer input)
Fails if the configured nuclear or membrane channels are missing from the image.
"""

import argparse
import os
import sys

import pandas as pd

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(script_dir, "../utils/"))
from readwrite import MultiplexImage, write_ome_tiff


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--image", required=True, help="qptiff or OME-TIFF image")
    parser.add_argument("--nuclear_channel", required=True)
    parser.add_argument("--membrane_channels", nargs="+", required=True)
    parser.add_argument("--pixel_size_fallback", type=float, required=True, help="µm/px used when the image metadata has none")
    parser.add_argument("--out_channels", required=True)
    parser.add_argument("--out_dapi", required=True)
    parser.add_argument("--out_segmentation_input", required=True)
    parser.add_argument("--threads", type=int, default=1)
    return parser.parse_args()


def main():
    args = parse_args()

    with MultiplexImage(args.image) as img:
        print(f"image: {args.image}\nshape (C, Y, X): {img.shape}, dtype: {img.dtype}")

        # all requested channels must exist before anything is written
        wanted = [args.nuclear_channel] + args.membrane_channels
        missing = [c for c in wanted if c not in img.channel_names]
        if missing:
            sys.exit(f"channels {missing} not in image; available: {img.channel_names}")

        # pixel size: image metadata first, config fallback otherwise
        pixel_size = img.pixel_size_um or args.pixel_size_fallback
        if img.pixel_size_um is None:
            print(f"no pixel size in image metadata, using fallback {pixel_size} µm/px")
        elif abs(img.pixel_size_um - args.pixel_size_fallback) / args.pixel_size_fallback > 0.01:
            print(f"WARNING: image pixel size {img.pixel_size_um:.5f} µm/px differs from the configured "
                  f"{args.pixel_size_fallback} µm/px; using the image value")
        print(f"pixel size: {pixel_size:.5f} µm/px")

        # channel table
        pd.DataFrame({"channel_index": range(img.n_channels), "marker_name": img.channel_names}).to_csv(
            args.out_channels, index=False)

        # nuclear channel alone, for the alignment to Xenium
        write_ome_tiff(args.out_dapi, [img.read_channel(args.nuclear_channel, args.threads)],
                       shape=img.shape[1:], dtype=img.dtype, channel_names=[args.nuclear_channel],
                       pixel_size_um=pixel_size)
        print(f"written: {args.out_dapi}")

        # nuclear + membrane channels, one plane in memory at a time
        planes = (img.read_channel(c, args.threads) for c in wanted)
        write_ome_tiff(args.out_segmentation_input, planes, shape=(len(wanted),) + img.shape[1:],
                       dtype=img.dtype, channel_names=wanted, pixel_size_um=pixel_size)
        print(f"written: {args.out_segmentation_input} ({', '.join(wanted)})")


if __name__ == "__main__":
    main()
