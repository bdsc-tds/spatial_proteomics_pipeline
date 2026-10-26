"""Channels VirTues reads, resampled to its resolution, clipped and log-transformed as spora standardizes.

For every channel with a UniProt ID in the panel's marker table: the full-resolution plane with its pixel
value offset removed (+100 on Xenium protein channels) is resampled to resolution_um by area averaging
(cv2.INTER_AREA), clipped at its upper_quantile inside the tissue and log1p-transformed. The tissue is the
smoothed nuclear signal (channels mapped to Histone H3, P68431) above its Otsu threshold, dilated by 20 µm.

Writes the (channels, Y, X) float16 stack as .npy and, per channel, the clipping quantile and the pixel
count, sum and sum of squares of the transformed values inside the tissue (for the panel's mean / std).
"""

import argparse
import os
import sys

import cv2
import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from skimage.filters import gaussian, threshold_otsu
from skimage.morphology import disk

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(script_dir, "../utils/"))
from readwrite import MultiplexImage, resolve_pixel_size

# UniProt ID of the nuclear stains (DAPI, Hoechst, Iridium), as in spora
NUCLEAR_UNIPROT = "P68431"


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--image", required=True, help="multiplexed image (qptiff or OME-TIFF)")
    parser.add_argument("--markers", required=True, help="marker table of the panel (channel_name, uniprot_id)")
    parser.add_argument("--pixel_size_fallback", type=float, required=True, help="µm/px used when the image metadata has none")
    parser.add_argument("--resolution_um", type=float, default=1.0, help="µm/px VirTues reads")
    parser.add_argument("--upper_quantile", type=float, default=0.99, help="per-channel clipping quantile inside the tissue")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--out_image", required=True, help="(channels, Y, X) float16 .npy")
    parser.add_argument("--out_stats", required=True, help="per-channel statistics (CSV)")
    return parser.parse_args()


def tissue_mask(nuclear, resolution_um):
    """Tissue pixels: smoothed log nuclear signal above its Otsu threshold, dilated by 20 µm, holes filled."""
    smooth = gaussian(np.log1p(nuclear), sigma=5 / resolution_um)
    mask = smooth > threshold_otsu(smooth)
    mask = ndi.binary_dilation(mask, structure=disk(round(20 / resolution_um)))
    return ndi.binary_fill_holes(mask)


def main():
    args = parse_args()
    cv2.setNumThreads(args.threads)

    # channels VirTues reads: those with a UniProt ID, in marker table order
    markers = pd.read_csv(args.markers, dtype=str)
    markers = markers[markers["uniprot_id"].notna()].reset_index(drop=True)
    markers["uniprot_id"] = markers["uniprot_id"].str.strip()
    print(f"{len(markers)} channels with a UniProt ID in {args.markers}", flush=True)

    # each channel resampled to resolution_um by area averaging, offsets removed
    with MultiplexImage(args.image) as img:
        pixel_size = resolve_pixel_size(img, args.pixel_size_fallback)
        scale = pixel_size / args.resolution_um
        height, width = img.shape[1:]
        out_shape = (round(height * scale), round(width * scale))
        print(f"image {img.shape[1:]} at {pixel_size:.4f} µm/px -> {out_shape} at {args.resolution_um} µm/px", flush=True)
        stack = np.empty((len(markers), *out_shape), dtype=np.float32)
        for c, name in enumerate(markers["channel_name"]):
            plane = img.read_channel(name, args.threads, subtract_offset=True).astype(np.float32)
            stack[c] = cv2.resize(plane, (out_shape[1], out_shape[0]), interpolation=cv2.INTER_AREA)
            del plane
            print(f"channel {c + 1}/{len(markers)}: {name}", flush=True)

    # tissue from the nuclear channels (whole image when the panel maps none)
    nuclear = markers["uniprot_id"] == NUCLEAR_UNIPROT
    if nuclear.any():
        tissue = tissue_mask(stack[nuclear.to_numpy()].max(axis=0), args.resolution_um)
    else:
        print("WARNING: no channel mapped to Histone H3 (P68431): statistics over the whole image")
        tissue = np.ones(out_shape, dtype=bool)
    print(f"tissue: {tissue.mean() * 100:.1f}% of the image", flush=True)

    # clip each channel at its upper quantile inside the tissue, log1p, and sum the tissue values
    rows = []
    for c, row in markers.iterrows():
        upper = float(np.quantile(stack[c][tissue], args.upper_quantile))
        np.clip(stack[c], 0, upper, out=stack[c])
        np.log1p(stack[c], out=stack[c])
        values = stack[c][tissue].astype(np.float64)
        rows.append({"channel_name": row["channel_name"], "uniprot_id": row["uniprot_id"], "upper_quantile_value": upper,
                     "n_pixels": values.size, "sum": values.sum(), "sum_sq": np.square(values).sum()})

    np.save(args.out_image, stack.astype(np.float16))
    pd.DataFrame(rows).to_csv(args.out_stats, index=False)
    print(f"written: {args.out_image} {stack.shape}, {args.out_stats}")


if __name__ == "__main__":
    main()
