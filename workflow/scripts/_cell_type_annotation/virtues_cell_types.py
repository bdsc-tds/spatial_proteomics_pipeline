"""VirTues cell type of every cell of a mask: the majority class of its pixels in the VirTues class map.

The class map (segmentation_virtues) holds the most likely cell type of every pixel of the image resampled
to 1 µm/px; each full-resolution pixel of the mask takes the class of the 1 µm pixel it falls in, so small
cells keep all their pixels. The mask is either a cell label mask of a segmentation method (tiled TIFF) or
10x's cell labels of a Xenium Ranger segmentation (its cells.zarr.zip, read with read_xenium_sample).

Writes one row per cell: cell_id (mask label, or 10x cell_id), cell_type, confidence (pixel fraction of
that class) and the pixel fraction of every class.
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


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--class_map", required=True, help="VirTues class map (tiled TIFF, 1 µm/px)")
    parser.add_argument("--classes", required=True, help="class index -> cell type (classes.csv)")
    parser.add_argument("--mask", required=True, help="cell label mask (tiled TIFF) or a Xenium cells.zarr.zip")
    parser.add_argument("--block_rows", type=int, default=2048, help="mask rows processed at once")
    parser.add_argument("--out_csv", required=True)
    return parser.parse_args()


def read_mask(path):
    """Lazily read label mask and, for Xenium, the 10x cell_id of every label (None for our masks)."""
    if os.path.basename(path) == "cells.zarr.zip":
        from xenium_readwrite import read_xenium_sample
        sdata = read_xenium_sample(os.path.dirname(path), cells_labels=True, cells_table=True)
        labels = sdata.labels["cell_labels"]["scale0"].image.data
        obs = sdata.tables["table"].obs
        cell_ids = pd.Series(obs["cell_id"].astype(str).to_numpy(), index=obs["cell_labels"].to_numpy())
        return labels, cell_ids
    return zarr.open(tifffile.imread(path, aszarr=True), mode="r"), None


def main():
    args = parse_args()

    # class map (small at 1 µm/px) and the classes it can hold
    class_map = tifffile.imread(args.class_map)
    # keep_default_na=False: the "None" class must stay a string, not NaN
    classes = pd.read_csv(args.classes, keep_default_na=False)
    n_classes = int(classes["class_index"].max()) + 1

    # mask and its largest label (Xenium: from the cell table, saving a pass over the mask)
    mask, cell_ids = read_mask(args.mask)
    height, width = mask.shape
    blocks = range(0, height, args.block_rows)
    if cell_ids is not None:
        n_labels = int(cell_ids.index.max()) + 1
    else:
        n_labels = max(int(np.asarray(mask[y : y + args.block_rows]).max()) for y in blocks) + 1
    print(f"mask {mask.shape}, {n_labels - 1} labels; class map {class_map.shape}", flush=True)

    # 1 µm pixel each full-resolution row / column falls in
    rows = np.minimum((np.arange(height) * class_map.shape[0] / height).astype(np.int64), class_map.shape[0] - 1)
    cols = np.minimum((np.arange(width) * class_map.shape[1] / width).astype(np.int64), class_map.shape[1] - 1)

    # pixels of every (label, class) pair, one row block at a time
    counts = np.zeros(n_labels * n_classes, dtype=np.int64)
    for y in blocks:
        labels = np.asarray(mask[y : y + args.block_rows]).astype(np.int64)
        block_classes = class_map[rows[y : y + labels.shape[0]]][:, cols]
        counts += np.bincount((labels * n_classes + block_classes).ravel(), minlength=counts.size)
    counts = counts.reshape(n_labels, n_classes)[:, classes["class_index"].to_numpy()]

    # majority class per cell, its pixel fraction as confidence
    cells = np.flatnonzero(counts[1:].sum(axis=1) > 0) + 1
    fractions = counts[cells] / counts[cells].sum(axis=1, keepdims=True)
    names = classes["cell_type"].to_numpy()
    table = pd.DataFrame({
        "cell_id": cell_ids.loc[cells].to_numpy() if cell_ids is not None else cells,
        "cell_type": names[fractions.argmax(axis=1)],
        "confidence": fractions.max(axis=1),
    })
    table = pd.concat([table, pd.DataFrame(fractions, columns=[f"fraction_{n}" for n in names])], axis=1)

    table.to_csv(args.out_csv, index=False)
    print(table["cell_type"].value_counts().to_string())
    print(f"written: {args.out_csv} ({len(table)} cells)")


if __name__ == "__main__":
    main()
