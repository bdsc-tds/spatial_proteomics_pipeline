"""Import the 10x cells of a Xenium segmentation, warped onto the proteomics image grid outside the pipeline.

Checks the warped cell, nuclei and cytoplasm label masks against the image the pipeline quantifies (same Y, X;
int32) and their label -> 10x cell_id table (unique labels, every label of the cell mask listed), then copies the
four files into the method's segmentation folder. The masks are not modified: they are quantified like any other
segmentation method's masks.
"""

import argparse
import os
import shutil
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
    parser.add_argument("--image", required=True, help="image the masks must match (raw or background-subtracted)")
    parser.add_argument("--cells", required=True, help="warped cell label mask (tiled TIFF, 0 = background)")
    parser.add_argument("--nuclei", required=True, help="warped nucleus label mask")
    parser.add_argument("--cytoplasm", required=True, help="warped cytoplasm label mask")
    parser.add_argument("--cell_ids", required=True, help="label -> 10x cell_id table (CSV: label, cell_id)")
    parser.add_argument("--out_cells", required=True)
    parser.add_argument("--out_nuclei", required=True)
    parser.add_argument("--out_cytoplasm", required=True)
    parser.add_argument("--out_cell_ids", required=True)
    parser.add_argument("--block_rows", type=int, default=2048, help="mask rows read at once")
    return parser.parse_args()


def mask_labels(path, block_rows):
    """Labels (> 0) present in a label mask, read in row blocks."""
    # mask read lazily, one row block of tiles at a time
    mask = zarr.open(tifffile.imread(path, aszarr=True), mode="r")
    labels = set()
    for y in range(0, mask.shape[0], block_rows):
        labels.update(np.unique(np.asarray(mask[y : y + block_rows])).tolist())
    labels.discard(0)
    return labels


def main():
    args = parse_args()

    # Y, X of the image the masks are quantified on
    with MultiplexImage(args.image) as img:
        image_yx = tuple(img.shape[1:])
    print(f"image: {args.image}, Y, X = {image_yx}", flush=True)

    # every mask on the image grid, as int32 labels
    for path in (args.cells, args.nuclei, args.cytoplasm):
        with tifffile.TiffFile(path) as tif:
            shape, dtype = tuple(tif.series[0].shape), tif.series[0].dtype
        print(f"{path}: shape {shape}, {dtype}", flush=True)
        if shape != image_yx:
            sys.exit(f"{path}: shape {shape} differs from the image {image_yx}")
        if dtype != np.int32:
            sys.exit(f"{path}: dtype {dtype}, expected int32")

    # label -> 10x cell_id table: unique labels, covering every cell of the cell mask
    cell_ids = pd.read_csv(args.cell_ids)
    if not {"label", "cell_id"} <= set(cell_ids.columns):
        sys.exit(f"{args.cell_ids}: needs columns label, cell_id (has {list(cell_ids.columns)})")
    if cell_ids["label"].duplicated().any():
        sys.exit(f"{args.cell_ids}: {int(cell_ids['label'].duplicated().sum())} duplicated labels")
    cells = mask_labels(args.cells, args.block_rows)
    unlisted = cells - set(cell_ids["label"].tolist())
    if unlisted:
        sys.exit(f"{len(unlisted)} labels of the cell mask are not in {args.cell_ids}")
    print(f"{len(cells)} cells in the mask, {len(cell_ids)} rows in the cell_id table "
          f"({len(cell_ids) - len(cells)} cells not on the image)", flush=True)

    # copy of the four files into the method's segmentation folder
    for src, dst in ((args.cells, args.out_cells), (args.nuclei, args.out_nuclei),
                     (args.cytoplasm, args.out_cytoplasm), (args.cell_ids, args.out_cell_ids)):
        shutil.copyfile(src, dst)
        print(f"copied: {src} -> {dst}", flush=True)


if __name__ == "__main__":
    main()
