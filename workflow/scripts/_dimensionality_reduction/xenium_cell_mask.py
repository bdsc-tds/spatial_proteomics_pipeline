"""10x's cells of a Xenium Ranger segmentation as a label mask, for the embedding methods.

Reads the cell labels and cell table of the segmentation (its cells.zarr.zip folder, read with
read_xenium_sample) and writes the cell label mask on the morphology image's pixel grid as a tiled TIFF, like
the segmentation methods' masks, with the 10x cell_id of every label (cell_ids.csv). The GPU environment of the
embedding scripts cannot read SpatialData, hence this separate step.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(script_dir, "../utils/"))
sys.path.insert(0, os.path.join(script_dir, "../_segmentation/"))
from label_masks import write_label_mask
from xenium_readwrite import read_xenium_sample


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seg_dir", required=True, help="Xenium Ranger segmentation folder (holds cells.zarr.zip)")
    parser.add_argument("--out_mask", required=True, help="cell label mask (tiled TIFF)")
    parser.add_argument("--out_cell_ids", required=True, help="label -> 10x cell_id table (CSV)")
    return parser.parse_args()


def main():
    args = parse_args()

    # 10x cell labels and the cell table (label of every cell_id)
    sdata = read_xenium_sample(args.seg_dir, cells_labels=True, cells_table=True)
    obs = sdata.tables["table"].obs
    cell_ids = pd.DataFrame({"label": obs["cell_labels"].to_numpy().astype(np.int64),
                             "cell_id": obs["cell_id"].astype(str).to_numpy()})
    print(f"{len(cell_ids)} cells in the table", flush=True)

    # label mask in memory (int32, ~9 GB for a full Xenium slide), written as a tiled TIFF
    labels = np.asarray(sdata.labels["cell_labels"]["scale0"].image.data.compute())
    print(f"cell labels {labels.shape}", flush=True)
    write_label_mask(args.out_mask, labels)
    cell_ids.to_csv(args.out_cell_ids, index=False)
    print(f"written: {args.out_cell_ids}")


if __name__ == "__main__":
    main()
