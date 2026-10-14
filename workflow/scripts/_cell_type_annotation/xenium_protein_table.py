"""Write the 10x protein intensities of one Xenium run as a quantify_cells feature table.

Reads the run with read_xenium_sample(gex_only=False) and keeps its "Protein Expression" features
(prot_ prefix removed), for the cells of the run's own segmentation. Writes one row per cell:
cell_id, centroid (px and µm), area (px and µm²) and one column per protein, so that run_tacit.R and
plot_tacit.py read it like a quantify_cells feature table.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import scipy.sparse

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(script_dir, "../utils/"))
from xenium_readwrite import read_xenium_sample, read_xenium_specs


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xenium_dir", required=True, help="Xenium output folder (experiment.xenium, cells.parquet)")
    parser.add_argument("--out_csv", required=True, help="output feature table")
    return parser.parse_args()


def main(args):
    # protein modality of the Xenium cells
    mdata = read_xenium_sample(args.xenium_dir, cells_table=True, anndata=True, gex_only=False)
    prot = mdata["prot"]
    markers = list(prot.var_names.str.removeprefix("prot_"))
    print(f"{prot.n_obs} cells, {len(markers)} proteins: {', '.join(markers)}", flush=True)

    # Xenium morphology pixel size, for the px centroids and areas
    pixel_size = read_xenium_specs(os.path.join(args.xenium_dir, "experiment.xenium"))["pixel_size"]

    # cell metadata in the quantify_cells columns
    obs = prot.obs
    meta = pd.DataFrame({
        "cell_id": obs["cell_id"].astype(str).to_numpy(),
        "x_px": obs["x_centroid"].to_numpy() / pixel_size,
        "y_px": obs["y_centroid"].to_numpy() / pixel_size,
        "area_px": obs["cell_area"].to_numpy() / pixel_size ** 2,
        "x_um": obs["x_centroid"].to_numpy(),
        "y_um": obs["y_centroid"].to_numpy(),
        "area_um2": obs["cell_area"].to_numpy(),
    })

    # one column per protein
    X = prot.X.toarray() if scipy.sparse.issparse(prot.X) else np.asarray(prot.X)
    table = pd.concat([meta, pd.DataFrame(X, columns=markers)], axis=1)

    os.makedirs(os.path.dirname(os.path.abspath(args.out_csv)), exist_ok=True)
    table.to_csv(args.out_csv, index=False)
    print("written:", args.out_csv)


if __name__ == "__main__":
    main(parse_args())
