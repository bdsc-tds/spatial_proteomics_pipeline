"""Bring the proteomics cells into Xenium coordinates, match them to Xenium cells and plot the result.

The proteomics centroids (pixels of the proteomics image) are mapped with the alignment matrix to
Xenium morphology pixels, then to µm with the Xenium pixel size from experiment.xenium. Cells are
matched as mutual nearest neighbours within --max_dist µm. Writes:
  - the table of matched cells (proteomics cell_id, Xenium cell_id, distance in µm)
  - an overlay of Xenium and aligned proteomics centroids
  - a bar plot of the number of Xenium, proteomics, matched and unmatched proteomics cells
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from plotnine import (aes, element_rect, element_text, geom_col, geom_text, ggplot, scale_fill_manual,
                      scale_y_continuous, theme, theme_classic, xlab, ylab)

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(script_dir, "../utils/"))
from utils import apply_affine, mutual_nearest_neighbours
from xenium_readwrite import read_xenium_sample, read_xenium_specs


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--feature_table", required=True, help="quantify_cells output (cell_id, x_px, y_px, ...)")
    parser.add_argument("--affine", required=True, help="3x3 matrix, proteomics px -> Xenium px")
    parser.add_argument("--xenium_dir", required=True, help="Xenium output folder (cells table, experiment.xenium)")
    parser.add_argument("--max_dist", type=float, default=5, help="max distance (µm) between matched cells")
    parser.add_argument("--sample_id", default="", help="used in the plot titles")
    parser.add_argument("--out_matches", required=True)
    parser.add_argument("--out_overlay", required=True)
    parser.add_argument("--out_counts", required=True)
    return parser.parse_args()


def main():
    args = parse_args()

    # proteomics centroids in Xenium µm
    xenium_pixel_size = read_xenium_specs(os.path.join(args.xenium_dir, "experiment.xenium"))["pixel_size"]
    cells = pd.read_csv(args.feature_table, usecols=["cell_id", "x_px", "y_px"])
    matrix = np.loadtxt(args.affine, delimiter=",")
    aligned_um = apply_affine(cells[["x_px", "y_px"]].to_numpy(), matrix) * xenium_pixel_size
    print(f"{len(cells)} proteomics cells, Xenium pixel size {xenium_pixel_size} µm/px")

    # Xenium centroids (µm), from the cells table read as in analysis_repo
    xenium = read_xenium_sample(args.xenium_dir, anndata=True).obs[["cell_id", "x_centroid", "y_centroid"]]
    xenium_um = xenium[["x_centroid", "y_centroid"]].to_numpy()
    print(f"{len(xenium)} Xenium cells")

    # mutual nearest neighbours within max_dist
    pairs = mutual_nearest_neighbours(aligned_um, xenium_um, args.max_dist)
    matches = pd.DataFrame({
        "proteomics_cell_id": cells["cell_id"].to_numpy()[pairs["index_a"]],
        "xenium_cell_id": xenium["cell_id"].to_numpy()[pairs["index_b"]],
        "dist_um": pairs["dist"].to_numpy(),
    })
    matches.to_csv(args.out_matches, index=False)
    print(f"{len(matches)} matched cells within {args.max_dist} µm")

    # overlay of the centroids (many points: matplotlib)
    plt.figure(figsize=(8, 8))
    plt.scatter(xenium_um[:, 0], xenium_um[:, 1], s=0.5, alpha=0.4, label="Xenium", color="#19618A")
    plt.scatter(aligned_um[:, 0], aligned_um[:, 1], s=0.5, alpha=0.4, label="proteomics (aligned)", color="#09B39C")
    plt.gca().set_aspect("equal")
    plt.gca().invert_yaxis()
    plt.title(f"Aligned centroids {args.sample_id}", fontsize=9)
    plt.xlabel("x (µm)")
    plt.ylabel("y (µm)")
    plt.legend(markerscale=10)
    plt.savefig(args.out_overlay, dpi=200, bbox_inches="tight")
    plt.close()

    # counts of Xenium, proteomics, matched and unmatched proteomics cells
    counts = pd.DataFrame({
        "modality": ["Xenium", "proteomics", "matched", "unmatched"],
        "ncells": [len(xenium), len(cells), len(matches), len(cells) - len(matches)],
    })
    counts["modality"] = pd.Categorical(counts["modality"], categories=counts["modality"], ordered=True)
    print(counts.to_string(index=False))
    colors = {"Xenium": "#19618A", "proteomics": "#09B39C", "matched": "#FFB300", "unmatched": "#D32F2F"}
    p = (
        ggplot(counts, aes(x="modality", y="ncells", fill="modality"))
        + geom_col()
        + geom_text(aes(label="ncells"), va="bottom", nudge_y=counts["ncells"].max() * 0.01, size=9,
                    format_string="{:,}")
        + theme_classic()
        + theme(panel_border=element_rect(fill=None, color="black"),
                axis_text_x=element_text(angle=45, hjust=1, color="black"),
                legend_position="none",
                figure_size=(4, 4))
        # headroom above the tallest bar for its label
        + scale_y_continuous(expand=[0, 0, 0.01, 0], limits=(0, counts["ncells"].max() * 1.1))
        + scale_fill_manual(values=colors)
        + xlab("")
        + ylab("Number of cells")
    )
    p.save(args.out_counts, dpi=300, verbose=False)
    print(f"written: {args.out_matches}, {args.out_overlay}, {args.out_counts}")


if __name__ == "__main__":
    main()
