"""Bring the proteomics cells into Xenium coordinates, match them to Xenium cells and plot the result.

The proteomics centroids (pixels of the proteomics image) are mapped with the alignment matrix to
Xenium morphology pixels, then to µm with the Xenium pixel size from experiment.xenium. Cells are
matched as mutual nearest neighbours within --max_dist µm. Reading, alignment and matching use the
helpers shared with analysis_repo (multimodal_readwrite, xenium_readwrite), so the matches are
the same as those of read_xenium_proteomics_sample there. Writes:
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
from plotting import savefig, savefig_plotnine, set_mpl_rcparams
from multimodal_readwrite import apply_affine, mesmer_feature_table_to_adata, xenium_pixel_size
from xenium_readwrite import match_cells_spatial_mnn, read_xenium_sample


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


def align_proteomics_cells(feature_table, affine, xenium_dir):
    """Proteomics cells (feature table metadata) and their centroids in Xenium µm, as in
    read_xenium_proteomics_sample (analysis_repo). xenium_dir can be a 10x output folder or a
    Proseg raw_results folder (the pixel size then comes from the matching raw Xenium run)."""
    # Xenium pixel size from experiment.xenium
    pixel_size = xenium_pixel_size(xenium_dir)
    # proteomics centroids: proteomics px -> Xenium px -> Xenium µm
    cells = mesmer_feature_table_to_adata(feature_table).obs
    matrix = np.loadtxt(affine, delimiter=",")
    aligned_um = apply_affine(cells[["x_px", "y_px"]].to_numpy(), matrix) * pixel_size
    print(f"{len(cells)} proteomics cells, Xenium pixel size {pixel_size} µm/px")
    return cells, aligned_um


def read_xenium_cells(xenium_dir):
    """Xenium cells (obs, indexed by cell ID) with centroids in µm, read as in analysis_repo
    (10x output folder or Proseg raw_results folder)."""
    xenium = read_xenium_sample(xenium_dir, anndata=True, gex_only=True).obs
    print(f"{len(xenium)} Xenium cells ({xenium_dir})")
    return xenium


def match_cells(cells, aligned_um, xenium, max_dist):
    """Mutual nearest neighbours within max_dist µm, as in read_xenium_proteomics_sample."""
    # aligned proteomics centroids indexed like the proteomics cells
    prot_xy = pd.DataFrame(aligned_um, index=cells.index, columns=["x_centroid", "y_centroid"])
    pairs = match_cells_spatial_mnn(xenium, prot_xy, max_dist=max_dist)
    return pd.DataFrame({
        "proteomics_cell_id": cells.loc[pairs["xr_cell"], "cell_id"].to_numpy(),
        "xenium_cell_id": pairs["proseg_cell"].to_numpy(),
        "dist_um": pairs["dist"].to_numpy(),
    })


def main():
    args = parse_args()
    set_mpl_rcparams()

    # proteomics centroids in Xenium µm
    cells, aligned_um = align_proteomics_cells(args.feature_table, args.affine, args.xenium_dir)

    # Xenium centroids (µm)
    xenium = read_xenium_cells(args.xenium_dir)
    xenium_um = xenium[["x_centroid", "y_centroid"]].to_numpy()

    # mutual nearest neighbours within max_dist
    matches = match_cells(cells, aligned_um, xenium, args.max_dist)
    matches.to_csv(args.out_matches, index=False)
    print(f"{len(matches)} matched cells within {args.max_dist} µm")

    # overlay of the centroids (many points: matplotlib)
    plt.figure(figsize=(8, 8))
    plt.scatter(xenium_um[:, 0], xenium_um[:, 1], s=0.05, alpha=0.4, linewidths=0, label="Xenium", color="#19618A")
    plt.scatter(aligned_um[:, 0], aligned_um[:, 1], s=0.05, alpha=0.4, linewidths=0, label="proteomics (aligned)", color="#09B39C")
    plt.gca().set_aspect("equal")
    plt.gca().invert_yaxis()
    plt.title(f"Aligned centroids {args.sample_id}", fontsize=9)
    plt.xlabel("x (µm)")
    plt.ylabel("y (µm)")
    plt.legend(markerscale=30)
    savefig(None, args.out_overlay, dpi=200)

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
                text=element_text(font="Helvetica"),
                axis_text_x=element_text(angle=45, hjust=1, color="black"),
                legend_position="none",
                figure_size=(4, 4))
        # headroom above the tallest bar for its label
        + scale_y_continuous(expand=[0, 0, 0.01, 0], limits=(0, counts["ncells"].max() * 1.1))
        + scale_fill_manual(values=colors)
        + xlab("")
        + ylab("Number of cells")
    )
    savefig_plotnine(p, counts, args.out_counts)
    print(f"written: {args.out_matches}, {args.out_overlay}, {args.out_counts}")


if __name__ == "__main__":
    main()
