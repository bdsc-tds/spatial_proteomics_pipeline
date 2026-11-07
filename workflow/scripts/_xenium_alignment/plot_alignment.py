"""Bring the proteomics cells into Xenium coordinates, match them to Xenium cells and plot the result.

The proteomics centroids (pixels of the proteomics image) are mapped with the alignment matrix to
Xenium morphology pixels, then to µm with the Xenium pixel size from experiment.xenium. Cells are
matched as mutual nearest neighbours within --max_dist µm. Reading, alignment and matching use the
helpers of utils/multimodal_readwrite.py and utils/xenium_readwrite.py, so the matches are the same
as those of read_xenium_proteomics_sample. Writes:
  - the table of matched cells (proteomics cell_id, Xenium cell_id, distance in µm)
  - an overlay of Xenium and aligned proteomics centroids
  - per modality (Xenium, aligned proteomics), the matched and the unmatched cells side by side
  - a bar plot of the number of cells of each Xenium segmentation (--xenium_segs, then the default
    10x_mm_5um one), of proteomics, matched and unmatched proteomics cells
  - a stacked bar of the fraction of proteomics cells matched / unmatched to a Xenium cell
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import anndata as ad
import squidpy as sq
from plotnine import (aes, element_rect, element_text, geom_col, geom_text, ggplot, ggtitle, position_stack,
                      scale_fill_manual, scale_y_continuous, theme, theme_classic, xlab, ylab)

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
    parser.add_argument("--xenium_segs", nargs="*", default=[],
                        help="name=output dir of other Xenium segmentations (10x cells.parquet or Proseg cell-metadata.csv.gz), cell numbers only")
    parser.add_argument("--max_dist", type=float, default=5, help="max distance (µm) between matched cells")
    parser.add_argument("--sample_id", default="", help="used in the plot titles")
    parser.add_argument("--out_matches", required=True)
    parser.add_argument("--out_overlay", required=True)
    parser.add_argument("--out_matched_xenium", required=True)
    parser.add_argument("--out_matched_proteomics", required=True)
    parser.add_argument("--out_counts", required=True)
    parser.add_argument("--out_fraction_matched", required=True)
    return parser.parse_args()


def align_proteomics_cells(feature_table, affine, xenium_dir):
    """Proteomics cells (feature table metadata) and their centroids in Xenium µm, as in
    read_xenium_proteomics_sample (multimodal_readwrite). xenium_dir can be a 10x output folder or a
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
    """Xenium cells (obs, indexed by cell ID) with centroids in µm, read with read_xenium_sample
    (10x output folder or Proseg raw_results folder)."""
    xenium = read_xenium_sample(xenium_dir, anndata=True, gex_only=True).obs
    print(f"{len(xenium)} Xenium cells ({xenium_dir})")
    return xenium


def count_xenium_cells(seg_dir):
    """Number of cells of a Xenium segmentation output: 10x cells.parquet or Proseg cell-metadata.csv.gz."""
    # 10x / Xenium Ranger: row count from the parquet metadata
    cells_parquet = os.path.join(seg_dir, "cells.parquet")
    if os.path.isfile(cells_parquet):
        return pq.ParquetFile(cells_parquet).metadata.num_rows
    # Proseg: one row per cell
    return len(pd.read_csv(os.path.join(seg_dir, "cell-metadata.csv.gz"), usecols=[0]))


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


def plot_matched_spatial(xy, matched, name, sample_id, out_file):
    """Cells of one modality (centroids in Xenium µm): matched cells highlighted (left) and unmatched
    cells highlighted (right), the other cells in gray, with squidpy's spatial scatter so both panels
    keep the µm aspect ratio."""
    gray = "#D9D9D9"
    n_matched = int(matched.sum())
    n_unmatched = len(xy) - n_matched
    # figure height from the tissue's height/width ratio, so the equal-aspect panels fill the figure
    extent = np.ptp(xy, axis=0)
    panel_width = 6
    fig, axes = plt.subplots(1, 2, figsize=(2 * panel_width + 2, panel_width * extent[1] / extent[0] + 1),
                             layout="constrained")
    panels = [(axes[0], matched, "matched", "unmatched", "#2166AC",
               f"matched: {n_matched:,} / {len(xy):,} ({n_matched / len(xy):.1%})"),
              (axes[1], ~matched, "unmatched", "matched", "#E66101",
               f"unmatched: {n_unmatched:,} / {len(xy):,} ({n_unmatched / len(xy):.1%})")]
    for ax, highlight, label, other, color, title in panels:
        # minimal AnnData with the centroids (µm) in obsm["spatial"] and the match status in obs
        status = pd.Categorical(np.where(highlight, label, other), categories=[label, other])
        adata = ad.AnnData(obs=pd.DataFrame({"status": status}, index=[str(i) for i in range(len(xy))]))
        adata.obsm["spatial"] = xy
        # many points: squidpy scatter (no image), gray cells drawn first, highlighted cells on top
        order = np.argsort(highlight, kind="stable")
        sq.pl.spatial_scatter(adata[order], shape=None, color="status", img=False, size=0.5,
                              palette=ListedColormap([color, gray]), ax=ax, frameon=True,
                              legend_loc="right margin")
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("x (µm)")
        ax.set_ylabel("y (µm)")
    fig.suptitle(f"{name} cells {sample_id}", fontsize=10)
    savefig(None, out_file, dpi=200)


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

    # matched / unmatched cells, one figure per modality
    xenium_matched = xenium.index.isin(matches["xenium_cell_id"])
    prot_matched = cells["cell_id"].isin(matches["proteomics_cell_id"]).to_numpy()
    plot_matched_spatial(xenium_um, xenium_matched, "Xenium", args.sample_id, args.out_matched_xenium)
    plot_matched_spatial(aligned_um, prot_matched, "Proteomics (aligned)", args.sample_id,
                         args.out_matched_proteomics)

    # cells of the other Xenium segmentations, name=dir
    other_segs = dict(s.split("=", 1) for s in args.xenium_segs)
    other_ncells = {seg: count_xenium_cells(d) for seg, d in sorted(other_segs.items())}

    # counts of every Xenium segmentation, proteomics, matched and unmatched proteomics cells (matches
    # to the default 10x_mm_5um segmentation)
    counts = pd.DataFrame({
        "modality": [*other_ncells, "10x_mm_5um", "proteomics", "matched", "unmatched"],
        "ncells": [*other_ncells.values(), len(xenium), len(cells), len(matches), len(cells) - len(matches)],
    })
    counts["modality"] = pd.Categorical(counts["modality"], categories=counts["modality"], ordered=True)
    print(counts.to_string(index=False))
    colors = {"10x_mm_5um": "#19618A", "proteomics": "#09B39C", "matched": "#2166AC", "unmatched": "#E66101"}
    # other Xenium segmentations: shades of purple
    other_colors = ["#7B1FA2", "#BA68C8", "#4A148C", "#CE93D8"]
    colors.update({seg: other_colors[i % len(other_colors)] for i, seg in enumerate(other_ncells)})
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

    # fraction of proteomics cells matched / unmatched to a Xenium cell, one stacked bar
    fractions = pd.DataFrame({
        "cells": "proteomics",
        "status": pd.Categorical(["unmatched", "matched"], categories=["unmatched", "matched"]),
        "ncells": [len(cells) - len(matches), len(matches)],
    })
    fractions["fraction"] = fractions["ncells"] / len(cells)
    fractions["label"] = [f"{f:.1%}\n({n:,})" for f, n in zip(fractions["fraction"], fractions["ncells"])]
    p = (
        ggplot(fractions, aes(x="cells", y="fraction", fill="status"))
        + geom_col(position="stack", width=0.8)
        + geom_text(aes(label="label"), position=position_stack(vjust=0.5), size=7, color="white")
        + theme_classic()
        + theme(panel_border=element_rect(fill=None, color="black"),
                text=element_text(font="Helvetica"),
                axis_text_x=element_text(color="black"),
                figure_size=(3, 4))
        + scale_y_continuous(expand=[0, 0, 0.01, 0], limits=(0, 1), labels=lambda v: [f"{x:.0%}" for x in v])
        + scale_fill_manual(values={"matched": "#2166AC", "unmatched": "#E66101"})
        + ggtitle(args.sample_id)
        + xlab("")
        + ylab("Fraction of cells matched to Xenium")
    )
    savefig_plotnine(p, fractions, args.out_fraction_matched)
    print(f"written: {args.out_matches}, {args.out_overlay}, {args.out_counts}, {args.out_fraction_matched}")


if __name__ == "__main__":
    main()
