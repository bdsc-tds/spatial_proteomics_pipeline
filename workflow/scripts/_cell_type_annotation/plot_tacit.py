"""Plot the TACIT cell types of one segmentation method over all samples.

Inputs, one per sample in --sample_ids order (tissue/panel/donor/sample): the quantify_cells feature
table, the run_tacit.R labels and kNN fractions, and the signature CSV of the sample's panel.

Cell categories: unique (not mixed), mixed (in a mixed microcluster, labelled by TACIT's kNN),
unassigned (Others) and stringent pure (positive for its final cell type only with breakpoint 4 as
cell type score threshold).

Summary figures (folder of --out_fractions):
  - cell type and unique / mixed / unassigned fractions per sample, one facet per tissue
  - fraction of unique and of stringent pure cells per cell type across samples
  - fraction of mixed cells per pair of positive cell types, with the final cell types of each pair
  - signature weights, one heatmap per panel
Per sample figures (folder of its --out_sample_fractions):
  - the same bar plots for the sample alone, plus the mixed cells positive for 3+ cell types
  - mean marker z-score per cell type (as given to TACIT), all cells and unique cells only
  - kNN neighbour fractions of the mixed cells
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from plotnine import (aes, coord_flip, element_blank, element_rect, element_text, facet_grid, facet_wrap,
                      geom_boxplot, geom_col, geom_jitter, geom_text, geom_violin, ggplot, guide_legend, guides, labs,
                      scale_color_manual, scale_fill_manual, scale_x_discrete, scale_y_continuous, scale_y_log10,
                      theme, theme_classic)

# PyComplexHeatmap still calls the removed matplotlib.cm.get_cmap in some versions
import matplotlib.cm
if not hasattr(matplotlib.cm, "get_cmap"):
    matplotlib.cm.get_cmap = matplotlib.colormaps.get_cmap
import PyComplexHeatmap as pch

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(script_dir, "../utils/"))
from plotting import TACIT_mm_cols, panel_ratios, tissue_cols, savefig, savefig_plotnine, set_mpl_rcparams

OTHERS = "Others"
OTHERS_COLOR = "#BDBDBD"
ASSIGNMENT_COLORS = {"unique": "#4A6FA5", "mixed": "#C08A4E", "unassigned": OTHERS_COLOR}
# mixed cells positive for three or more cell types, one bar / box in the pairs plots
MULTI = "3+ cell types"
PAIR_COLORS = {"2 cell types": "#BFBFBF", MULTI: "#5F5F5F"}
# heatmap z-scores clipped to +-Z_CLIP, as Z_CLIP in run_tacit.R
Z_CLIP = 10
# mean z-score per cell type and marker clipped to +-HEATMAP_CLIP for the heatmap colours
HEATMAP_CLIP = 2
# most frequent cell type pairs shown in the pairs plots
N_PAIRS = 20


def sample_label(sample_id, sep="/"):
    """Short sample name for the plots: the region of a Xenium output folder
    (output-{instrument}__{slide}__{region}__{date}__{time}), otherwise donor/sample (sep="\\n" for two rows)."""
    parts = sample_id.split("/")[2:]
    xenium = parts[-1].split("__")
    if parts[-1].startswith("output-") and len(xenium) == 5:
        return xenium[2]
    return sep.join(parts)


def empty_figure(text, title, out):
    """Placeholder figure with a centered message."""
    plt.figure(figsize=(6, 2))
    plt.text(0.5, 0.5, text, ha="center", va="center", fontsize=13)
    plt.axis("off")
    plt.title(title)
    savefig(None, out)


def stacked_by_tissue(counts, fill, values, breaks, fill_title, out):
    """Stacked bar of the fraction of cells per sample, one facet per tissue."""
    n_samples = counts["sample_id"].nunique()
    p = (
        ggplot(counts, aes("sample_id", "fraction", fill=fill))
        + geom_col(width=0.8)
        + facet_grid(cols="tissue", scales="free_x", space="free_x")
        + scale_x_discrete(labels=lambda ids: [sample_label(s, sep="\n") for s in ids])
        + scale_fill_manual(values=values, breaks=breaks)
        + labs(x="", y="fraction of cells", fill=fill_title)
        + theme_classic()
        + theme(panel_border=element_rect(fill=None, color="black"),
                text=element_text(font="Helvetica"),
                axis_text_x=element_text(rotation=90, ha="center", va="top"),
                figure_size=(0.3 * n_samples + 4, 5))
        + scale_y_continuous(expand=[0, 0, 0.01, 0])
    )
    savefig_plotnine(p, None, out)


def pure_per_cell_type(cells, ct_order):
    """Number and fraction of pure cells per sample and final cell type (Others left out)."""
    assigned = cells[cells["cell_type"] != OTHERS]
    pure = (assigned.groupby(["tissue", "sample_id", "cell_type"])["pure"].agg(["sum", "mean"])
            .rename(columns={"sum": "pure_cells", "mean": "fraction"}).reset_index())
    present = [ct for ct in ct_order if ct in set(pure["cell_type"])]
    pure["cell_type"] = pd.Categorical(pure["cell_type"], categories=present)
    return pure, present


def pairs_top(mixed, col, order, ct_order, palette, title):
    """Stacked bar of the final cell types of the mixed cells of each value of col (pairs or combinations,
    in order), x axis hidden: the top of a composition with a plot of the same pairs below."""
    # fraction of the final cell types per pair, unit scaled, cell types in signature order
    labels = pd.crosstab(mixed[col], mixed["cell_type"], normalize="index")
    labels = labels.reindex(index=order, columns=[ct for ct in ct_order if ct in labels.columns], fill_value=0)
    long = labels.rename_axis(index="pair", columns="cell_type").stack().rename("fraction").reset_index()
    long["pair"] = pd.Categorical(long["pair"], categories=order)
    long["cell_type"] = pd.Categorical(long["cell_type"], categories=list(labels.columns))

    # first signature cell type on top of each bar
    p = (
        ggplot(long, aes("pair", "fraction", fill="cell_type"))
        + geom_col(width=0.8)
        + scale_fill_manual(values=palette, breaks=list(labels.columns))
        + labs(x="", y="final cell type\nfraction", fill="final TACIT cell type", title=title)
        + theme_classic()
        + theme(panel_border=element_rect(fill=None, color="black"),
                text=element_text(font="Helvetica"),
                axis_text_x=element_blank(), axis_ticks_major_x=element_blank(),
                figure_size=(0.32 * len(order) + 5, 9))
        + scale_y_continuous(expand=[0, 0, 0.01, 0])
    )
    return p


def sample_bars(data, x, order, fill_values, y_title, title, y_scale=None):
    """Horizontal bar per category of one sample, categories top to bottom in order."""
    data = data.assign(**{x: pd.Categorical(data[x], categories=order[::-1])})
    p = (
        ggplot(data, aes(x, "value", fill=x))
        + geom_col(width=0.8, show_legend=False)
        + coord_flip()
        + scale_fill_manual(values=fill_values)
        + labs(x="", y=y_title, title=title)
        + theme_classic()
        + theme(panel_border=element_rect(fill=None, color="black"),
                text=element_text(font="Helvetica"),
                figure_size=(6, 0.25 * len(order) + 1.5))
        + (y_scale if y_scale is not None else scale_y_continuous(expand=[0, 0, 0.01, 0]))
    )
    return p


def sample_mixed_counts(cells, col, ct_order, palette, title, empty_text, out):
    """Mixed cells of one sample per value of col (pairs or combinations, N_PAIRS most frequent, MULTI
    last), log scale, below the final TACIT cell types of each."""
    counts = cells[col].value_counts()
    if counts.empty:
        empty_figure(empty_text, title, out)
        return
    order = [p for p in counts.index if p != MULTI][:N_PAIRS] + ([MULTI] if MULTI in counts.index else [])
    counts = counts.reindex(order).rename("cells").rename_axis("pair").reset_index()
    counts["pair"] = pd.Categorical(counts["pair"], categories=order)
    # dark for 3+ cell types
    counts["n_types"] = ["2 cell types" if p != MULTI and p.count(" + ") == 1 else MULTI for p in order]

    # mixed cells per pair / combination (log scale), count above each bar
    p_bot = (
        ggplot(counts, aes("pair", "cells"))
        + geom_col(aes(fill="n_types"), width=0.8, show_legend=False)
        + geom_text(aes(label="cells"), angle=90, va="bottom", nudge_y=0.05, size=9, family="Helvetica")
        + scale_fill_manual(values=PAIR_COLORS)
        + labs(x="", y="mixed cells (log scale)")
        + theme_classic()
        + theme(panel_border=element_rect(fill=None, color="black"),
                text=element_text(font="Helvetica"),
                axis_text_x=element_text(rotation=45, ha="right"),
                figure_size=(0.32 * len(order) + 5, 9))
        + scale_y_log10(limits=(0.8, counts["cells"].max() * 5), expand=[0, 0, 0, 0])
    )
    savefig_plotnine(panel_ratios(pairs_top(cells, col, order, ct_order, palette, title) / p_bot, [1, 2]), None, out)


def zero_white_cmap(vmin):
    """RdBu_r stretched so that vmin (< 0) is its darkest blue, 0 white and 1 its darkest red (white to red
    when vmin is 0)."""
    base = matplotlib.colormaps["RdBu_r"]
    reds = base(np.linspace(0.5, 1, 128))
    # RdBu_r's midpoint is light grey: pure white at 0
    reds[0] = (1, 1, 1, 1)
    if vmin >= 0:
        return matplotlib.colors.LinearSegmentedColormap.from_list("white_red", reds)
    # position of 0 between vmin and 1, blues below it and reds above it
    zero = -vmin / (1 - vmin)
    blues = base(np.linspace(0, 0.5, 128))
    blues[-1] = (1, 1, 1, 1)
    stops = list(zip(np.linspace(0, zero, 128), blues)) + list(zip(np.linspace(zero, 1, 128)[1:], reds[1:]))
    return matplotlib.colors.LinearSegmentedColormap.from_list("blue_white_red", stops)


def plot_heatmap(features, cell_types, signature, cofactor, ct_order, palette, title, out, keep=None):
    """Mean z-score of every signature marker per cell type, markers split by whether they are used in the
    signature. With keep (boolean mask), only those cells are averaged (z-scores still over all cells)."""
    markers = list(signature.columns)
    # arcsinh, then z-score per marker over all cells clipped to +-Z_CLIP, as run_tacit.R (a constant marker is 0)
    norm = np.arcsinh(features[markers] / cofactor)
    norm = ((norm - norm.mean()) / norm.std()).clip(-Z_CLIP, Z_CLIP).replace([np.inf, -np.inf], np.nan).fillna(0)
    if keep is not None:
        keep = np.asarray(keep, dtype=bool)
        norm, cell_types = norm[keep], cell_types[keep]
    # mean z-score per cell type and marker, rows in signature order
    means = norm.groupby(cell_types.to_numpy()).mean()
    means = means.loc[[ct for ct in ct_order if ct in means.index]]
    # mean z-scores clipped to +-HEATMAP_CLIP, colour scale over that range (white at 0)
    means = means.clip(-HEATMAP_CLIP, HEATMAP_CLIP)
    vmax = HEATMAP_CLIP

    # number of cells per cell type as a row annotation
    n_cells = cell_types.value_counts().loc[means.index].to_frame("N cells")
    right = pch.HeatmapAnnotation(**{"N cells": pch.anno_barplot(n_cells, colors="grey", legend=False)}, axis=0,
                                  verbose=0)
    # used / not used in the signature as a column annotation, columns split by it when both occur
    used = (signature.fillna(0) != 0).any(axis=0)
    usage = pd.Series(np.where(used, "yes", "no"), index=markers)
    top = pch.HeatmapAnnotation(df=usage.to_frame("Used for phenotyping"),
                                colors={"Used for phenotyping": {"yes": "#4D4D4D", "no": "#D9D9D9"}},
                                label_kws={"visible": False}, axis=1, verbose=0)
    split = usage if usage.nunique() > 1 else None
    # cell type colour as a row annotation
    ct = pd.DataFrame({"cell type": means.index}, index=means.index)
    left = pch.HeatmapAnnotation(df=ct, colors={"cell type": {c: palette[c] for c in means.index}},
                                 label_kws={"visible": False}, plot_legend=False, axis=0, verbose=0)

    plt.figure(figsize=(0.2 * len(markers) + 2.5, 0.16 * len(means) + 1.6))
    cm = pch.ClusterMapPlotter(data=means, row_cluster=False, col_cluster=True, col_dendrogram=False,
                               col_split=split, col_split_order=["yes", "no"], col_split_gap=1.5,
                               top_annotation=top, left_annotation=left, right_annotation=right,
                               show_rownames=True, show_colnames=True,
                               row_names_side="left", cmap=zero_white_cmap(-1), vmin=-vmax, vmax=vmax,
                               label="mean z-score", legend_kws={"extend": "both"}, legend_hpad=3,
                               xticklabels_kws={"labelrotation": 45},
                               xlabel="Protein", ylabel="TACIT automated cell types", ylabel_kws={"labelpad": 110},
                               verbose=0)
    # PyComplexHeatmap hides the row names with left and right annotations: write them next to the colour strip
    strip = cm.left_annotation.axes[0, 0]
    strip.yaxis.set_visible(True)
    strip.set_yticks(cm.heatmap_axes[0, 0].get_yticks(), labels=means.index)
    strip.tick_params(axis="y", left=False, labelleft=True, labelsize=plt.rcParams["ytick.labelsize"])
    # colorbar title centered above the bar, 0.8 inch right of the cell number bars
    cbar = next(c for c in cm.cbars if isinstance(c, matplotlib.colorbar.Colorbar))
    cbar.set_label("")
    cbar.ax.set_title("mean\nz-score", loc="center", fontsize=13, pad=8)
    bars = cm.right_annotation.axes[0, -1].get_position()
    pos = cbar.ax.get_position()
    cbar.ax.set_position([bars.x1 + 0.8 / plt.gcf().get_figwidth(), pos.y0, pos.width, pos.height])
    # used for phenotyping legend below the colorbar, no frame
    for legend in cm.cbars:
        if isinstance(legend, matplotlib.legend.Legend):
            legend.set_frame_on(False)
            legend.set_bbox_to_anchor((cbar.ax.get_position().x0 - 0.3 / plt.gcf().get_figwidth(),
                                       cbar.ax.get_position().y0 - 0.15), transform=plt.gcf().transFigure)
            legend.set_loc("upper left")
    plt.suptitle(title, y=1.02)
    savefig(None, out)


# ---------------------------------------------------------------- arguments
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--sample_ids", nargs="+", required=True, help="tissue/panel/donor/sample of every sample")
parser.add_argument("--feature_tables", nargs="+", required=True, help="quantify_cells feature tables")
parser.add_argument("--cell_types", nargs="+", required=True, help="run_tacit.R outputs")
parser.add_argument("--knn_fractions", nargs="+", required=True, help="run_tacit.R kNN fractions of the mixed cells")
parser.add_argument("--signatures", nargs="+", required=True, help="signature CSV of each sample's panel")
parser.add_argument("--cofactor", type=float, required=True, help="arcsinh cofactor")
parser.add_argument("--threads", type=int, default=1, help="unused")
parser.add_argument("--out_fractions", required=True,
                    help="stacked bar plot of the cell type fractions; the other summary figures go to its folder")
parser.add_argument("--out_sample_fractions", nargs="+", required=True,
                    help="per sample, cell type fractions; the sample's other figures go to its folder")
args = parser.parse_args()

# every per-sample list must follow --sample_ids
for name in ["feature_tables", "cell_types", "knn_fractions", "signatures", "out_sample_fractions"]:
    if len(getattr(args, name)) != len(args.sample_ids):
        raise ValueError(f"--{name}: {len(getattr(args, name))} values for {len(args.sample_ids)} samples")

set_mpl_rcparams()

# ---------------------------------------------------------------- inputs
signatures = [pd.read_csv(path).set_index("cell_type") for path in args.signatures]
labels = {s: pd.read_csv(path) for s, path in zip(args.sample_ids, args.cell_types)}

# cell types of all signatures in order of appearance, then Others
ct_order = list(dict.fromkeys(ct for sig in signatures for ct in sig.index)) + [OTHERS]
# fixed colours for the known cell types, tab20 / husl for any other
unknown = [ct for ct in ct_order[:-1] if ct not in TACIT_mm_cols]
colors = sns.color_palette("tab20" if len(unknown) <= 20 else "husl", len(unknown)).as_hex()
palette = {ct: TACIT_mm_cols[ct] for ct in ct_order if ct in TACIT_mm_cols} | dict(zip(unknown, colors))
palette[OTHERS] = OTHERS_COLOR

# labels of all samples in one table, with the tissue and the assignment (unique / mixed / unassigned)
df = pd.concat([lab.assign(sample_id=s, tissue=s.split("/")[0]) for s, lab in labels.items()], ignore_index=True)
df["positive_cell_types"] = df["positive_cell_types"].fillna("")
df["assignment"] = np.where(df["cell_type"] == OTHERS, "unassigned", np.where(df["mixed"], "mixed", "unique"))

# mixed cells: combination of positive cell types (signature order) and pair (MULTI for 3+)
rank = {ct: i for i, ct in enumerate(ct_order)}
positives = df["positive_cell_types"].str.split("::")
df["combination"] = [" + ".join(sorted(p, key=lambda ct: rank.get(ct, len(rank)))) if m else None
                     for p, m in zip(positives, df["mixed"].astype(bool))]
df["pair"] = [(c if len(p) == 2 else MULTI) if c is not None else None for c, p in zip(df["combination"], positives)]

# stringent pure cells: positive for their final cell type only with breakpoint 4 as cell type score threshold
df["pure"] = df["positive_cell_types_group4"].fillna("").eq(df["cell_type"]) & (df["cell_type"] != OTHERS)

# samples and tissues in input order; fixed colours for the known tissues, Dark2 for any other
sample_order = list(dict.fromkeys(df["sample_id"]))
tissue_order = list(dict.fromkeys(df["tissue"]))
unknown = [t for t in tissue_order if t not in tissue_cols]
tissues = ({t: tissue_cols[t] for t in tissue_order if t in tissue_cols}
           | dict(zip(unknown, sns.color_palette("Dark2", len(unknown)).as_hex())))

summary_dir = os.path.dirname(args.out_fractions)

# ---------------------------------------------------------------- summary: cell type fractions per sample
counts = df.groupby(["tissue", "sample_id", "cell_type"]).size().rename("cells").reset_index()
counts["fraction"] = counts["cells"] / counts.groupby("sample_id")["cells"].transform("sum")
present = [ct for ct in ct_order if ct in set(counts["cell_type"])]
counts["sample_id"] = pd.Categorical(counts["sample_id"], categories=sample_order)
counts["tissue"] = pd.Categorical(counts["tissue"], categories=tissue_order)
counts["cell_type"] = pd.Categorical(counts["cell_type"], categories=present)
stacked_by_tissue(counts, "cell_type", palette, present, "TACIT cell type", args.out_fractions)

# ---------------------------------------------------------------- summary: unique / mixed / unassigned per sample
counts = df.groupby(["tissue", "sample_id", "assignment"]).size().rename("cells").reset_index()
counts["fraction"] = counts["cells"] / counts.groupby("sample_id")["cells"].transform("sum")
counts["sample_id"] = pd.Categorical(counts["sample_id"], categories=sample_order)
counts["tissue"] = pd.Categorical(counts["tissue"], categories=tissue_order)
counts["assignment"] = pd.Categorical(counts["assignment"], categories=list(ASSIGNMENT_COLORS))
stacked_by_tissue(counts, "assignment", ASSIGNMENT_COLORS, list(ASSIGNMENT_COLORS), "TACIT assignment",
                  os.path.join(summary_dir, "mixed_fraction.png"))

# ---------------------------------------------------------------- summary: unique fraction per cell type
# fraction of the cells of each final cell type that are not mixed, one point per sample (Others are never mixed)
assigned = df[df["cell_type"] != OTHERS]
frac = (assigned.assign(unique=~assigned["mixed"].astype(bool))
        .groupby(["tissue", "sample_id", "cell_type"])["unique"].mean().rename("fraction").reset_index())
present = [ct for ct in ct_order if ct in set(frac["cell_type"])]
frac["cell_type"] = pd.Categorical(frac["cell_type"], categories=present)
p = (
    ggplot(frac, aes("cell_type", "fraction"))
    + geom_boxplot(aes(fill="cell_type"), outlier_shape="", width=0.7)
    + geom_jitter(aes(color="tissue"), width=0.2, height=0, size=1.2, alpha=0.8, random_state=0)
    + scale_fill_manual(values=palette)
    + scale_color_manual(values=tissues, breaks=list(tissues))
    + guides(fill=guide_legend(order=1), color=guide_legend(order=2, override_aes={"size": 3, "alpha": 1}))
    + labs(x="", y="fraction of unique cells (not mixed)", fill="TACIT cell type", color="tissue")
    + theme_classic()
    + theme(panel_border=element_rect(fill=None, color="black"),
            text=element_text(font="Helvetica"),
            axis_text_x=element_text(rotation=45, ha="right"),
            figure_size=(0.45 * len(present) + 4, 5))
    + scale_y_continuous(expand=[0, 0, 0.01, 0], limits=(0, 1))
)
savefig_plotnine(p, None, os.path.join(summary_dir, "unique_fraction_per_cell_type.png"))

# ---------------------------------------------------------------- summary: pure fraction per cell type
# fraction of pure cells per cell type across samples, below the number of pure cells (all samples, log scale)
pure, present = pure_per_cell_type(df, ct_order)
size = (0.45 * len(present) + 4, 6)

# pure cells per cell type over all samples (log scale, empty bars at 0.8), count above each bar, x axis hidden
totals = pure.groupby("cell_type", observed=False)["pure_cells"].sum().rename("cells").reset_index()
totals["bar"] = totals["cells"].clip(lower=0.8)
p_top = (
    ggplot(totals, aes("cell_type", "bar", fill="cell_type"))
    + geom_col(width=0.7, show_legend=False)
    + geom_text(aes(label="cells"), angle=90, va="bottom", nudge_y=0.05, size=9, family="Helvetica")
    + scale_fill_manual(values=palette)
    + labs(x="", y="pure cells\n(log scale)")
    + theme_classic()
    + theme(panel_border=element_rect(fill=None, color="black"),
            text=element_text(font="Helvetica"),
            axis_text_x=element_blank(), axis_ticks_major_x=element_blank(),
            figure_size=size)
    + scale_y_log10(limits=(0.8, max(totals["cells"].max(), 1) * 300), expand=[0, 0, 0, 0])
)

# one box per cell type, one jittered point per sample coloured by tissue
p_bot = (
    ggplot(pure, aes("cell_type", "fraction"))
    + geom_boxplot(aes(fill="cell_type"), outlier_shape="", width=0.7, show_legend=False)
    + geom_jitter(aes(color="tissue"), width=0.2, height=0, size=1.2, alpha=0.8, random_state=0)
    + scale_fill_manual(values=palette)
    + scale_color_manual(values=tissues, breaks=list(tissues))
    + guides(color=guide_legend(override_aes={"size": 3, "alpha": 1}))
    + labs(x="", y="fraction of pure cells\n(breakpoint 4 as cell type score threshold)", color="tissue")
    + theme_classic()
    + theme(panel_border=element_rect(fill=None, color="black"),
            text=element_text(font="Helvetica"),
            axis_text_x=element_text(rotation=45, ha="right"),
            figure_size=size)
    + scale_y_continuous(expand=[0, 0, 0.01, 0], limits=(0, 1))
)
savefig_plotnine(panel_ratios(p_top / p_bot, [1, 3]), None,
                 os.path.join(summary_dir, "stringent_pure_fraction_per_cell_type.png"))

# ---------------------------------------------------------------- summary: mixed cells per cell type pair
# fraction of each sample's mixed cells per pair (N_PAIRS most frequent over all samples, 3+ last), log scale
mixed = df[df["pair"].notna()]
totals = mixed["pair"].value_counts()
order = [p for p in totals.index if p != MULTI][:N_PAIRS] + ([MULTI] if MULTI in totals.index else [])
counts = pd.crosstab(mixed["sample_id"], mixed["pair"])
frac = counts.div(counts.sum(axis=1), axis=0).reindex(columns=order, fill_value=0)
frac = frac.stack().rename("fraction").reset_index()
# samples without a pair left out (log scale)
frac = frac[frac["fraction"] > 0]
frac["tissue"] = frac["sample_id"].str.split("/").str[0]
frac["pair"] = pd.Categorical(frac["pair"], categories=order)
frac["n_types"] = np.where(frac["pair"] == MULTI, MULTI, "2 cell types")

# one box per pair filled by 2 / 3+ cell types, one jittered point per sample coloured by tissue
p_bot = (
    ggplot(frac, aes("pair", "fraction"))
    + geom_boxplot(aes(fill="n_types"), outlier_shape="", width=0.7, show_legend=False)
    + geom_jitter(aes(color="tissue"), width=0.2, height=0, size=1.2, alpha=0.8, random_state=0)
    + scale_x_discrete(limits=order)
    + scale_fill_manual(values=PAIR_COLORS)
    + scale_color_manual(values=tissues, breaks=list(tissues))
    + guides(color=guide_legend(override_aes={"size": 3, "alpha": 1}))
    + labs(x="", y="fraction of mixed cells\n(log scale)", color="tissue")
    + theme_classic()
    + theme(panel_border=element_rect(fill=None, color="black"),
            text=element_text(font="Helvetica"),
            axis_text_x=element_text(rotation=45, ha="right"),
            figure_size=(0.32 * len(order) + 5, 9))
    + scale_y_log10()
)
savefig_plotnine(panel_ratios(pairs_top(mixed, "pair", order, ct_order, palette, "") / p_bot, [1, 2]), None,
                 os.path.join(summary_dir, "mixed_cell_type_pairs.png"))

# ---------------------------------------------------------------- summary: signature weights per panel
for path in dict.fromkeys(args.signatures):
    panel = os.path.splitext(os.path.basename(path))[0]
    title = f"{panel}: TACIT signature"
    out = os.path.join(summary_dir, f"signature_{panel}.png")
    weights = pd.read_csv(path).set_index("cell_type").fillna(0).astype(float)
    if weights.empty:
        empty_figure("no cell types in the signature", title, out)
        continue

    # used markers ordered by the first cell type they define (block diagonal), then the unused ones
    used = (weights != 0).any(axis=0)
    first = weights.loc[:, used].ne(0).to_numpy().argmax(axis=0)
    order = list(weights.columns[used][np.argsort(first, kind="stable")]) + list(weights.columns[~used])
    weights = weights[order]

    # used for phenotyping as a column annotation (columns split by it), cell type colour as a row annotation
    usage = pd.Series(np.where(used[order], "yes", "no"), index=order)
    top = pch.HeatmapAnnotation(df=usage.to_frame("Used for phenotyping"),
                                colors={"Used for phenotyping": {"yes": "#4D4D4D", "no": "#D9D9D9"}},
                                label_kws={"visible": False}, axis=1, verbose=0)
    split = usage if usage.nunique() > 1 else None
    ct = pd.DataFrame({"cell type": weights.index}, index=weights.index)
    left = pch.HeatmapAnnotation(df=ct, colors={"cell type": {c: palette[c] for c in weights.index}},
                                 label_kws={"visible": False}, plot_legend=False, axis=0, verbose=0)

    # non-zero weights written in their cells; 0 white, 1 red, the lowest negative weight blue
    annot = weights.map(lambda w: f"{w:g}" if w != 0 else "")
    vmin, vmax = min(weights.min().min(), 0), max(weights.max().max(), 1)

    plt.figure(figsize=(0.2 * len(order) + 2.5, 0.2 * len(weights) + 1.6))
    cm = pch.ClusterMapPlotter(data=weights, row_cluster=False, col_cluster=False,
                               col_split=split, col_split_order=["yes", "no"], col_split_gap=1.5,
                               top_annotation=top, left_annotation=left, show_rownames=True, show_colnames=True,
                               row_names_side="left", cmap=zero_white_cmap(vmin / vmax), vmin=vmin, vmax=vmax,
                               annot=annot, fmt=None,
                               annot_kws={"fontsize": 6}, label="weight", legend_kws={"extend": "neither"},
                               xticklabels_kws={"labelrotation": 45}, xlabel="Protein", ylabel="cell type",
                               ylabel_kws={"labelpad": 110}, linewidths=0.3, linecolor="#E0E0E0", verbose=0)
    # cell type names next to the colour strip instead of right of the last block
    for heatmap_ax in cm.heatmap_axes.flat:
        heatmap_ax.tick_params(axis="y", right=False, labelright=False)
    strip = cm.left_annotation.axes[0, 0]
    strip.yaxis.set_visible(True)
    strip.set_yticks(cm.heatmap_axes[0, 0].get_yticks(), labels=weights.index)
    strip.tick_params(axis="y", left=False, labelleft=True, labelsize=plt.rcParams["ytick.labelsize"])
    # colorbar title centered above the bar, ticks at the lowest weight / 0 / highest, 0.6 inch right of the heatmap
    cbar = next(c for c in cm.cbars if isinstance(c, matplotlib.colorbar.Colorbar))
    cbar.set_label("")
    cbar.ax.set_title("weight", loc="center", fontsize=13, pad=8)
    cbar.set_ticks([vmin, 0, vmax] if vmin < 0 else [0, vmax / 2, vmax])
    right = cm.heatmap_axes[0, -1].get_position().x1
    cbar_pos = cbar.ax.get_position()
    cbar.ax.set_position([right + 0.6 / plt.gcf().get_figwidth(), cbar_pos.y0, cbar_pos.width, cbar_pos.height])
    # used for phenotyping legend below the colorbar, no frame
    for legend in cm.cbars:
        if isinstance(legend, matplotlib.legend.Legend):
            legend.set_frame_on(False)
            legend.set_bbox_to_anchor((cbar.ax.get_position().x0 - 0.3 / plt.gcf().get_figwidth(),
                                       cbar.ax.get_position().y0 - 0.1), transform=plt.gcf().transFigure)
            legend.set_loc("upper left")
    plt.suptitle(title, y=1.02)
    savefig(None, out)

# ---------------------------------------------------------------- per sample
for i, sample_id in enumerate(args.sample_ids):
    print(f"{sample_id} ({i + 1}/{len(args.sample_ids)})", flush=True)
    sample_dir = os.path.dirname(args.out_sample_fractions[i])
    title = sample_label(sample_id)
    sub = df[df["sample_id"] == sample_id]
    unit_scale = scale_y_continuous(expand=[0, 0, 0.01, 0], limits=(0, 1))

    # fraction of each cell type
    frac = sub["cell_type"].value_counts(normalize=True).rename("value").rename_axis("cell_type").reset_index()
    present = [ct for ct in ct_order if ct in set(frac["cell_type"])]
    savefig_plotnine(sample_bars(frac, "cell_type", present, palette, "fraction of cells", title), None,
                     args.out_sample_fractions[i])

    # fraction of unique, mixed and unassigned cells
    frac = sub["assignment"].value_counts(normalize=True).rename("value").rename_axis("assignment").reset_index()
    order = [a for a in ASSIGNMENT_COLORS if a in set(frac["assignment"])]
    p = sample_bars(frac, "assignment", order, ASSIGNMENT_COLORS, "fraction of cells", title, y_scale=unit_scale)
    savefig_plotnine(p, None, os.path.join(sample_dir, "mixed_fraction.png"))

    # unique fraction per final cell type (Others are never mixed)
    assigned = sub[sub["cell_type"] != OTHERS]
    frac = (~assigned["mixed"].astype(bool)).groupby(assigned["cell_type"]).mean().rename("value").reset_index()
    present = [ct for ct in ct_order if ct in set(frac["cell_type"])]
    p = sample_bars(frac, "cell_type", present, palette, "fraction of unique cells (not mixed)", title,
                    y_scale=unit_scale)
    savefig_plotnine(p, None, os.path.join(sample_dir, "unique_fraction_per_cell_type.png"))

    # mixed cells per pair (3+ cell types in one bar) and per combination of 3+ cell types
    sample_mixed_counts(sub[sub["pair"].notna()], "pair", ct_order, palette, title, "no mixed cells",
                        os.path.join(sample_dir, "mixed_cell_type_pairs.png"))
    sample_mixed_counts(sub[sub["pair"] == MULTI], "combination", ct_order, palette, title,
                        "no mixed cells positive for 3+ cell types", os.path.join(sample_dir, "mixed_cell_type_3plus.png"))

    # fraction of pure cells per cell type, next to the number of pure cells (log scale), cell types top to bottom
    pure, present = pure_per_cell_type(sub, ct_order)
    p_frac = sample_bars(pure.rename(columns={"fraction": "value"}), "cell_type", present, palette,
                         "fraction of pure cells\n(breakpoint 4 as cell type score threshold)", title,
                         y_scale=unit_scale)
    # number of pure cells (log scale, empty bars at 0.8), count right of each bar, cell type labels only on the left
    counts = pure.assign(cell_type=pd.Categorical(pure["cell_type"], categories=present[::-1]),
                         bar=pure["pure_cells"].clip(lower=0.8))
    # every second decade as break, the narrow panel has no room for all of them
    upper = max(counts["pure_cells"].max(), 1) * 30
    breaks = [10 ** k for k in range(0, int(np.log10(upper)) + 1, 2)]
    p_n = (
        ggplot(counts, aes("cell_type", "bar", fill="cell_type"))
        + geom_col(width=0.8, show_legend=False)
        + geom_text(aes(label="pure_cells"), ha="left", nudge_y=0.05, size=8, family="Helvetica")
        + coord_flip()
        + scale_fill_manual(values=palette)
        + labs(x="", y="pure cells\n(log scale)")
        + theme_classic()
        + theme(panel_border=element_rect(fill=None, color="black"),
                text=element_text(font="Helvetica"),
                axis_text_y=element_blank(), axis_ticks_major_y=element_blank(),
                figure_size=(7, 0.25 * len(present) + 1.5))
        + scale_y_log10(limits=(0.8, upper), breaks=breaks, expand=[0, 0, 0, 0])
    )
    savefig_plotnine(panel_ratios(p_frac | p_n, [3, 1]), None,
                     os.path.join(sample_dir, "stringent_pure_fraction_per_cell_type.png"))

    # kNN neighbour fractions of the mixed cells resolved by kNN (unresolved ones are Others), long format;
    # a cell type outside a cell's mix (NA) has no neighbours
    knn = pd.read_csv(args.knn_fractions[i])
    frac_cols = [c for c in knn.columns if c not in ("cell_id", "knn_k")]
    knn = knn.merge(labels[sample_id][["cell_id", "cell_type"]], on="cell_id", how="left", validate="1:1")
    knn = knn[knn["cell_type"] != OTHERS]
    long = knn.melt(id_vars=["cell_id", "cell_type"], value_vars=frac_cols, var_name="neighbour",
                    value_name="fraction")
    long["fraction"] = long["fraction"].fillna(0)
    long["cell_type"] = pd.Categorical(long["cell_type"], categories=[ct for ct in ct_order if ct in set(long["cell_type"])])
    long["neighbour"] = pd.Categorical(long["neighbour"], categories=[ct for ct in ct_order if ct in frac_cols])
    own = long["neighbour"].astype(str) == long["cell_type"].astype(str)

    # neighbours of every other cell type, one facet per final cell type (up to 4 columns)
    out = os.path.join(sample_dir, "knn_neighbour_fractions.png")
    other = long[~own]
    if other.empty:
        empty_figure("no mixed cells resolved by kNN", title, out)
    else:
        ncol = min(4, other["cell_type"].nunique())
        nrow = -(-other["cell_type"].nunique() // ncol)
        p = (
            ggplot(other, aes("neighbour", "fraction", fill="neighbour"))
            + geom_violin(scale="width", width=0.9, size=0.2, show_legend=False)
            + geom_boxplot(width=0.15, outlier_shape="", fill="white", size=0.3)
            + facet_wrap("cell_type", ncol=ncol, scales="free_x")
            + scale_fill_manual(values=palette)
            + labs(x="neighbour cell type", y="fraction of k nearest clean neighbours",
                   title=f"{title}: mixed cells per final TACIT cell type")
            + theme_classic()
            + theme(panel_border=element_rect(fill=None, color="black"),
                    text=element_text(font="Helvetica"),
                    axis_text_x=element_text(rotation=45, ha="right"),
                    figure_size=(4 * ncol, 3.5 * nrow + 1))
            + scale_y_continuous(expand=[0, 0, 0.01, 0], limits=(0, 1))
        )
        savefig_plotnine(p, None, out)

    # neighbours of the final cell type, per final cell type
    out = os.path.join(sample_dir, "knn_assigned_fraction.png")
    if not own.any():
        empty_figure("no mixed cells resolved by kNN", title, out)
    else:
        p = (
            ggplot(long[own], aes("cell_type", "fraction", fill="cell_type"))
            + geom_boxplot(width=0.7, outlier_size=0.5, show_legend=False)
            + scale_fill_manual(values=palette)
            + labs(x="final TACIT cell type", y="fraction of k nearest clean neighbours\nof the final cell type",
                   title=f"{title}: mixed cells")
            + theme_classic()
            + theme(panel_border=element_rect(fill=None, color="black"),
                    text=element_text(font="Helvetica"),
                    axis_text_x=element_text(rotation=45, ha="right"),
                    figure_size=(0.45 * long.loc[own, "cell_type"].nunique() + 2, 5))
            + scale_y_continuous(expand=[0, 0, 0.01, 0], limits=(0, 1))
        )
        savefig_plotnine(p, None, out)

    # TACIT labels in feature table order (run_tacit.R keeps it)
    features = pd.read_csv(args.feature_tables[i])
    cells = features[["cell_id"]].merge(labels[sample_id], on="cell_id", how="left", validate="1:1")
    cell_types = cells["cell_type"]
    if cell_types.isna().any():
        raise ValueError(f"{sample_id}: {cell_types.isna().sum()} cells without a TACIT label")

    # mean marker z-score per cell type, all cells and unique cells only (not mixed, not unassigned)
    plot_heatmap(features, cell_types, signatures[i], args.cofactor, ct_order, palette,
                 f"{title}: mean marker z-score per TACIT cell type",
                 os.path.join(sample_dir, "mean_intensity_heatmap.png"))
    unique = (cell_types != OTHERS) & ~cells["mixed"].astype(bool)
    plot_heatmap(features, cell_types, signatures[i], args.cofactor, ct_order, palette,
                 f"{title}: mean marker z-score per TACIT cell type (unique cells)",
                 os.path.join(sample_dir, "mean_intensity_heatmap_unique.png"), keep=unique)
