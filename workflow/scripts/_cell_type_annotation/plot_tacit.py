"""Plot the TACIT cell types of one segmentation method over all samples.

Reads, per sample (same order in every list): the feature table of quantify_cells, the TACIT labels
of run_tacit.R and the signature CSV of the sample's panel. The tissue of a sample is the first part
of its sample ID (tissue/panel/donor/sample).
Writes, over all samples (in the folder of --out_fractions):
  - a stacked bar plot of the fraction of each cell type per sample, one facet per tissue
  - a stacked bar plot of the fraction of unique (TACIT not mixed), mixed and unassigned (Others) cells
    per sample, one facet per tissue
  - a boxplot of the fraction of the cells of each cell type that are unique, i.e. not reassigned from
    a mixed microcluster, across samples (one point per sample, coloured by tissue)
  - per signature (panel), a heatmap of the weight of every marker per cell type, markers that define a
    cell type and the other markers in separate blocks
  - a boxplot of the fraction of the mixed cells of a sample per pair of positive cell types (the
    N_PAIRS most frequent pairs over all samples), mixed cells positive for three or more cell types
    grouped in one box, across samples (one point per sample, coloured by tissue), below a stacked bar
    of the final TACIT cell types of the mixed cells of each pair (all samples pooled, unit scaled)
Writes, per sample (in the folder of its --out_sample_fractions):
  - a bar plot of the fraction of each cell type
  - a bar plot of the fraction of unique, mixed and unassigned cells
  - a bar plot of the fraction of the cells of each cell type that are unique
  - a bar plot of the number of mixed cells per pair of positive cell types (the N_PAIRS most frequent
    pairs of the sample, log scale), mixed cells positive for three or more cell types in one bar,
    below a stacked bar of the final TACIT cell types of the mixed cells of each pair (unit scaled)
  - the same for the mixed cells positive for three or more cell types, per combination of positive cell
    types (the N_PAIRS most frequent combinations of the sample)
  - a heatmap of the mean normalized intensity (arcsinh(x / cofactor), z-score within the sample) of
    every signature marker per cell type, markers split (and clustered separately) by whether their
    signature weight defines a cell type
  - the same heatmap with the unique cells only (z-scores still over all cells of the sample)
  - the UMAP of preprocess_protein (same cofactor) on the signature markers, coloured by cell type
"""

import argparse
import os
import sys


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample_ids", nargs="+", required=True, help="tissue/panel/donor/sample of every sample")
    parser.add_argument("--feature_tables", nargs="+", required=True, help="quantify_cells feature tables")
    parser.add_argument("--cell_types", nargs="+", required=True, help="run_tacit.R outputs")
    parser.add_argument("--signatures", nargs="+", required=True, help="signature CSV of each sample's panel")
    parser.add_argument("--cofactor", type=float, required=True, help="arcsinh cofactor")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--out_fractions", required=True,
                        help="stacked bar plot of the cell type fractions; the other summary figures go to its folder")
    parser.add_argument("--out_sample_fractions", nargs="+", required=True,
                        help="per sample, cell type fractions; the sample's other figures go to its folder")
    args = parser.parse_args()

    # every per-sample list must follow --sample_ids
    n = len(args.sample_ids)
    for name in ["feature_tables", "cell_types", "signatures", "out_sample_fractions"]:
        if len(getattr(args, name)) != n:
            raise ValueError(f"--{name}: {len(getattr(args, name))} values for {n} samples")
    return args


# numba (K-NN graph, UMAP) must not use more cores than the job has: set before scanpy is imported
ARGS = parse_args()
os.environ.setdefault("NUMBA_NUM_THREADS", str(ARGS.threads))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from plotnine import (aes, coord_flip, element_rect, element_text, facet_grid, geom_boxplot, geom_col, geom_jitter,
                      ggplot, guide_legend, guides, labs, scale_color_manual, scale_fill_manual, scale_x_discrete,
                      scale_y_continuous, theme, theme_classic)

# PyComplexHeatmap still calls the removed matplotlib.cm.get_cmap in some versions
import matplotlib.cm
if not hasattr(matplotlib.cm, "get_cmap"):
    matplotlib.cm.get_cmap = matplotlib.colormaps.get_cmap
import PyComplexHeatmap as pch

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(script_dir, "../utils/"))
from multimodal_readwrite import mesmer_feature_table_to_adata
from plotting import TACIT_mm_cols, tissue_cols, savefig, savefig_plotnine, set_mpl_rcparams
from preprocessing import preprocess_protein

OTHERS = "Others"
OTHERS_COLOR = "#BDBDBD"

# unique (TACIT not mixed), mixed and unassigned (Others) cells
ASSIGNMENT_COLORS = {"unique": "#4A6FA5", "mixed": "#C08A4E", "unassigned": OTHERS_COLOR}
# mixed cells positive for three or more cell types, one bar / box in the pairs plots
MULTI = "3+ cell types"
PAIR_COLORS = {"2 cell types": "#BFBFBF", MULTI: "#5F5F5F"}
# most frequent cell type pairs shown in the pairs plots
N_PAIRS = 20


def sample_label(sample_id, sep="/"):
    """Short sample name for the plots: donor/sample (sep="\\n" for two rows)."""
    return sep.join(sample_id.split("/")[2:])


def cell_type_palette(signatures):
    """Cell type order (signature order, then Others) and one colour per cell type."""
    # union of the cell types of all signatures, in order of appearance
    order = list(dict.fromkeys(ct for sig in signatures for ct in sig.index))
    # fixed colours for the known cell types, tab20 / husl for any other signature cell type
    unknown = [ct for ct in order if ct not in TACIT_mm_cols]
    colors = sns.color_palette("tab20" if len(unknown) <= 20 else "husl", len(unknown)).as_hex()
    palette = {ct: TACIT_mm_cols[ct] for ct in order if ct in TACIT_mm_cols} | dict(zip(unknown, colors))
    palette[OTHERS] = OTHERS_COLOR
    return order + [OTHERS], palette


def all_labels(labels, ct_order):
    """TACIT labels of all samples in one table, with the tissue, the assignment (unique / mixed /
    unassigned) and, for mixed cells, the combination of positive cell types (signature order) and the pair
    (the combination of two cell types, MULTI for 3+)."""
    df = pd.concat([lab.assign(sample_id=s, tissue=s.split("/")[0]) for s, lab in labels.items()], ignore_index=True)
    df["positive_cell_types"] = df["positive_cell_types"].fillna("")
    df["assignment"] = np.where(df["cell_type"] == OTHERS, "unassigned", np.where(df["mixed"], "mixed", "unique"))

    # combination and pair label of the mixed cells
    rank = {ct: i for i, ct in enumerate(ct_order)}
    positives = df["positive_cell_types"].str.split("::")
    df["combination"] = [" + ".join(sorted(p, key=lambda ct: rank.get(ct, len(rank)))) if m else None
                         for p, m in zip(positives, df["mixed"].astype(bool))]
    df["pair"] = [(c if len(p) == 2 else MULTI) if c is not None else None for c, p in zip(df["combination"], positives)]
    return df


def tissue_palette(df):
    """One colour per tissue, in order of appearance."""
    tissues = list(dict.fromkeys(df["tissue"]))
    # fixed colours for the known tissues, Dark2 for any other tissue
    unknown = [t for t in tissues if t not in tissue_cols]
    colors = sns.color_palette("Dark2", len(unknown)).as_hex()
    return {t: tissue_cols[t] for t in tissues if t in tissue_cols} | dict(zip(unknown, colors))


def stacked_by_tissue(counts, fill, values, breaks, fill_title, out):
    """Stacked vertical bar per sample (fraction of cells), one facet per tissue (width by its number
    of samples), donor and sample on two rows below each bar."""
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


def summary_fractions(df, ct_order, palette, out):
    """Stacked bar of the fraction of each cell type, one bar per sample, one facet per tissue."""
    # cells and fraction per sample and cell type
    counts = df.groupby(["tissue", "sample_id", "cell_type"]).size().rename("cells").reset_index()
    counts["fraction"] = counts["cells"] / counts.groupby("sample_id")["cells"].transform("sum")

    # samples and tissues in the given order; cell types stacked top to bottom in signature order
    present = [ct for ct in ct_order if ct in set(counts["cell_type"])]
    counts["sample_id"] = pd.Categorical(counts["sample_id"], categories=list(dict.fromkeys(df["sample_id"])))
    counts["tissue"] = pd.Categorical(counts["tissue"], categories=list(dict.fromkeys(df["tissue"])))
    counts["cell_type"] = pd.Categorical(counts["cell_type"], categories=present)
    stacked_by_tissue(counts, "cell_type", palette, present, "TACIT cell type", out)


def summary_mixed(df, out):
    """Stacked bar of the fraction of unique, mixed and unassigned cells, one bar per sample, one facet per tissue."""
    # cells and fraction per sample and assignment
    counts = df.groupby(["tissue", "sample_id", "assignment"]).size().rename("cells").reset_index()
    counts["fraction"] = counts["cells"] / counts.groupby("sample_id")["cells"].transform("sum")

    # samples and tissues in the given order; unique on top
    counts["sample_id"] = pd.Categorical(counts["sample_id"], categories=list(dict.fromkeys(df["sample_id"])))
    counts["tissue"] = pd.Categorical(counts["tissue"], categories=list(dict.fromkeys(df["tissue"])))
    counts["assignment"] = pd.Categorical(counts["assignment"], categories=list(ASSIGNMENT_COLORS))
    stacked_by_tissue(counts, "assignment", ASSIGNMENT_COLORS, list(ASSIGNMENT_COLORS), "TACIT assignment", out)


def boxplot_across_samples(data, x, fill, fill_values, fill_title, tissues, y_title, y_max=1):
    """Boxplot of one value per sample for every x category, boxes filled by fill and jittered points
    (one per sample) coloured by tissue, y axis from 0 to y_max (None: data maximum)."""
    p = (
        ggplot(data, aes(x, "fraction"))
        + geom_boxplot(aes(fill=fill), outlier_shape="", width=0.7)
        + geom_jitter(aes(color="tissue"), width=0.2, height=0, size=1.2, alpha=0.8, random_state=0)
        + scale_fill_manual(values=fill_values)
        + scale_color_manual(values=tissues, breaks=list(tissues))
        + guides(fill=guide_legend(order=1), color=guide_legend(order=2, override_aes={"size": 3, "alpha": 1}))
        + labs(x="", y=y_title, fill=fill_title, color="tissue")
        + theme_classic()
        + theme(panel_border=element_rect(fill=None, color="black"),
                text=element_text(font="Helvetica"),
                axis_text_x=element_text(rotation=45, ha="right"),
                figure_size=(0.45 * data[x].nunique() + 4, 5))
        + scale_y_continuous(expand=[0, 0, 0.01, 0], limits=(0, y_max))
    )
    return p


def summary_unique(df, ct_order, palette, out):
    """Fraction of the cells of each cell type that are unique (not mixed) across samples."""
    # unique fraction per sample and final cell type (Others are never mixed)
    assigned = df[df["cell_type"] != OTHERS]
    frac = (assigned.assign(unique=~assigned["mixed"].astype(bool))
            .groupby(["tissue", "sample_id", "cell_type"])["unique"].mean().rename("fraction").reset_index())

    # cell types left to right in signature order
    present = [ct for ct in ct_order if ct in set(frac["cell_type"])]
    frac["cell_type"] = pd.Categorical(frac["cell_type"], categories=present)
    p = boxplot_across_samples(frac, "cell_type", "cell_type", palette, "TACIT cell type", tissue_palette(df),
                               "fraction of unique cells (not mixed)")
    savefig_plotnine(p, None, out)


def pair_labels(mixed, col, order, ct_order):
    """Fraction of the final TACIT cell types (columns, signature order) of the mixed cells of every value of
    col (rows, in order: pairs or combinations), unit scaled per row."""
    wide = pd.crosstab(mixed[col], mixed["cell_type"], normalize="index")
    return wide.reindex(index=order, columns=[ct for ct in ct_order if ct in wide.columns], fill_value=0)


def pairs_figure(order, draw_counts, labels, palette, title, out):
    """Vertical panels aligned on the pairs: a stacked bar of the final TACIT cell types of the mixed cells
    of each pair above draw_counts(ax, x), twice as tall."""
    x = np.arange(len(order))
    fig, (ax_top, ax_bot) = plt.subplots(2, 1, sharex=True, figsize=(0.32 * len(order) + 4, 5.5),
                                         gridspec_kw={"height_ratios": [1, 2], "hspace": 0.08})

    # final cell types stacked with the first signature cell type on top
    bottom = np.zeros(len(order))
    for ct in labels.columns[::-1]:
        ax_top.bar(x, labels[ct], bottom=bottom, width=0.8, color=palette[ct], label=ct)
        bottom += labels[ct].to_numpy()
    ax_top.set(ylim=(0, 1), xlim=(-0.6, len(order) - 0.4))
    ax_top.set_ylabel("final cell type\nfraction", fontsize=13)
    ax_top.set_title(title, fontsize=14)
    ax_top.tick_params(axis="x", length=0)
    # legend in signature order (reverse of the stacking order), right of the stacked bars
    handles, names = ax_top.get_legend_handles_labels()
    ax_top.legend(handles[::-1], names[::-1], title="final TACIT cell type", bbox_to_anchor=(1.01, 1),
                  loc="upper left", frameon=False, fontsize=10, title_fontsize=11, labelspacing=0.3)

    # counts / fractions per pair below, pair names on the x axis
    draw_counts(ax_bot, x)
    ax_bot.set_xticks(x, order, rotation=45, ha="right", rotation_mode="anchor")
    for ax in (ax_top, ax_bot):
        ax.tick_params(labelsize=12)
    savefig(None, out)


def summary_pairs(df, ct_order, palette, out):
    """Fraction of the mixed cells of a sample per pair of positive cell types (N_PAIRS most frequent
    over all samples, 3+ cell types in one box) across samples, log scale, below the final TACIT cell
    types of the mixed cells of each pair (all samples pooled)."""
    mixed = df[df["pair"].notna()]

    # N_PAIRS most frequent pairs over all samples left to right, the 3+ box last
    totals = mixed["pair"].value_counts()
    order = [p for p in totals.index if p != MULTI][:N_PAIRS] + ([MULTI] if MULTI in totals.index else [])

    # fraction of each sample's mixed cells per pair; samples without a pair left out (log scale)
    counts = pd.crosstab(mixed["sample_id"], mixed["pair"])
    frac = counts.div(counts.sum(axis=1), axis=0).reindex(columns=order, fill_value=0)
    frac = frac.stack().rename("fraction").reset_index()
    frac = frac[frac["fraction"] > 0]
    frac["tissue"] = frac["sample_id"].str.split("/").str[0]
    tissues = tissue_palette(df)

    def draw_counts(ax, x):
        # one box per pair, filled by 2 / 3+ cell types
        boxes = ax.boxplot([frac.loc[frac["pair"] == p, "fraction"] for p in order], positions=x, widths=0.7,
                           patch_artist=True, showfliers=False, medianprops={"color": "black"})
        for box, p in zip(boxes["boxes"], order):
            box.set(facecolor=PAIR_COLORS[MULTI if p == MULTI else "2 cell types"], edgecolor="black")
        # one jittered point per sample, coloured by tissue
        rng = np.random.default_rng(0)
        pos = {p: i for i, p in enumerate(order)}
        for tissue, color in tissues.items():
            sub = frac[frac["tissue"] == tissue]
            ax.scatter(sub["pair"].map(pos) + rng.uniform(-0.2, 0.2, len(sub)), sub["fraction"], s=12,
                       color=color, alpha=0.8, linewidths=0, zorder=3, label=tissue)
        ax.set_yscale("log")
        ax.set_ylabel("fraction of mixed cells\n(log scale)", fontsize=13)

        # tissue legend at the bottom right, clear of the cell type legend above
        ax.legend(title="tissue", bbox_to_anchor=(1.01, 0), loc="lower left", frameon=False, fontsize=10,
                  title_fontsize=11, markerscale=2)

    pairs_figure(order, draw_counts, pair_labels(mixed, "pair", order, ct_order), palette, "", out)


def sample_bars(data, x, order, fill_values, x_title, y_title, title, y_scale=None):
    """Horizontal bar per category of one sample, categories top to bottom in order."""
    data = data.assign(**{x: pd.Categorical(data[x], categories=order[::-1])})
    p = (
        ggplot(data, aes(x, "value", fill=x))
        + geom_col(width=0.8, show_legend=False)
        + coord_flip()
        + scale_fill_manual(values=fill_values)
        + labs(x=x_title, y=y_title, title=title)
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
        plt.figure(figsize=(6, 2))
        plt.text(0.5, 0.5, empty_text, ha="center", va="center", fontsize=13)
        plt.axis("off")
        plt.title(title)
        savefig(None, out)
        return
    order = [p for p in counts.index if p != MULTI][:N_PAIRS] + ([MULTI] if MULTI in counts.index else [])
    counts = counts.reindex(order)

    def draw_counts(ax, x):
        # mixed cells per pair / combination (log scale), dark for 3+ cell types, count above each bar
        colors = [PAIR_COLORS["2 cell types"] if p != MULTI and p.count(" + ") == 1 else PAIR_COLORS[MULTI]
                  for p in order]
        ax.bar(x, counts, width=0.8, color=colors)
        ax.set_yscale("log")
        ax.set_ylim(0.8, counts.max() * 5)
        for xi, v in zip(x, counts):
            ax.text(xi, v * 1.15, str(v), ha="center", va="bottom", rotation=90, fontsize=10)
        ax.set_ylabel("mixed cells (log scale)", fontsize=13)

    pairs_figure(order, draw_counts, pair_labels(cells, col, order, ct_order), palette, title, out)


def sample_plots(sub, ct_order, palette, out_fractions, out_mixed, out_unique, out_pairs, out_multi):
    """Per sample bar plots: cell type fractions, unique / mixed / unassigned fractions, unique fraction
    per cell type, mixed cells per cell type pair and per combination of 3+ cell types."""
    title = sample_label(sub["sample_id"].iloc[0])

    # fraction of each cell type
    frac = sub["cell_type"].value_counts(normalize=True).rename("value").rename_axis("cell_type").reset_index()
    present = [ct for ct in ct_order if ct in set(frac["cell_type"])]
    savefig_plotnine(sample_bars(frac, "cell_type", present, palette, "", "fraction of cells", title), None,
                     out_fractions)

    # fraction of unique, mixed and unassigned cells
    frac = sub["assignment"].value_counts(normalize=True).rename("value").rename_axis("assignment").reset_index()
    order = [a for a in ASSIGNMENT_COLORS if a in set(frac["assignment"])]
    p = sample_bars(frac, "assignment", order, ASSIGNMENT_COLORS, "", "fraction of cells", title,
                    y_scale=scale_y_continuous(expand=[0, 0, 0.01, 0], limits=(0, 1)))
    savefig_plotnine(p, None, out_mixed)

    # unique fraction per final cell type (Others are never mixed)
    assigned = sub[sub["cell_type"] != OTHERS]
    frac = (~assigned["mixed"].astype(bool)).groupby(assigned["cell_type"]).mean().rename("value").reset_index()
    present = [ct for ct in ct_order if ct in set(frac["cell_type"])]
    p = sample_bars(frac, "cell_type", present, palette, "", "fraction of unique cells (not mixed)", title,
                    y_scale=scale_y_continuous(expand=[0, 0, 0.01, 0], limits=(0, 1)))
    savefig_plotnine(p, None, out_unique)

    # mixed cells per pair (N_PAIRS most frequent, 3+ cell types in one bar) and per combination of 3+ cell types
    sample_mixed_counts(sub[sub["pair"].notna()], "pair", ct_order, palette, title, "no mixed cells", out_pairs)
    sample_mixed_counts(sub[sub["pair"] == MULTI], "combination", ct_order, palette, title,
                        "no mixed cells positive for 3+ cell types", out_multi)


def normalize(features, markers, cofactor):
    """arcsinh(x / cofactor) and z-score per marker within the sample, as run_tacit.R (without its shift to min 0)."""
    norm = np.arcsinh(features[markers] / cofactor)
    return (norm - norm.mean()) / norm.std(ddof=1)


def plot_heatmap(features, cell_types, signature, cofactor, ct_order, palette, title, out, keep=None):
    """Mean normalized intensity of every signature marker per cell type, cell type colours on the left,
    cell numbers on the right; markers split (and clustered separately) by whether they define a cell type
    in the signature (a non-empty, non-zero weight), shown as a column annotation. With keep (boolean mask
    of the cells), only those cells are averaged, z-scores still over all cells of the sample."""
    markers = list(signature.columns)
    # z-score over all cells, then only the kept cells
    norm = normalize(features, markers, cofactor)
    if keep is not None:
        keep = np.asarray(keep, dtype=bool)
        norm, cell_types = norm[keep], cell_types[keep]
    # mean z-score per cell type and marker, rows in signature order
    means = norm.groupby(cell_types.to_numpy()).mean()
    means = means.loc[[ct for ct in ct_order if ct in means.index]]

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
    # cell type colour (same palette as the fractions and UMAP plots) as a row annotation
    ct = pd.DataFrame({"cell type": means.index}, index=means.index)
    left = pch.HeatmapAnnotation(df=ct, colors={"cell type": {c: palette[c] for c in means.index}},
                                 label_kws={"visible": False}, plot_legend=False, axis=0, verbose=0)

    plt.figure(figsize=(0.2 * len(markers) + 2.5, 0.16 * len(means) + 1.6))
    # rectangular colorbar (PyComplexHeatmap defaults to extend="both", i.e. arrow ends)
    cm = pch.ClusterMapPlotter(data=means, row_cluster=False, col_cluster=True, col_dendrogram=False,
                               col_split=split, col_split_order=["yes", "no"], col_split_gap=1.5,
                               top_annotation=top, left_annotation=left, right_annotation=right,
                               show_rownames=True, show_colnames=True,
                               row_names_side="left", cmap="RdBu_r", center=0, vmin=-2, vmax=2,
                               label="mean z-score", legend_kws={"extend": "neither"}, legend_hpad=3,
                               xticklabels_kws={"labelrotation": 45},
                               xlabel="Protein", ylabel="TACIT automated cell types", ylabel_kws={"labelpad": 110},
                               verbose=0)
    # PyComplexHeatmap hides the row names when there are left and right annotations: write them
    # next to the cell type colour strip
    strip = cm.left_annotation.axes[0, 0]
    strip.yaxis.set_visible(True)
    strip.set_yticks(cm.heatmap_axes[0, 0].get_yticks(), labels=means.index)
    strip.tick_params(axis="y", left=False, labelleft=True, labelsize=plt.rcParams["ytick.labelsize"])
    # colorbar title centered above the bar instead of along its side
    cbar = next(c for c in cm.cbars if isinstance(c, matplotlib.colorbar.Colorbar))
    cbar.set_label("")
    cbar.ax.set_title("mean z-score", loc="center", fontsize=13, pad=8)
    # colorbar 0.8 inch right of the cell number bars, room for its centered title (PyComplexHeatmap leaves more)
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


def zero_white_cmap(vmin):
    """RdBu_r stretched so that vmin (< 0) is its darkest blue, 0 white and 1 its darkest red (white to red
    when vmin is 0), for a linear scale from vmin to the maximum weight."""
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


def plot_signature(signature, palette, title, out):
    """Signature weights (cell type x marker), cell type colours on the left, markers that define a cell
    type (non-empty, non-zero weight) and the other markers in separate blocks."""
    weights = signature.fillna(0).astype(float)
    if weights.empty:
        plt.figure(figsize=(6, 2))
        plt.text(0.5, 0.5, "no cell types in the signature", ha="center", va="center", fontsize=13)
        plt.axis("off")
        plt.title(title)
        savefig(None, out)
        return

    # used markers ordered by the first cell type they define (block diagonal), then the unused ones
    used = (weights != 0).any(axis=0)
    first = weights.loc[:, used].ne(0).to_numpy().argmax(axis=0)
    order = list(weights.columns[used][np.argsort(first, kind="stable")]) + list(weights.columns[~used])
    weights = weights[order]

    # used for phenotyping as a column annotation, columns split by it when both occur
    usage = pd.Series(np.where(used[order], "yes", "no"), index=order)
    top = pch.HeatmapAnnotation(df=usage.to_frame("Used for phenotyping"),
                                colors={"Used for phenotyping": {"yes": "#4D4D4D", "no": "#D9D9D9"}},
                                label_kws={"visible": False}, axis=1, verbose=0)
    split = usage if usage.nunique() > 1 else None
    # cell type colour as a row annotation
    ct = pd.DataFrame({"cell type": weights.index}, index=weights.index)
    left = pch.HeatmapAnnotation(df=ct, colors={"cell type": {c: palette[c] for c in weights.index}},
                                 label_kws={"visible": False}, plot_legend=False, axis=0, verbose=0)
    # non-zero weights written in their cells
    annot = weights.map(lambda w: f"{w:g}" if w != 0 else "")
    # 0 white, 1 red, the lowest negative weight blue
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
    # cell type names next to the cell type colour strip instead of right of the last block
    for ax in cm.heatmap_axes.flat:
        ax.tick_params(axis="y", right=False, labelright=False)
    strip = cm.left_annotation.axes[0, 0]
    strip.yaxis.set_visible(True)
    strip.set_yticks(cm.heatmap_axes[0, 0].get_yticks(), labels=weights.index)
    strip.tick_params(axis="y", left=False, labelleft=True, labelsize=plt.rcParams["ytick.labelsize"])
    # colorbar title centered above the bar, 0.4 inch right of the heatmap
    cbar = next(c for c in cm.cbars if isinstance(c, matplotlib.colorbar.Colorbar))
    cbar.set_label("")
    cbar.ax.set_title("weight", loc="center", fontsize=13, pad=8)
    # 0 always labelled, with the lowest and highest weight
    cbar.set_ticks([vmin, 0, vmax] if vmin < 0 else [0, vmax / 2, vmax])
    right = cm.heatmap_axes[0, -1].get_position().x1
    pos = cbar.ax.get_position()
    cbar.ax.set_position([right + 0.6 / plt.gcf().get_figwidth(), pos.y0, pos.width, pos.height])
    # used for phenotyping legend below the colorbar, no frame
    for legend in cm.cbars:
        if isinstance(legend, matplotlib.legend.Legend):
            legend.set_frame_on(False)
            legend.set_bbox_to_anchor((cbar.ax.get_position().x0 - 0.3 / plt.gcf().get_figwidth(),
                                       cbar.ax.get_position().y0 - 0.1), transform=plt.gcf().transFigure)
            legend.set_loc("upper left")
    plt.suptitle(title, y=1.02)
    savefig(None, out)


def umap_coordinates(feature_table, markers, cofactor):
    """preprocess_protein UMAP of one sample on the signature markers."""
    adata = mesmer_feature_table_to_adata(feature_table, dapi_to_obs=False)

    # UMAP on exactly the signature markers (exclude=[]: no preprocess_protein default exclusions)
    processed = preprocess_protein(adata, arcsinh_cofactor=cofactor, include=markers, exclude=[])
    return pd.DataFrame({"cell_id": adata.obs["cell_id"].to_numpy(),
                         "UMAP1": processed.obsm["X_umap"][:, 0], "UMAP2": processed.obsm["X_umap"][:, 1]})


def plot_umap(coords, cell_types, ct_order, palette, title, out):
    """UMAP of one sample coloured by TACIT cell type."""
    df = coords.assign(cell_type=cell_types.to_numpy())
    hue_order = [ct for ct in ct_order if ct in set(df["cell_type"])]

    # square UMAP, legend on the right
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.set_box_aspect(1)
    sns.scatterplot(data=df, x="UMAP1", y="UMAP2", hue="cell_type", hue_order=hue_order, palette=palette, s=0.5,
                    linewidth=0, rasterized=True, ax=ax)
    ax.set(xlabel="UMAP1", ylabel="UMAP2", xticks=[], yticks=[])
    ax.set_title(title, fontsize=13)
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", markerscale=8, frameon=False, title="TACIT cell type")
    savefig(None, out)


def main(args):
    set_mpl_rcparams()

    # signatures (cell_type x marker) and TACIT labels of every sample
    signatures = [pd.read_csv(path).set_index("cell_type") for path in args.signatures]
    labels = {s: pd.read_csv(path) for s, path in zip(args.sample_ids, args.cell_types)}
    ct_order, palette = cell_type_palette(signatures)
    df = all_labels(labels, ct_order)

    # summaries over all samples: cell type and unique / mixed / unassigned fractions per tissue,
    # unique fraction per cell type and mixed cell fraction per cell type pair across samples
    summary_dir = os.path.dirname(args.out_fractions)
    summary_fractions(df, ct_order, palette, args.out_fractions)
    summary_mixed(df, os.path.join(summary_dir, "mixed_fraction.png"))
    summary_unique(df, ct_order, palette, os.path.join(summary_dir, "unique_fraction_per_cell_type.png"))
    summary_pairs(df, ct_order, palette, os.path.join(summary_dir, "mixed_cell_type_pairs.png"))

    # weights of every distinct signature (signature_{panel}.png)
    for path in dict.fromkeys(args.signatures):
        panel = os.path.splitext(os.path.basename(path))[0]
        plot_signature(pd.read_csv(path).set_index("cell_type"), palette, f"{panel}: TACIT signature",
                       os.path.join(summary_dir, f"signature_{panel}.png"))

    for i, sample_id in enumerate(args.sample_ids):
        print(f"{sample_id} ({i + 1}/{len(args.sample_ids)})", flush=True)
        sample_dir = os.path.dirname(args.out_sample_fractions[i])

        # bar plots of the sample's labels
        sample_plots(df[df["sample_id"] == sample_id], ct_order, palette, args.out_sample_fractions[i],
                     os.path.join(sample_dir, "mixed_fraction.png"),
                     os.path.join(sample_dir, "unique_fraction_per_cell_type.png"),
                     os.path.join(sample_dir, "mixed_cell_type_pairs.png"),
                     os.path.join(sample_dir, "mixed_cell_type_3plus.png"))

        features = pd.read_csv(args.feature_tables[i])
        markers = list(signatures[i].columns)

        # TACIT labels in feature table order (run_tacit.R keeps it)
        cells = features[["cell_id"]].merge(labels[sample_id], on="cell_id", how="left", validate="1:1")
        cell_types = cells["cell_type"]
        if cell_types.isna().any():
            raise ValueError(f"{sample_id}: {cell_types.isna().sum()} cells without a TACIT label")

        # mean marker intensity per cell type
        plot_heatmap(features, cell_types, signatures[i], args.cofactor, ct_order, palette,
                     f"{sample_label(sample_id)}: mean marker intensity per TACIT cell type",
                     os.path.join(sample_dir, "mean_intensity_heatmap.png"))

        # the same with the unique cells only (not mixed, not unassigned)
        unique = (cell_types != OTHERS) & ~cells["mixed"].astype(bool)
        plot_heatmap(features, cell_types, signatures[i], args.cofactor, ct_order, palette,
                     f"{sample_label(sample_id)}: mean marker intensity per TACIT cell type (unique cells)",
                     os.path.join(sample_dir, "mean_intensity_heatmap_unique.png"), keep=unique)

        # # preprocess_protein UMAP coloured by cell type
        # coords = umap_coordinates(args.feature_tables[i], markers, args.cofactor)
        # # donor, sample and cell number on separate lines when the title is long
        # title = f"{sample_label(sample_id)} ({len(coords)} cells)"
        # if len(title) > 50:
        #     title = f"{sample_label(sample_id, sep='/' + chr(10))}\n({len(coords)} cells)"
        # plot_umap(coords, cell_types, ct_order, palette, title, os.path.join(sample_dir, "umap.png"))


if __name__ == "__main__":
    main(ARGS)
