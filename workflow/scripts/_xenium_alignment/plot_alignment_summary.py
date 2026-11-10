"""Summarise the alignment to Xenium over all samples: fraction of Xenium cells matched to a proteomics cell.

Reads the cell counts tables written by plot_alignment.py ({seg_method}/counts.csv, matches to the
default Xenium segmentation, named by --xenium_seg: nuclear_expansion, or multimodal with the 10x
multimodal cell segmentation kit), one per alignment method, max_width, sample and segmentation
method; the wildcards are recovered from the paths, relative to --root (OUT/xenium_alignment):
  {align_method}/max_width_{max_width}/{sample_id}/figures/{seg_method}/counts.csv
For every sample whose Xenium run also has a Proseg segmentation
(xenium/processed/segmentation/proseg/<raw run path>/raw_results), the proteomics cells are matched
to the Proseg cells as well, with the same functions as plot_alignment.py. Samples without a Proseg
folder are skipped.
Writes:
  - the table of matched fractions (one row per counts table)
  - per Xenium segmentation (--xenium_seg, proseg) and proteomics segmentation method: a bar plot of
    the matched fraction per sample, one bar per max_width, below the number of Xenium and
    proteomics cells of each sample (fraction_matched_per_sample_{xenium_seg}_{seg_method}.png)
  - per Xenium segmentation: the same at --summary_max_width only, one bar per proteomics
    segmentation method (fraction_matched_per_seg_method_{xenium_seg}.png)
  - a box plot of the matched fraction across samples per max_width
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from mizani.labels import label_comma
from plotnine import (aes, element_blank, element_line, element_rect, element_text, facet_wrap, geom_boxplot, geom_col, geom_jitter,
                      ggplot, position_dodge, scale_fill_manual, scale_y_continuous, theme, theme_classic, xlab, ylab)

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(script_dir, "../utils/"))
from experiments import get_samples, load_experiments
from plotting import FONT, savefig, savefig_plotnine, set_mpl_rcparams
from plot_alignment import align_proteomics_cells, match_cells, read_xenium_cells

MAX_WIDTH_COLORS = ["#19618A", "#09B39C", "#FFB300", "#D32F2F", "#7B1FA2", "#5D4037"]
# Xenium (neutral grey), then the proteomics segmentation methods (sorted) in complementary blue / orange,
# distinct from the max_width colours
MODALITY_COLORS = ["#4D4D4D", "#4C72E0", "#F28E2B", "#E15FA4", "#76B041"]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", nargs="+", required=True, help="plot_alignment counts tables")
    parser.add_argument("--root", required=True, help="OUT/xenium_alignment, the counts paths are parsed from there")
    parser.add_argument("--seg_root", required=True, help="OUT/segmentation, with {seg_method}/{sample_id}/feature_table.csv")
    parser.add_argument("--experiments", required=True, help="experiments.yml, for the Xenium run of each sample")
    parser.add_argument("--max_dist", type=float, default=5, help="max distance (µm) between matched cells")
    parser.add_argument("--xenium_seg", default="multimodal",
                        help="name of the default Xenium segmentation of the counts tables (e.g. nuclear_expansion)")
    parser.add_argument("--summary_max_width", default="codex",
                        help="max_width of the plot comparing the proteomics segmentation methods")
    parser.add_argument("--out_table", required=True)
    parser.add_argument("--out_figures_dir", required=True, help="per-sample bar plots are written here")
    parser.add_argument("--out_boxplot", required=True)
    return parser.parse_args()


def read_counts(path, root, xenium_seg):
    """One row: wildcards of the counts table and its number of Xenium / matched cells."""
    # wildcards from the path relative to root
    parts = os.path.relpath(path, root).split(os.sep)
    align_method, max_width, seg_method = parts[0], parts[1].removeprefix("max_width_"), parts[-2]
    sample_id = "/".join(parts[2:-3])
    # cell numbers per modality
    ncells = pd.read_csv(path, index_col=0).set_index("modality")["ncells"]
    return {
        "align_method": align_method,
        "max_width": max_width,
        "sample_id": sample_id,
        "seg_method": seg_method,
        "xenium_seg": xenium_seg,
        "n_xenium": ncells["10x_mm_5um"],
        "n_proteomics": ncells["proteomics"],
        "n_matched": ncells["matched"],
    }


def proseg_dir(xenium_dir):
    """Proseg raw_results of a Xenium run (xenium/raw/... -> xenium/processed/segmentation/proseg/...)."""
    raw_base, run = xenium_dir.rsplit("/raw/", 1)
    return f"{raw_base}/processed/segmentation/proseg/{run}/raw_results"


def proseg_counts(df, samples, args):
    """Rows as read_counts, matching the proteomics cells to the Proseg cells of each Xenium run, for
    every alignment method / max_width / segmentation method of df (the default Xenium segmentation rows)."""
    rows = []
    for sample_id, df_sample in df.groupby("sample_id"):
        # skip samples whose Xenium run has no Proseg segmentation
        pdir = proseg_dir(samples[sample_id]["xenium"])
        if not os.path.isdir(pdir):
            print(f"no Proseg segmentation for {sample_id} ({pdir}), skipped")
            continue

        # Proseg cells, read once per sample
        xenium = read_xenium_cells(pdir)

        # aligned proteomics cells matched to the Proseg cells, per alignment / max_width / segmentation
        for r in df_sample.itertuples():
            feature_table = os.path.join(args.seg_root, r.seg_method, sample_id, "feature_table.csv")
            affine = os.path.join(args.root, r.align_method, f"max_width_{r.max_width}", sample_id, "affine_matrix.csv")
            cells, aligned_um = align_proteomics_cells(feature_table, affine, pdir)
            matches = match_cells(cells, aligned_um, xenium, args.max_dist)
            rows.append({
                "align_method": r.align_method,
                "max_width": r.max_width,
                "sample_id": sample_id,
                "seg_method": r.seg_method,
                "xenium_seg": "proseg",
                "n_xenium": len(xenium),
                "n_proteomics": len(cells),
                "n_matched": len(matches),
            })
    print(f"Proseg segmentation found for {len({r['sample_id'] for r in rows})}/{df['sample_id'].nunique()} samples")
    return rows


def stacked_barplot(ncells, fraction, fill, fill_colors, modality_colors, bar_width, out):
    """Number of Xenium / proteomics cells per sample (ncells: sample, ncells, modality) stacked on the
    matched fraction per sample (fraction: sample, fraction_matched, fill), one bar per fill value;
    the title is the "method" column of ncells."""
    # number of Xenium and proteomics (segmented) cells per sample
    p_ncells = (
        ggplot(ncells, aes(x="sample", y="ncells", fill="modality"))
        + geom_col(position=position_dodge(width=0.8), width=0.8)
        + facet_wrap("method", ncol=1)
        + theme_classic()
        + theme(panel_border=element_rect(fill=None, color="black"),
                text=element_text(font=FONT),
                axis_text_x=element_blank(),
                axis_ticks_major_x=element_blank(),
                axis_title_y=element_text(size=14),
                legend_title=element_text(size=18),
                legend_text=element_text(size=16),
                strip_text=element_text(size=18),
                figure_size=(bar_width, 3))
        + scale_y_continuous(expand=[0, 0, 0.01, 0], labels=label_comma())
        + scale_fill_manual(values=modality_colors, name="")
        + xlab("")
        + ylab("Number of\ncells")
    )

    # matched fraction per sample, one bar per fill value (strip already on the counts plot)
    p_fraction = (
        ggplot(fraction, aes(x="sample", y="fraction_matched", fill=fill))
        + geom_col(position=position_dodge(width=0.8), width=0.8)
        + theme_classic()
        + theme(panel_border=element_rect(fill=None, color="black"),
                text=element_text(font=FONT),
                axis_text_x=element_text(angle=45, hjust=1, color="black"),
                axis_title_y=element_text(size=14),
                legend_title=element_text(size=18),
                legend_text=element_text(size=16),
                figure_size=(bar_width, 4))
        + scale_y_continuous(expand=[0, 0, 0.01, 0], limits=(0, 1))
        + scale_fill_manual(values=fill_colors, name=fill)
        + xlab("")
        + ylab("Fraction of Xenium\ncells matched")
    )
    # counts stacked on top of the matched fractions, the composed figure sized as both plots together
    stacked = (p_ncells / p_fraction) & theme(figure_size=(bar_width, 7))
    savefig_plotnine(stacked, None, out)


def modality_cells(df, modality_colors):
    """Number of Xenium cells and of proteomics cells per segmentation method, per sample (one row per
    sample / modality), from rows of a single max_width."""
    # Xenium cells (the same for every segmentation method) and proteomics cells per segmentation method
    xenium = (df.groupby(["method", "sample"], observed=True)["n_xenium"].first().reset_index()
              .rename(columns={"n_xenium": "ncells"}).assign(modality="Xenium"))
    proteomics = (df[["method", "sample", "seg_method", "n_proteomics"]]
                  .rename(columns={"n_proteomics": "ncells", "seg_method": "modality"}))
    ncells = pd.concat([xenium, proteomics], ignore_index=True)
    # legend order of modality_colors, only the modalities present
    ncells["modality"] = (pd.Categorical(ncells["modality"], categories=list(modality_colors), ordered=True)
                          .remove_unused_categories())
    return ncells


def plot_per_sample(df, widths, colors, modality_colors, out):
    """One method combination: number of Xenium / proteomics cells per sample, stacked on the matched
    fraction per sample, one bar per max_width."""
    # width grows with the number of bars, under plotnine's 25 inch limit
    bar_width = min(24, max(4, 2 + 0.1 * df["sample"].nunique() * len(widths)))
    # cell numbers are the same for every max_width
    ncells = modality_cells(df[df["max_width"] == str(widths[0])], modality_colors)
    stacked_barplot(ncells, df, "max_width", colors, modality_colors, bar_width, out)


def plot_per_seg_method(df, max_width, modality_colors, out):
    """Number of Xenium / proteomics cells per sample (one bar per proteomics segmentation method),
    stacked on the matched fraction per sample, one bar per segmentation method, at one max_width."""
    df = df[df["max_width"] == max_width].copy()
    n_seg = df["seg_method"].nunique()
    bar_width = min(24, max(4, 2 + 0.1 * df["sample"].nunique() * (n_seg + 1)))
    # one title for all segmentation methods
    df["method"] = df["align_method"] + " / Xenium " + df["xenium_seg_label"] + f" / max_width {max_width}"
    ncells = modality_cells(df, modality_colors)
    df["seg_method"] = pd.Categorical(df["seg_method"], categories=list(modality_colors)[1:], ordered=True)
    seg_colors = {m: c for m, c in modality_colors.items() if m != "Xenium"}
    stacked_barplot(ncells, df, "seg_method", seg_colors, modality_colors, bar_width, out)


def placeholder(text, out):
    """Empty figure with a message, where the output is still expected."""
    plt.figure(figsize=(4, 2))
    plt.text(0.5, 0.5, text, ha="center", va="center")
    plt.axis("off")
    savefig(None, out)


def main():
    args = parse_args()
    set_mpl_rcparams()

    # matches to the default Xenium multimodal segmentation, per counts table
    df = pd.DataFrame([read_counts(p, args.root, args.xenium_seg) for p in args.counts])

    # matches to the Proseg Xenium segmentation, where it exists
    samples = get_samples(load_experiments(args.experiments), args.experiments)
    df = pd.concat([df, pd.DataFrame(proseg_counts(df, samples, args))], ignore_index=True)

    # matched fraction of Xenium cells
    df["fraction_matched"] = df["n_matched"] / df["n_xenium"]
    df = df.sort_values(["align_method", "seg_method", "xenium_seg", "sample_id", "max_width"]).reset_index(drop=True)
    df.to_csv(args.out_table, index=False)
    print(df.to_string(index=False))

    # short sample labels (sample name), the full sample_id only where the name is ambiguous
    names = pd.Series({s: s.split("/")[-1] for s in df["sample_id"].unique()})
    names[names.duplicated(keep=False)] = names.index[names.duplicated(keep=False)]
    df["sample"] = df["sample_id"].map(names)
    # categorical max_width (ordered numerically, "codex" last) for the fill / x axis, one facet per method
    # combination (alignment / proteomics segmentation / Xenium segmentation)
    widths = sorted(df["max_width"].unique(), key=lambda w: (not w.isdigit(), int(w) if w.isdigit() else 0))
    df["max_width"] = pd.Categorical(df["max_width"].astype(str), categories=[str(w) for w in widths], ordered=True)
    df["xenium_seg_label"] = df["xenium_seg"].str.replace("_", " ")
    df["method"] = df["align_method"] + " / " + df["seg_method"] + " / Xenium " + df["xenium_seg_label"]
    colors = dict(zip([str(w) for w in widths], MAX_WIDTH_COLORS * len(widths)))
    # Xenium then the proteomics segmentation methods, the same colours in every plot
    modalities = ["Xenium"] + sorted(df["seg_method"].unique())
    modality_colors = dict(zip(modalities, MODALITY_COLORS * len(modalities)))
    n_methods = df["method"].nunique()
    # samples ordered by their matched fraction at the smallest max_width (mean over methods), highest first
    order = (df[df["max_width"] == str(widths[0])].groupby("sample")["fraction_matched"].mean()
             .sort_values(ascending=False).index)
    df["sample"] = pd.Categorical(df["sample"], categories=order, ordered=True)

    # bar plots per Xenium segmentation, samples in the same order in all of them
    for xenium_seg in [args.xenium_seg, "proseg"]:
        df_seg = df[df["xenium_seg"] == xenium_seg]
        # one per proteomics segmentation method: one bar per max_width
        for seg_method in sorted(df["seg_method"].unique()):
            out = os.path.join(args.out_figures_dir, f"fraction_matched_per_sample_{xenium_seg}_{seg_method}.png")
            df_method = df_seg[df_seg["seg_method"] == seg_method]
            # no sample with this segmentation (e.g. no Proseg run): placeholder figure
            if df_method.empty:
                placeholder(f"no Xenium {xenium_seg} segmentation", out)
                continue
            plot_per_sample(df_method, widths, colors, modality_colors, out)

        # proteomics segmentation methods side by side, at summary_max_width
        out = os.path.join(args.out_figures_dir, f"fraction_matched_per_seg_method_{xenium_seg}.png")
        if args.summary_max_width not in df_seg["max_width"].astype(str).values:
            placeholder(f"no Xenium {xenium_seg} segmentation at max_width {args.summary_max_width}", out)
            continue
        plot_per_seg_method(df_seg, args.summary_max_width, modality_colors, out)

    # box plot: distribution of the matched fraction across samples, per max_width
    # strip label on two lines so it fits the narrow facets
    df["method_strip"] = df["align_method"] + " / " + df["seg_method"] + "\nXenium " + df["xenium_seg_label"]
    p = (
        ggplot(df, aes(x="max_width", y="fraction_matched", fill="max_width"))
        + geom_boxplot(outlier_shape="", alpha=0.6)
        + geom_jitter(width=0.15, height=0, size=1.5)
        + facet_wrap("method_strip", ncol=n_methods)
        + theme_classic()
        + theme(panel_border=element_rect(fill=None, color="black"),
                text=element_text(font=FONT),
                axis_text_x=element_text(color="black"),
                panel_grid_major=element_line(color="lightgrey", size=0.4),
                legend_position="none",
                figure_size=(max(3, 1 * len(widths)) * n_methods, 4))
        + scale_y_continuous(expand=[0, 0, 0.01, 0], limits=(0, 1))
        + scale_fill_manual(values=colors)
        + xlab("max_width (px)")
        + ylab("Fraction of Xenium cells matched")
    )
    savefig_plotnine(p, None, args.out_boxplot)
    print(f"written: {args.out_table}, {args.out_figures_dir}/fraction_matched_per_*.png, {args.out_boxplot}")


if __name__ == "__main__":
    main()
