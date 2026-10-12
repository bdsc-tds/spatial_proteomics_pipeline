"""Summarise the alignment to Xenium over all samples: fraction of Xenium cells matched to a proteomics cell.

Reads the cell counts tables written by plot_alignment.py ({seg_method}/counts.csv, matches to the
default Xenium multimodal segmentation), one per alignment method, max_width, sample and segmentation
method; the wildcards are recovered from the paths, relative to --root (OUT/xenium_alignment):
  {align_method}/max_width_{max_width}/{sample_id}/figures/{seg_method}/counts.csv
For every sample whose Xenium run also has a Proseg segmentation
(xenium/processed/segmentation/proseg/<raw run path>/raw_results), the proteomics cells are matched
to the Proseg cells as well, with the same functions as plot_alignment.py. Samples without a Proseg
folder are skipped. Each alignment method / proteomics segmentation / Xenium segmentation
(multimodal, proseg) combination is one facet.
Writes:
  - the table of matched fractions (one row per counts table)
  - a bar plot of the matched fraction per sample, one bar per max_width, below the number of
    Xenium and proteomics cells of each sample; one file per Xenium segmentation (multimodal, proseg)
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
from plotting import savefig, savefig_plotnine, set_mpl_rcparams
from plot_alignment import align_proteomics_cells, match_cells, read_xenium_cells

MAX_WIDTH_COLORS = ["#19618A", "#09B39C", "#FFB300", "#D32F2F", "#7B1FA2", "#5D4037"]
# Xenium / proteomics segmentation, distinct from the max_width colours of the plot below
MODALITY_COLORS = ["#7B1FA2", "#FFB300"]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", nargs="+", required=True, help="plot_alignment counts tables")
    parser.add_argument("--root", required=True, help="OUT/xenium_alignment, the counts paths are parsed from there")
    parser.add_argument("--seg_root", required=True, help="OUT/segmentation, with {seg_method}/{sample_id}/feature_table.csv")
    parser.add_argument("--experiments", required=True, help="experiments.yml, for the Xenium run of each sample")
    parser.add_argument("--max_dist", type=float, default=5, help="max distance (µm) between matched cells")
    parser.add_argument("--out_table", required=True)
    parser.add_argument("--out_barplot_multimodal", required=True, help="bar plot, Xenium multimodal segmentation")
    parser.add_argument("--out_barplot_proseg", required=True, help="bar plot, Proseg segmentation")
    parser.add_argument("--out_boxplot", required=True)
    return parser.parse_args()


def read_counts(path, root):
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
        "xenium_seg": "multimodal",
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
    every alignment method / max_width / segmentation method of df (the multimodal rows)."""
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


def plot_per_sample(df, widths, colors, out):
    """Number of Xenium / proteomics cells per sample, stacked on the matched fraction per sample (one
    bar per max_width), one facet per method combination."""
    n_methods = df["method"].nunique()
    # width grows with the number of bars, under plotnine's 25 inch limit
    bar_width = min(24, max(4, 2 + 0.1 * df["sample"].nunique() * len(widths)))

    # number of Xenium and proteomics (segmented) cells per sample, the same for every max_width
    ncells = (df[df["max_width"] == str(widths[0])]
              .melt(id_vars=["method", "seg_method", "sample"], value_vars=["n_xenium", "n_proteomics"],
                    var_name="modality", value_name="ncells"))
    ncells["modality"] = ncells["modality"].map({"n_xenium": "Xenium"}).fillna(ncells["seg_method"])
    modalities = ["Xenium"] + sorted(df["seg_method"].unique())
    ncells["modality"] = pd.Categorical(ncells["modality"], categories=modalities, ordered=True)
    p_ncells = (
        ggplot(ncells, aes(x="sample", y="ncells", fill="modality"))
        + geom_col(position=position_dodge(width=0.8), width=0.8)
        + facet_wrap("method", ncol=1)
        + theme_classic()
        + theme(panel_border=element_rect(fill=None, color="black"),
                text=element_text(font="Helvetica"),
                axis_text_x=element_blank(),
                axis_ticks_major_x=element_blank(),
                axis_title_y=element_text(size=14),
                legend_title=element_text(size=14),
                legend_text=element_text(size=12),
                strip_text=element_text(size=14),
                figure_size=(bar_width, 3 * n_methods))
        + scale_y_continuous(expand=[0, 0, 0.01, 0], labels=label_comma())
        + scale_fill_manual(values=dict(zip(modalities, MODALITY_COLORS * len(modalities))), name="")
        + xlab("")
        + ylab("Number of\ncells")
    )

    # bar plot: matched fraction per sample, one bar per max_width (strips already on the counts plot)
    p_fraction = (
        ggplot(df, aes(x="sample", y="fraction_matched", fill="max_width"))
        + geom_col(position=position_dodge(width=0.8), width=0.8)
        + facet_wrap("method", ncol=1)
        + theme_classic()
        + theme(panel_border=element_rect(fill=None, color="black"),
                text=element_text(font="Helvetica"),
                axis_text_x=element_text(angle=45, hjust=1, color="black"),
                axis_title_y=element_text(size=14),
                legend_title=element_text(size=14),
                legend_text=element_text(size=12),
                strip_text=element_blank(),
                strip_background=element_blank(),
                figure_size=(bar_width, 4 * n_methods))
        + scale_y_continuous(expand=[0, 0, 0.01, 0], limits=(0, 1))
        + scale_fill_manual(values=colors, name="max_width")
        + xlab("")
        + ylab("Fraction of Xenium\ncells matched")
    )
    # counts stacked on top of the matched fractions, the composed figure sized as both plots together
    stacked = (p_ncells / p_fraction) & theme(figure_size=(bar_width, 7 * n_methods))
    savefig_plotnine(stacked, None, out)


def main():
    args = parse_args()
    set_mpl_rcparams()

    # matches to the default Xenium multimodal segmentation, per counts table
    df = pd.DataFrame([read_counts(p, args.root) for p in args.counts])

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
    df["method"] = df["align_method"] + " / " + df["seg_method"] + " / Xenium " + df["xenium_seg"]
    colors = dict(zip([str(w) for w in widths], MAX_WIDTH_COLORS * len(widths)))
    n_methods = df["method"].nunique()
    # samples ordered by their matched fraction at the smallest max_width (mean over methods), highest first
    order = (df[df["max_width"] == str(widths[0])].groupby("sample")["fraction_matched"].mean()
             .sort_values(ascending=False).index)
    df["sample"] = pd.Categorical(df["sample"], categories=order, ordered=True)

    # bar plots: one per Xenium segmentation, samples in the same order in both
    out_barplots = {"multimodal": args.out_barplot_multimodal, "proseg": args.out_barplot_proseg}
    for xenium_seg, out in out_barplots.items():
        df_seg = df[df["xenium_seg"] == xenium_seg]
        # no sample with this segmentation (e.g. no Proseg run): placeholder figure, the output is still expected
        if df_seg.empty:
            plt.figure(figsize=(4, 2))
            plt.text(0.5, 0.5, f"no Xenium {xenium_seg} segmentation", ha="center", va="center")
            plt.axis("off")
            savefig(None, out)
            continue
        plot_per_sample(df_seg, widths, colors, out)

    # box plot: distribution of the matched fraction across samples, per max_width
    p = (
        ggplot(df, aes(x="max_width", y="fraction_matched", fill="max_width"))
        + geom_boxplot(outlier_shape="", alpha=0.6)
        + geom_jitter(width=0.15, height=0, size=1.5)
        + facet_wrap("method", ncol=n_methods)
        + theme_classic()
        + theme(panel_border=element_rect(fill=None, color="black"),
                text=element_text(font="Helvetica"),
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
    print(f"written: {args.out_table}, {args.out_barplot_multimodal}, {args.out_barplot_proseg}, {args.out_boxplot}")


if __name__ == "__main__":
    main()
