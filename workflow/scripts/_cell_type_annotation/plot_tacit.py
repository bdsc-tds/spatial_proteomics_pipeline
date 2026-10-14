"""Plot the TACIT cell types of one segmentation method over all samples.

Reads, per sample (same order in every list): the feature table of quantify_cells, the TACIT labels
of run_tacit.R and the signature CSV of the sample's panel.
Writes:
  - a stacked bar plot of the fraction of each cell type per sample
  - per sample, a heatmap of the mean normalized intensity (arcsinh(x / cofactor), z-score within
    the sample) of every signature marker per cell type
  - per sample, the UMAP of preprocess_protein (same cofactor) coloured by cell type, once on the
    signature markers and once on all markers except --exclude
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
    parser.add_argument("--exclude", nargs="*", default=["DAPI"], help="channels left out of the all-markers UMAP")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--out_fractions", required=True, help="bar plot of the cell type fractions")
    parser.add_argument("--out_heatmaps", nargs="+", required=True, help="one heatmap per sample")
    parser.add_argument("--out_umaps", nargs="+", required=True, help="one UMAP figure per sample")
    args = parser.parse_args()

    # every per-sample list must follow --sample_ids
    n = len(args.sample_ids)
    for name in ["feature_tables", "cell_types", "signatures", "out_heatmaps", "out_umaps"]:
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
from plotnine import (aes, coord_flip, element_rect, element_text, geom_col, ggplot, labs, scale_fill_manual,
                      scale_y_continuous, theme, theme_classic)

# PyComplexHeatmap still calls the removed matplotlib.cm.get_cmap in some versions
import matplotlib.cm
if not hasattr(matplotlib.cm, "get_cmap"):
    matplotlib.cm.get_cmap = matplotlib.colormaps.get_cmap
import PyComplexHeatmap as pch

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(script_dir, "../utils/"))
from multimodal_readwrite import mesmer_feature_table_to_adata
from plotting import savefig, savefig_plotnine, set_mpl_rcparams
from preprocessing import preprocess_protein

OTHERS = "Others"
OTHERS_COLOR = "#BDBDBD"


def sample_label(sample_id):
    """Short sample name for the plots: donor/sample."""
    return "/".join(sample_id.split("/")[2:])


def cell_type_palette(signatures):
    """Cell type order (signature order, then Others) and one colour per cell type."""
    # union of the cell types of all signatures, in order of appearance
    order = list(dict.fromkeys(ct for sig in signatures for ct in sig.index))
    colors = sns.color_palette("tab20" if len(order) <= 20 else "husl", len(order)).as_hex()
    palette = dict(zip(order, colors))
    palette[OTHERS] = OTHERS_COLOR
    return order + [OTHERS], palette


def plot_fractions(labels, ct_order, palette, out):
    """Stacked bar of the fraction of each cell type, one bar per sample."""
    # cells and fraction per sample and cell type
    counts = (pd.concat([lab.assign(sample=sample_label(s)) for s, lab in labels.items()])
              .groupby(["sample", "cell_type"]).size().rename("cells").reset_index())
    counts["fraction"] = counts["cells"] / counts.groupby("sample")["cells"].transform("sum")

    # sample order as given; cell types stacked in signature order
    counts["sample"] = pd.Categorical(counts["sample"], categories=[sample_label(s) for s in labels][::-1])
    counts["cell_type"] = pd.Categorical(counts["cell_type"], categories=[ct for ct in ct_order[::-1]
                                                                          if ct in set(counts["cell_type"])])

    p = (
        ggplot(counts, aes("sample", "fraction", fill="cell_type"))
        + geom_col(width=0.8)
        + coord_flip()
        + scale_fill_manual(values=palette, breaks=[ct for ct in ct_order if ct in set(counts["cell_type"])])
        + labs(x="", y="fraction of cells", fill="TACIT cell type")
        + theme_classic()
        + theme(panel_border=element_rect(fill=None, color="black"),
                text=element_text(font="Helvetica"),
                figure_size=(8, 0.4 * len(labels) + 2))
        + scale_y_continuous(expand=[0, 0, 0.01, 0])
    )
    savefig_plotnine(p, None, out)


def normalize(features, markers, cofactor):
    """arcsinh(x / cofactor) and z-score per marker within the sample, as run_tacit.R (without its shift to min 0)."""
    norm = np.arcsinh(features[markers] / cofactor)
    return (norm - norm.mean()) / norm.std(ddof=1)


def plot_heatmap(features, cell_types, markers, cofactor, ct_order, title, out):
    """Mean normalized intensity of every signature marker per cell type, cell numbers on the right."""
    # mean z-score per cell type and marker, rows in signature order
    norm = normalize(features, markers, cofactor)
    means = norm.groupby(cell_types.to_numpy()).mean()
    means = means.loc[[ct for ct in ct_order if ct in means.index]]

    # number of cells per cell type as a row annotation
    n_cells = cell_types.value_counts().loc[means.index].to_frame("cells")
    right = pch.HeatmapAnnotation(cells=pch.anno_barplot(n_cells, colors="grey", legend=False), axis=0, verbose=0)

    plt.figure(figsize=(0.3 * len(markers) + 3, 0.3 * len(means) + 2))
    pch.ClusterMapPlotter(data=means, row_cluster=False, col_cluster=True, col_dendrogram=True,
                          right_annotation=right, show_rownames=True, show_colnames=True,
                          row_names_side="left", cmap="RdBu_r", center=0, vmin=-2, vmax=2,
                          label="mean z-score", xticklabels_kws={"labelrotation": 90}, verbose=0)
    plt.suptitle(title, y=1.02)
    savefig(None, out)


def umap_coordinates(feature_table, markers, cofactor, exclude):
    """preprocess_protein UMAP of one sample, on the signature markers and on all markers but exclude."""
    adata = mesmer_feature_table_to_adata(feature_table, dapi_to_obs=False)
    coords = pd.DataFrame({"cell_id": adata.obs["cell_id"].to_numpy()})

    # one UMAP per marker set, each from the raw intensities
    for name, include in [("signature", markers), ("all", None)]:
        print(f"preprocess_protein: {name} markers", flush=True)
        processed = preprocess_protein(adata.copy(), arcsinh_cofactor=cofactor, include=include, exclude=exclude)
        coords[[f"UMAP1_{name}", f"UMAP2_{name}"]] = processed.obsm["X_umap"]
    return coords


def plot_umaps(coords, cell_types, ct_order, palette, title, out):
    """UMAPs of one sample side by side, coloured by TACIT cell type."""
    df = coords.assign(cell_type=cell_types.to_numpy())
    hue_order = [ct for ct in ct_order if ct in set(df["cell_type"])]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, name in zip(axes, ["signature", "all"]):
        sns.scatterplot(data=df, x=f"UMAP1_{name}", y=f"UMAP2_{name}", hue="cell_type", hue_order=hue_order,
                        palette=palette, s=0.5, linewidth=0, rasterized=True,
                        legend="auto" if ax is axes[-1] else False, ax=ax)
        ax.set(title=f"{name} markers", xlabel="UMAP1", ylabel="UMAP2", xticks=[], yticks=[])
    axes[-1].legend(bbox_to_anchor=(1.01, 1), loc="upper left", markerscale=8, frameon=False, title="TACIT cell type")
    fig.suptitle(title)
    savefig(None, out)


def main(args):
    set_mpl_rcparams()

    # signatures (cell_type x marker) and TACIT labels of every sample
    signatures = [pd.read_csv(path).set_index("cell_type") for path in args.signatures]
    labels = {s: pd.read_csv(path) for s, path in zip(args.sample_ids, args.cell_types)}
    ct_order, palette = cell_type_palette(signatures)

    # fraction of each cell type per sample
    plot_fractions(labels, ct_order, palette, args.out_fractions)

    for i, sample_id in enumerate(args.sample_ids):
        print(f"{sample_id} ({i + 1}/{len(args.sample_ids)})", flush=True)
        features = pd.read_csv(args.feature_tables[i])
        markers = list(signatures[i].columns)

        # TACIT labels in feature table order (run_tacit.R keeps it)
        cell_types = features[["cell_id"]].merge(labels[sample_id], on="cell_id", how="left", validate="1:1")["cell_type"]
        if cell_types.isna().any():
            raise ValueError(f"{sample_id}: {cell_types.isna().sum()} cells without a TACIT label")

        # mean marker intensity per cell type
        plot_heatmap(features, cell_types, markers, args.cofactor, ct_order,
                     f"{sample_label(sample_id)}: mean marker intensity per TACIT cell type", args.out_heatmaps[i])

        # preprocess_protein UMAPs coloured by cell type
        coords = umap_coordinates(args.feature_tables[i], markers, args.cofactor, args.exclude)
        plot_umaps(coords, cell_types, ct_order, palette, f"{sample_label(sample_id)} ({len(coords)} cells)",
                   args.out_umaps[i])


if __name__ == "__main__":
    main(ARGS)
