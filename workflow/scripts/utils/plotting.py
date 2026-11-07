"""Shared plotting setup and save helpers, so every figure of the pipeline looks the same (Helvetica).

The font file is shipped in utils/fonts/ because Helvetica is not installed in the container.
Usage: call set_mpl_rcparams() once before plotting, then save with savefig() / savefig_plotnine().
plotnine plots should still set text=element_text(font="Helvetica") in theme(), since plotnine
themes do not inherit font.family from rcParams.
"""

import logging
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager

script_dir = os.path.dirname(os.path.abspath(__file__))

# fixed colours of some tissues (first part of the sample ID); any other tissue gets a default palette colour
tissue_cols = {"CRC": "#5696E8",
               "DLBCL": "#AA3F5DFF",
               "GBM": "#938467",
               "NSCLC": "#CC79A7"}

# fixed colours of the TACIT cell types of the Xenium hImmune protein signatures (config/example/signatures/),
# by lineage hue family; any other cell type gets a default palette colour:
# lymphoid = blues (B / plasma indigo), myeloid = greens, stromal = browns, endothelial = gold, epithelial = purples,
# malignant = reds. The trailing comments name the matching RNA cell type labels.
TACIT_mm_cols = {# lymphoid
                 "T/NK cell": "#4A6FA5",                     # T_CD4, T_CD8, NK
                 "B cell": "#5C5A93",                        # B_cell
                 "Plasma cell": "#38366B",                   # Plasma
                 # myeloid
                 "Macrophage": "#4E8C72",                    # MoMac
                 "Dendritic cell": "#82B39A",                # DC
                 "Neutrophil": "#2E5E4E",                    # Granulocyte
                 # stromal
                 "Fibroblast": "#C08A4E",                    # Fibroblast
                 "Smooth muscle": "#8B5E3C",                 # Muscle
                 "Endothelial": "#D9B168",                   # Endothelial
                 # epithelial
                 "Epithelial": "#7E6BA8",                    # Epithelial_absorptive
                 # malignant
                 "Malignant epithelial": "#9E4B54",          # Malignant
                 "Cycling malignant epithelial": "#C4787E"}  # Malignant_crc


# Plotting functions for consistency
def set_mpl_rcparams(
    params={
        "pdf.fonttype": 42,  # embed TrueType fonts
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "text.usetex": False,
        "font.family": "Helvetica",
        "savefig.transparent": False,
        "figure.figsize": (4, 4),
        "axes.titlesize": 16,
    },
):
    # add Helvetica as font
    font_path = os.path.join(script_dir, "fonts/Helvetica.ttf")
    font_manager.fontManager.addfont(font_path)

    # apply the shared matplotlib defaults
    mpl.rcParams.update(params)


def savefig(df, out_file, dpi=300, bbox_inches="tight", close=True, pdf=False):
    # prevents printing tons of stuff when saving
    logging.getLogger("fontTools").setLevel(logging.ERROR)
    logging.getLogger("fontTools.subset").setLevel(logging.ERROR)
    logging.getLogger("matplotlib").setLevel(logging.ERROR)

    # create the output folder if needed
    out_file = Path(out_file)
    out_file.parent.mkdir(exist_ok=True, parents=True)
    # to ensure df is not forgotten by mistake
    if df is not None:
        df.to_csv(Path(out_file).with_suffix(".csv"))
    # always save a png, plus a pdf on request
    if pdf:
        plt.savefig(Path(out_file).with_suffix(".pdf"), dpi=dpi, bbox_inches=bbox_inches)
    plt.savefig(Path(out_file).with_suffix(".png"), dpi=dpi, bbox_inches=bbox_inches)
    if close:
        plt.close()


def _patch_plotnine_layout():
    """Let LayoutTree.resize honour per-composition panel ratios (set by panel_ratios), once.

    plotnine (0.15) has no plot_layout(widths=/heights=): after drawing, its layout tree resets the gridspec
    of every composition so that all panels get the same size, which overrides any ratio set beforehand.
    """
    from plotnine._mpl.layout_manager._layout_tree import ColumnsTree, LayoutTree, RowsTree
    if getattr(LayoutTree, "_panel_ratios_patched", False):
        return
    create = LayoutTree.create

    def create_with_ratios(cmp, lookup_spaces):
        # carry the composition's ratios over to its layout tree (nested compositions recurse through here)
        tree = create(cmp, lookup_spaces)
        tree.panel_ratios = getattr(cmp, "_panel_ratios", None)
        return tree

    def new_sizes(tree, plot_sizes, panel_sizes):
        # panels split their mean size by the ratios (equal without ratios), space around the panels unchanged
        plot_sizes, panel_sizes = np.array(plot_sizes), np.array(panel_sizes)
        ratios = getattr(tree, "panel_ratios", None)
        ratios = np.ones(len(panel_sizes)) if ratios is None else ratios
        return panel_sizes.mean() * ratios / ratios.mean() + (plot_sizes - panel_sizes)

    def resize_columns(self):
        widths = new_sizes(self, self.plot_widths, self.panel_widths)
        self.gridspec.set_width_ratios(widths / widths.min())
        self.resize_sub_compositions()

    def resize_rows(self):
        heights = new_sizes(self, self.plot_heights, self.panel_heights)
        self.gridspec.set_height_ratios(heights / heights.max())
        self.resize_sub_compositions()

    LayoutTree.create = staticmethod(create_with_ratios)
    ColumnsTree.resize = resize_columns
    RowsTree.resize = resize_rows
    LayoutTree._panel_ratios_patched = True


def panel_ratios(cmp, ratios):
    """Panel size ratios of a plotnine composition: widths for p1 | p2, heights for p1 / p2.

    Apply it last (after any `+`/`&` on the composition). Each plot is still drawn by plotnine, matplotlib's
    gridspec sizes them, e.g. savefig_plotnine(panel_ratios(p_frac | p_counts, [3, 1]), None, out).
    """
    _patch_plotnine_layout()
    cmp._panel_ratios = np.asarray(ratios, dtype=float)
    return cmp


# function to save plots with plotnine
def savefig_plotnine(plot, df, out_file, dpi=300, bbox_inches="tight", width=None, height=None, close=True, pdf=False):
    # prevents printing tons of stuff when saving
    logging.getLogger("fontTools").setLevel(logging.ERROR)
    logging.getLogger("fontTools.subset").setLevel(logging.ERROR)
    logging.getLogger("matplotlib").setLevel(logging.ERROR)

    # create the output folder if needed
    out_file = Path(out_file)
    out_file.parent.mkdir(exist_ok=True, parents=True)
    # to ensure df is not forgotten by mistake
    if df is not None:
        df.to_csv(Path(out_file).with_suffix(".csv"))

    # always save a png, plus a pdf on request
    if pdf:
        plot.save(Path(out_file).with_suffix(".pdf"), dpi=dpi, bbox_inches=bbox_inches, width=width, height=height)
    plot.save(Path(out_file).with_suffix(".png"), dpi=dpi, bbox_inches=bbox_inches, width=width, height=height)
    if close:
        plt.close()
