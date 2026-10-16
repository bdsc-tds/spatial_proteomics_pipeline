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
from matplotlib import font_manager

script_dir = os.path.dirname(os.path.abspath(__file__))

# tissue (first part of the sample ID) colours: cols_tissue in analysis_repo's plotting_utils.py,
# NSCLC from its cols_cancer_type
tissue_cols = {"CRC": "#5696E8",
               "DLBCL": "#AA3F5DFF",
               "GBM": "#938467",
               "NSCLC": "#CC79A7"}

# TACIT cell types of the Xenium MM protein panel coloured as the RNA Level2 cell types (cols_level2_mm
# in analysis_repo's plotting_utils.py; lymphoid from Level3, cols_level3), same lineage hue families when
# there is no RNA equivalent:
# lymphoid = blues (B / plasma indigo), myeloid = greens, stromal = browns, endothelial = gold, epithelial = purples,
# malignant = reds
TACIT_mm_cols = {# lymphoid (cols_level3)
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
