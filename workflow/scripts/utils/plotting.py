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
