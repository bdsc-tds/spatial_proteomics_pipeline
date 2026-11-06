"""VirTues cell type of every cell of the VirTues mask, from the segmentation head's cell type logits.

Follows VirTues' segmentation demo (notebooks/3_demo_segmentation.ipynb), except that only the None
(background) logit is suppressed and Unknown is kept as a cell type: every 1 µm pixel takes its most likely
remaining class, and every cell takes the majority class of its pixels
(virtues.utils.segmentation.assign_cell_types). The mask is brought to the logits' 1 µm grid as in
virtues_embeddings.py, by centre sampling (each 1 µm pixel takes the label of the full-resolution pixel at its
centre); cells covering no 1 µm pixel centre get no cell type.

Writes one row per cell: cell_id (mask label), cell_type, confidence (pixel fraction of that class) and the pixel
fraction of every class.
"""

import argparse

import numpy as np
import pandas as pd
import tifffile

# classes suppressed before the argmax: the background class only (the demo also suppresses Unknown)
SUPPRESSED = ["None"]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--logits", required=True, help="cell type logits (classes x Y x X at 1 µm, segmentation_virtues)")
    parser.add_argument("--classes", required=True, help="class index -> cell type (classes.csv)")
    parser.add_argument("--mask", required=True, help="cell label mask on the image's pixel grid (tiled TIFF)")
    parser.add_argument("--out_csv", required=True)
    return parser.parse_args()


def centre_sample(mask, height_1um, width_1um):
    """Label mask on the 1 µm grid: each 1 µm pixel takes the label of the full-resolution pixel at its centre."""
    rows = ((np.arange(height_1um) + 0.5) * mask.shape[0] / height_1um).astype(np.int64)
    cols = ((np.arange(width_1um) + 0.5) * mask.shape[1] / width_1um).astype(np.int64)
    return mask[np.ix_(rows, cols)]


def main():
    args = parse_args()

    # logits and the classes kept (None suppressed); keep_default_na=False keeps "None" a string
    logits = tifffile.imread(args.logits)
    classes = pd.read_csv(args.classes, keep_default_na=False).sort_values("class_index")
    keep = classes.loc[~classes["cell_type"].isin(SUPPRESSED), "class_index"].to_numpy()
    names = classes.set_index("class_index").loc[keep, "cell_type"].to_numpy()

    # most likely kept class of every 1 µm pixel (index into keep, = argmax with the suppressed logits at -inf)
    pixel_class = logits[keep].argmax(axis=0)
    del logits

    # 1 µm mask on the logits' grid (centre sampling)
    mask = tifffile.imread(args.mask)
    n_cells = np.count_nonzero(np.bincount(mask.ravel())[1:])
    mask_1um = centre_sample(mask, *pixel_class.shape).astype(np.int64)
    del mask

    # pixels of every (label, class) pair, as assign_cell_types
    counts = np.bincount((mask_1um * len(keep) + pixel_class).ravel(), minlength=(mask_1um.max() + 1) * len(keep))
    counts = counts.reshape(-1, len(keep))

    # majority class per cell present at 1 µm, its pixel fraction as confidence
    cells = np.flatnonzero(counts[1:].sum(axis=1) > 0) + 1
    fractions = counts[cells] / counts[cells].sum(axis=1, keepdims=True)
    print(f"mask {mask_1um.shape} at 1 µm: {len(cells)} of {n_cells} cells kept", flush=True)
    table = pd.DataFrame({
        "cell_id": cells,
        "cell_type": names[fractions.argmax(axis=1)],
        "confidence": fractions.max(axis=1),
    })
    table = pd.concat([table, pd.DataFrame(fractions, columns=[f"fraction_{n}" for n in names])], axis=1)

    table.to_csv(args.out_csv, index=False)
    print(table["cell_type"].value_counts().to_string())
    print(f"written: {args.out_csv} ({len(table)} cells)")


if __name__ == "__main__":
    main()
