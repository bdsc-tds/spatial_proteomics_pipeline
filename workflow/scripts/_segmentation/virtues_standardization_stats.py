"""Mean and standard deviation of every VirTues channel over the tissue pixels of all samples of a panel.

Pools the per-sample pixel counts, sums and sums of squares written by virtues_prepare_image.py (spora's
"global" statistics, on the clipped and log1p-transformed values). Writes channel_name, uniprot_id, mean, std.
"""

import argparse

import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--channel_stats", nargs="+", required=True, help="channel_stats.csv of every sample of the panel")
    parser.add_argument("--out_csv", required=True)
    return parser.parse_args()


def main():
    args = parse_args()

    # every sample of a panel must have the same channels, in the same order
    tables = [pd.read_csv(p) for p in args.channel_stats]
    channels = tables[0]["channel_name"].tolist()
    for path, table in zip(args.channel_stats, tables):
        if table["channel_name"].tolist() != channels:
            raise ValueError(f"{path}: channels differ from {args.channel_stats[0]}")

    # pooled moments over all tissue pixels of all samples
    stats = pd.concat(tables).groupby(["channel_name", "uniprot_id"], sort=False)[["n_pixels", "sum", "sum_sq"]].sum()
    mean = stats["sum"] / stats["n_pixels"]
    std = np.sqrt(np.maximum(stats["sum_sq"] / stats["n_pixels"] - mean**2, 0))
    out = pd.DataFrame({"mean": mean, "std": std}).reset_index()

    # a constant channel cannot be standardized
    if (out["std"] <= 0).any():
        raise ValueError(f"channels with zero variance: {out.loc[out['std'] <= 0, 'channel_name'].tolist()}")

    out.to_csv(args.out_csv, index=False)
    print(out.to_string(index=False))
    print(f"written: {args.out_csv} ({len(tables)} samples)")


if __name__ == "__main__":
    main()
