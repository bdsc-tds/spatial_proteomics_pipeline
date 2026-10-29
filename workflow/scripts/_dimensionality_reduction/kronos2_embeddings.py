"""KRONOS2 embedding of every cell of a mask: the CLS token of a box centred on the cell, holding only its pixels.

Follows CORAL's KRONOS2 cell recipe (github.com/mahmoodlab/CORAL), at the image's native resolution: a box_size px
box around the rounded centroid of each cell (zero-padded at the image border), every pixel outside the cell set to
0 (neighbouring cells too), divided by the maximum of the image dtype, z-scored per marker with KRONOS2's
pretraining statistics (model.preprocess) and passed through KRONOS2 (Hugging Face MahmoodLab/KRONOS2, loaded from
model_dir); the embedding is the 768-d CLS token. Channels are those with a KRONOS2 marker name in the panel's
marker table, with their pixel value offset removed (+100 on Xenium protein channels).

The image is read one horizontal band of cells at a time (all channels of the band's rows), so the full-resolution
stack is never in memory.

Writes cell_embeddings.npz: cell_id (mask label, or the 10x cell_id of --cell_ids) and embedding (cells x 768,
float32), in mask label order.
"""

import argparse
import os
import sys
import tempfile
import time

import numpy as np
import pandas as pd
import tifffile
import torch
import zarr

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(script_dir, "../utils/"))
from readwrite import MultiplexImage


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--image", required=True, help="multiplexed image (qptiff or OME-TIFF)")
    parser.add_argument("--markers", required=True, help="marker table of the panel (channel_name, kronos2_marker)")
    parser.add_argument("--mask", required=True, help="cell label mask on the image's pixel grid (tiled TIFF)")
    parser.add_argument("--cell_ids", default="", help="label -> cell_id table (CSV); empty = cell_id is the label")
    parser.add_argument("--model_dir", required=True, help="local copy of Hugging Face MahmoodLab/KRONOS2")
    parser.add_argument("--box_size", type=int, default=64, help="px box around each cell centroid (CORAL: 64)")
    parser.add_argument("--band_rows", type=int, default=2048, help="image rows of cell centroids read at once")
    parser.add_argument("--batch_size", type=int, default=256)
    parser.add_argument("--out_npz", required=True)
    return parser.parse_args()


def load_model(model_dir, device):
    """KRONOS2 in eval mode, from its local copy, without network access or a Hugging Face token."""
    # trust_remote_code copies the model code to $HF_HOME/modules: a private, writable folder per job
    os.environ["HF_HOME"] = tempfile.mkdtemp(prefix="kronos2_hf_")
    os.environ["HF_HUB_OFFLINE"] = "1"
    from transformers import AutoModel

    return AutoModel.from_pretrained(model_dir, trust_remote_code=True, device=device)


def cell_centroids(mask, chunk_rows=2048):
    """Labels present in the mask, their rounded centroid (row, col) and pixel count, accumulated by row chunks."""
    n = int(mask.max()) + 1
    count, sum_r, sum_c = np.zeros(n), np.zeros(n), np.zeros(n)
    cols = np.arange(mask.shape[1], dtype=np.float64)
    for r0 in range(0, mask.shape[0], chunk_rows):
        # pixel count and coordinate sums of every label in this chunk of rows
        m = mask[r0:r0 + chunk_rows]
        lab = m.ravel()
        count += np.bincount(lab, minlength=n)
        sum_r += np.bincount(lab, weights=np.repeat(np.arange(r0, r0 + m.shape[0], dtype=np.float64), m.shape[1]),
                             minlength=n)
        sum_c += np.bincount(lab, weights=np.tile(cols, m.shape[0]), minlength=n)
    labels = np.nonzero(count[1:])[0] + 1
    rows = np.rint(sum_r[labels] / count[labels]).astype(np.int64)
    cols = np.rint(sum_c[labels] / count[labels]).astype(np.int64)
    return labels, rows, cols, count[labels].astype(np.int64)


def read_rows(img, channels, r0, r1):
    """(channels, r1 - r0, X) block of rows of the given channels, offsets removed, read from the tiles it covers."""
    block = None
    for k, name in enumerate(channels):
        c = img.channel_index(name)
        # only the tiles of rows r0:r1 are decoded (the store is not closed: it shares the image's file handle)
        plane = zarr.open(img._level.pages[c].aszarr(), mode="r")[r0:r1]
        offset = img.channel_offsets[c]
        if offset:
            np.maximum(plane, offset, out=plane)
            plane -= np.asarray(offset, dtype=plane.dtype)
        if block is None:
            block = np.empty((len(channels), *plane.shape), dtype=plane.dtype)
        block[k] = plane
    return block


def main():
    args = parse_args()
    device = "cuda" if torch.cuda.is_available() else "cpu"
    half = args.box_size // 2

    # channels KRONOS2 reads: those with a KRONOS2 marker name, in marker table order
    markers = pd.read_csv(args.markers, dtype=str)
    markers = markers[markers["kronos2_marker"].notna()].reset_index(drop=True)
    channels, names = markers["channel_name"].tolist(), markers["kronos2_marker"].str.strip().tolist()
    preferred_dapi = "DAPI" if "DAPI" in names else None
    print(f"{len(names)} channels with a KRONOS2 marker in {args.markers}: "
          + ", ".join(f"{c} -> {n}" for c, n in zip(channels, names)), flush=True)

    # model first: an unknown marker name fails before any image is read
    model = load_model(args.model_dir, device)
    model(torch.zeros(1, len(names), args.box_size, args.box_size, device=device), names)
    print(f"loaded KRONOS2 from {args.model_dir} on {device}", flush=True)

    with MultiplexImage(args.image) as img:
        # full-resolution mask, on the image's pixel grid
        height, width = img.shape[1:]
        mask = tifffile.imread(args.mask)
        if mask.shape != (height, width):
            raise ValueError(f"mask {mask.shape} is not on the image's pixel grid {(height, width)}")
        scale = np.iinfo(img.dtype).max if np.issubdtype(img.dtype, np.integer) else 1.0

        # rounded centroid of every cell
        labels, rows, cols, n_pixels = cell_centroids(mask)
        print(f"{len(labels)} cells, {img.dtype} image {(height, width)} (divided by {scale})", flush=True)

        embeddings = np.zeros((len(labels), model.config.embed_dim), dtype=np.float32)
        in_box = np.zeros(len(labels), dtype=np.int64)
        t0 = time.time()
        for b0 in range(0, height, args.band_rows):
            # cells whose centroid lies in this band of rows
            cells = np.nonzero((rows >= b0) & (rows < b0 + args.band_rows))[0]
            if len(cells) == 0:
                continue

            # the band's image rows and mask rows, with every box inside it (zero-padded at the image border)
            top = b0 - half
            r0, r1 = max(top, 0), min(b0 + args.band_rows - half + args.box_size, height)
            band = np.zeros((len(channels), args.band_rows + args.box_size, width + args.box_size), dtype=img.dtype)
            band[:, r0 - top:r1 - top, half:half + width] = read_rows(img, channels, r0, r1)
            band_mask = np.zeros(band.shape[1:], dtype=mask.dtype)
            band_mask[r0 - top:r1 - top, half:half + width] = mask[r0:r1]

            for s in range(0, len(cells), args.batch_size):
                # box of every cell of the batch, pixels outside the cell set to 0, divided by the dtype maximum
                batch = cells[s:s + args.batch_size]
                ys, xs = rows[batch] - top - half, cols[batch]
                crops = np.stack([band[:, y:y + args.box_size, x:x + args.box_size] for y, x in zip(ys, xs)])
                inside = np.stack([band_mask[y:y + args.box_size, x:x + args.box_size] for y, x in zip(ys, xs)])
                inside = inside == labels[batch][:, None, None]
                in_box[batch] = inside.sum(axis=(1, 2))
                patches = crops.astype(np.float32) * inside[:, None] / scale

                # per-marker z-score with the pretraining statistics, then the CLS token
                patches = model.preprocess(patches, names, preferred_dapi=preferred_dapi)
                cls = model(torch.from_numpy(patches).to(device), names)
                embeddings[batch] = cls.float().cpu().numpy()
            print(f"rows {b0}-{min(b0 + args.band_rows, height)} of {height}: {len(cells)} cells "
                  f"({time.time() - t0:.0f} s)", flush=True)
            del band, band_mask

    # part of each cell inside its box (cells larger than the box lose their edges)
    frac = in_box / n_pixels
    print(f"pixels inside the {args.box_size} px box: median {np.median(frac):.3f}, "
          f"{np.mean(frac < 0.95) * 100:.1f}% of cells below 0.95, {np.mean(frac < 0.8) * 100:.1f}% below 0.8",
          flush=True)

    # cell ids: mask labels, or the 10x cell_id of every label
    cell_id = labels
    if args.cell_ids:
        table = pd.read_csv(args.cell_ids).set_index("label")["cell_id"]
        cell_id = table.reindex(labels).fillna(pd.Series(labels.astype(str), index=labels)).to_numpy(str)
        print(f"{np.isin(labels, table.index).sum()} of {len(labels)} labels in {args.cell_ids}", flush=True)

    np.savez(args.out_npz, cell_id=cell_id, embedding=embeddings)
    print(f"written: {args.out_npz} ({len(labels)} cells x {embeddings.shape[1]})")


if __name__ == "__main__":
    main()
