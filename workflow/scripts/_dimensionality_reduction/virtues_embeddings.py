"""VirTues embedding of every cell of a mask: the overlap-weighted mean of the backbone's patch summary tokens.

Follows VirTues' cell tokens (virtues.utils.cell_tokens.compute_cell_tokens, notebooks/2_demo_cell_phenotyping):
the 1 µm stack, zero-padded by 120 px and then standardized and blurred (virtues_utils.standardized_stack), is cut
into 128 px crops at stride 42; the backbone (Hugging Face bunnelab/virtues) gives one 512-d patch summary token
per 8 x 8 px patch of every crop, and a cell's embedding is the mean of the tokens of every (crop, patch) it
overlaps, weighted by its number of 1 µm pixels in that patch. The mask is brought to the 1 µm grid by centre
sampling (each 1 µm pixel takes the label of the full-resolution pixel at its centre) and zero-padded like the
stack; cells covering no 1 µm pixel centre get no embedding.

The only difference to compute_cell_tokens, which cannot run on a whole slide (it moves every crop to the GPU at
once and keeps one token reference per (cell, patch, crop) in CPU RAM): the weighted sums are accumulated one
batch of crops at a time, in (cells, 512) arrays on the GPU.

Writes cell_embeddings.npz: cell_id (mask label, or the 10x cell_id of --cell_ids) and embedding (cells x 512,
float32), in mask label order.
"""

import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
import tifffile
import torch

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(script_dir, "../utils/"))
from readwrite import MultiplexImage
from virtues_utils import channel_ids, load_backbone, standardized_stack


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--image_1um", required=True, help="(channels, Y, X) stack of virtues_prepare_image")
    parser.add_argument("--channel_stats", required=True, help="channel table of the stack (virtues_prepare_image)")
    parser.add_argument("--standardization", required=True, help="panel mean / std (virtues_standardization_stats)")
    parser.add_argument("--embeddings", required=True, help="marker embeddings of the panel (virtues_marker_embeddings)")
    parser.add_argument("--image", required=True, help="full-resolution image, for the mask's pixel grid")
    parser.add_argument("--mask", required=True, help="cell label mask on the image's pixel grid (tiled TIFF)")
    parser.add_argument("--cell_ids", default="", help="label -> cell_id table (CSV); empty = cell_id is the label")
    parser.add_argument("--model_dir", required=True, help="folder holding (or receiving) the VirTues weights")
    parser.add_argument("--model", default="virtues-sp32", help="model folder of bunnelab/virtues on Hugging Face")
    parser.add_argument("--blur_sigma", type=float, default=1.0, help="std of the 3 x 3 Gaussian blur (0 = none)")
    parser.add_argument("--pad", type=int, default=120, help="zero padding of the 1 µm stack and mask (px, VirTues demo)")
    parser.add_argument("--tile_size", type=int, default=128)
    parser.add_argument("--stride", type=int, default=42)
    parser.add_argument("--batch_size", type=int, default=32)
    parser.add_argument("--out_npz", required=True)
    return parser.parse_args()


def crop_starts(size, tile_size, stride):
    """Start of every crop along one axis: a regular grid plus a last crop flush with the end (as VirTues)."""
    starts = list(range(0, size - tile_size + 1, stride))
    if starts[-1] != size - tile_size:
        starts.append(size - tile_size)
    return starts


def centre_sample(mask, height_1um, width_1um):
    """Label mask on the 1 µm grid: each 1 µm pixel takes the label of the full-resolution pixel at its centre."""
    rows = ((np.arange(height_1um) + 0.5) * mask.shape[0] / height_1um).astype(np.int64)
    cols = ((np.arange(width_1um) + 0.5) * mask.shape[1] / width_1um).astype(np.int64)
    return mask[np.ix_(rows, cols)]


def cell_embeddings(backbone, x, ids, mask, tile_size=128, patch_size=8, stride=42, batch_size=32, device="cuda"):
    """Overlap-weighted mean patch summary token of every label of mask.

    x: (C, Y, X) standardized 1 µm stack; ids: its channel identities; mask: (Y, X) label mask on the same grid.
    Returns the labels present in the mask and their embeddings (labels x D, float32), on the CPU.
    """
    # every crop of the stack, and the patch index of every pixel of a crop
    positions = [(i, j) for i in crop_starts(x.shape[1], tile_size, stride)
                 for j in crop_starts(x.shape[2], tile_size, stride)]
    g = tile_size // patch_size
    patch = torch.arange(tile_size, device=device) // patch_size
    patch = patch[:, None] * g + patch[None, :]
    tokens_per_batch = batch_size * g * g
    ids = ids.to(device)

    sums, weights = None, None
    n_labels = int(mask.max()) + 1
    t0 = time.time()
    for b, start in enumerate(range(0, len(positions), batch_size)):
        # patch summary tokens of one batch of crops: (crops * g * g, D)
        batch = positions[start:start + batch_size]
        crops = [torch.from_numpy(x[:, i:i + tile_size, j:j + tile_size].copy()).to(device) for i, j in batch]
        with torch.no_grad(), torch.amp.autocast("cuda", enabled=device != "cpu"):
            out = backbone.encoder.forward_list(crops, [ids] * len(crops), multiplex_mask=None)
        tokens = torch.stack(out.patch_summary_tokens).reshape(len(batch) * g * g, -1).float()
        if sums is None:
            sums = torch.zeros(n_labels, tokens.shape[1], dtype=torch.float32, device=device)
            weights = torch.zeros(n_labels, dtype=torch.float32, device=device)

        # key of every pixel of every crop: label * tokens_per_batch + token index (crop, patch)
        labels = torch.from_numpy(np.stack([mask[i:i + tile_size, j:j + tile_size] for i, j in batch])).to(device)
        token = torch.arange(len(batch), device=device)[:, None, None] * g * g + patch[None]
        keys = labels.long() * tokens_per_batch + token

        # pixels of every (cell, token) pair, added to the cell's weighted sum
        pairs, counts = torch.unique(keys, return_counts=True)
        cells, token = pairs // tokens_per_batch, pairs % tokens_per_batch
        keep = cells > 0
        cells, token, counts = cells[keep], token[keep], counts[keep].float()
        sums.index_add_(0, cells, counts[:, None] * tokens[token])
        weights.index_add_(0, cells, counts)

        if b % 200 == 0:
            print(f"crops {start + len(batch)}/{len(positions)} ({time.time() - t0:.0f} s)", flush=True)

    # weighted mean of every label present in the mask
    labels = torch.nonzero(weights > 0).squeeze(1)
    embeddings = sums[labels] / weights[labels, None]
    return labels.cpu().numpy(), embeddings.cpu().numpy()


def main():
    args = parse_args()

    # standardized, blurred, padded stack and the channel identities (row of each channel's protein in the table)
    x, channels = standardized_stack(args.image_1um, args.channel_stats, args.standardization, args.blur_sigma,
                                     pad=args.pad)
    embeddings = torch.load(args.embeddings)
    ids = channel_ids(embeddings, channels)

    # full-resolution mask, on the image's pixel grid
    with MultiplexImage(args.image) as img:
        height, width = img.shape[1:]
    mask = tifffile.imread(args.mask)
    if mask.shape != (height, width):
        raise ValueError(f"mask {mask.shape} is not on the image's pixel grid {(height, width)}")
    n_cells = np.count_nonzero(np.bincount(mask.ravel())[1:])

    # 1 µm mask on the stack's grid (centre sampling), zero-padded like the stack
    mask_1um = np.pad(centre_sample(mask, x.shape[1] - 2 * args.pad, x.shape[2] - 2 * args.pad), args.pad)
    del mask
    n_cells_1um = np.count_nonzero(np.bincount(mask_1um.ravel())[1:])
    print(f"mask {(height, width)} -> {mask_1um.shape} at 1 µm (padded): {n_cells_1um} of {n_cells} cells kept",
          flush=True)

    # cell embeddings, batch of crops by batch of crops on the GPU
    backbone = load_backbone(embeddings["embeddings"], args.model_dir, args.model)
    print(f"loaded {args.model} backbone", flush=True)
    labels, cell_tokens = cell_embeddings(backbone, x, ids, mask_1um, tile_size=args.tile_size, stride=args.stride,
                                          batch_size=args.batch_size)

    # cell ids: mask labels, or the 10x cell_id of every label
    cell_id = labels
    if args.cell_ids:
        table = pd.read_csv(args.cell_ids).set_index("label")["cell_id"]
        cell_id = table.reindex(labels).fillna(pd.Series(labels.astype(str), index=labels)).to_numpy(str)
        print(f"{np.isin(labels, table.index).sum()} of {len(labels)} labels in {args.cell_ids}", flush=True)

    np.savez(args.out_npz, cell_id=cell_id, embedding=cell_tokens.astype(np.float32))
    print(f"written: {args.out_npz} ({len(labels)} cells x {cell_tokens.shape[1]})")


if __name__ == "__main__":
    main()
