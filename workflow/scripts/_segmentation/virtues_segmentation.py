"""VirTues whole-cell segmentation, from the image prepared by virtues_prepare_image.py.

Standardizes and blurs the 1 µm stack (virtues_utils.standardized_stack), runs the VirTues backbone with its
segmentation head (Hugging Face bunnelab/virtues) tile by tile, and writes:
- the whole-cell mask on the full-resolution pixel grid of the image (nearest neighbour from 1 µm), as a tiled
  TIFF like the other methods' masks (VirTues predicts no nuclei);
- the head's cell type logits of every 1 µm pixel (classes x Y x X, float16, tiled TIFF), all classes kept
  (virtues_cell_types.py keeps them all too), and the class names (classes.csv).
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import tifffile
import torch

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_dir)
sys.path.insert(0, os.path.join(script_dir, "../utils/"))
from label_masks import write_label_mask
from readwrite import MultiplexImage
from virtues_utils import channel_ids, load_backbone, standardized_stack

# cell type classes of the segmentation head, in logit order (virtues notebooks/3_demo_segmentation.ipynb)
HEAD_CLASSES = ["None", "B cell", "CD4 T cell", "CD8 T cell", "Fibroblasts /Stroma", "Myeloid", "Tumor", "Unknown",
                "Vessel / Endo."]
# instance channels of the head, and fragments dropped from the instance mask (1 µm px, VirTues demo)
HEAD_DIM_OUT = 5
MIN_CELL_SIZE = 15


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--image_1um", required=True, help="(channels, Y, X) stack from virtues_prepare_image.py")
    parser.add_argument("--channel_stats", required=True, help="channel_stats.csv of the sample (channel order)")
    parser.add_argument("--standardization", required=True, help="mean / std of every channel of the panel")
    parser.add_argument("--embeddings", required=True, help="marker embeddings of the panel (.pt)")
    parser.add_argument("--image", required=True, help="full-resolution image, for the mask's pixel grid")
    parser.add_argument("--model_dir", required=True, help="folder holding (or receiving) the VirTues weights")
    parser.add_argument("--model", default="virtues-sp32", help="model folder of bunnelab/virtues on Hugging Face")
    parser.add_argument("--blur_sigma", type=float, default=1.0, help="std of the 3 x 3 Gaussian blur (0 = none)")
    parser.add_argument("--tile_size", type=int, default=128)
    parser.add_argument("--overlap", type=int, default=32)
    parser.add_argument("--batch_size", type=int, default=16)
    parser.add_argument("--out_cells", required=True)
    parser.add_argument("--out_logits", required=True)
    parser.add_argument("--out_classes", required=True)
    return parser.parse_args()


def load_model(embeddings, model_dir, model):
    """VirTues backbone with its segmentation head on the GPU; the weights are downloaded on first use."""
    from huggingface_hub import hf_hub_download
    from safetensors.torch import load_file
    from virtues.modules.segmentation.unet import VirtuesSegmentationHead

    backbone = load_backbone(embeddings, model_dir, model)

    # head: its file holds no backbone weights, so only virtues_model.* keys may be missing
    head = VirtuesSegmentationHead(virtues_model=backbone, dim_out=HEAD_DIM_OUT, num_celltypes=len(HEAD_CLASSES))
    state = load_file(hf_hub_download("bunnelab/virtues", f"{model}/segmentation-head.safetensors", local_dir=model_dir))
    missing, unexpected = head.load_state_dict(state, strict=False)

    # aliases of loaded tensors are not missing: instanseg shares pixel_classifier with instance_processor
    head_state = head.state_dict()
    loaded = {head_state[k].data_ptr() for k in state if k in head_state}
    missing = [k for k in missing if not k.startswith("virtues_model.") and head_state[k].data_ptr() not in loaded]
    if missing or unexpected:
        sys.exit(f"segmentation head does not match the model: missing {missing}, unexpected {unexpected}")
    print(f"loaded {model} backbone and segmentation head", flush=True)
    return head.cuda().eval()


def main():
    args = parse_args()
    from virtues.utils.segmentation import remove_small_cells

    # standardized, blurred stack and the channel identities (row of each channel's protein in the embedding table)
    x, channels = standardized_stack(args.image_1um, args.channel_stats, args.standardization, args.blur_sigma)
    embeddings = torch.load(args.embeddings)
    ids = channel_ids(embeddings, channels)

    # instance mask and cell type logits, tile by tile (VirTues stitches tissues larger than 1500 px itself)
    head = load_model(embeddings["embeddings"], args.model_dir, args.model)
    with torch.no_grad():
        instances, logits = head.segment_tissue(torch.from_numpy(x).cuda(), ids.cuda(), tile_size=args.tile_size,
                                                overlap=args.overlap, batch_size=args.batch_size)
    del x

    # cell type logits of every 1 µm pixel (all classes) and the class names
    tifffile.imwrite(args.out_logits, logits.cpu().numpy().astype(np.float16), tile=(1024, 1024), compression="zlib")
    del logits
    pd.DataFrame({"class_index": range(len(HEAD_CLASSES)), "cell_type": HEAD_CLASSES}).to_csv(args.out_classes,
                                                                                            index=False)

    # whole cells at 1 µm, small fragments dropped
    cells = remove_small_cells(instances.cpu().numpy(), min_cell_size=MIN_CELL_SIZE)

    # whole-cell mask on the image's pixel grid: each pixel takes the 1 µm pixel it falls in
    with MultiplexImage(args.image) as img:
        height, width = img.shape[1:]
    rows = np.minimum((np.arange(height) * cells.shape[0] / height).astype(np.int64), cells.shape[0] - 1)
    cols = np.minimum((np.arange(width) * cells.shape[1] / width).astype(np.int64), cells.shape[1] - 1)
    write_label_mask(args.out_cells, cells[rows][:, cols])


if __name__ == "__main__":
    main()
