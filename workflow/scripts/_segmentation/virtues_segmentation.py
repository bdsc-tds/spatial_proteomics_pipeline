"""VirTues whole-cell segmentation and per-pixel cell types, from the image prepared by virtues_prepare_image.py.

Standardizes the 1 µm stack with the panel's mean / std (spora standardization, second half), runs the
VirTues backbone with its segmentation head (Hugging Face bunnelab/virtues) tile by tile, and writes:
- the whole-cell mask on the full-resolution pixel grid of the image (nearest neighbour from 1 µm), as a
  tiled TIFF like the other methods' masks (VirTues predicts no nuclei);
- the class map: the most likely class of every 1 µm pixel (head class index, None = background and
  Unknown included), as a tiled TIFF, and the class names (classes.csv).
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

# classes of the segmentation head, in output order (virtues notebooks/3_demo_segmentation.ipynb)
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
    parser.add_argument("--tile_size", type=int, default=128)
    parser.add_argument("--overlap", type=int, default=32)
    parser.add_argument("--batch_size", type=int, default=16)
    parser.add_argument("--out_cells", required=True)
    parser.add_argument("--out_class_map", required=True)
    parser.add_argument("--out_classes", required=True)
    return parser.parse_args()


def load_model(embeddings, model_dir, model):
    """VirTues backbone with its segmentation head on the GPU; the weights are downloaded on first use."""
    from huggingface_hub import hf_hub_download
    from safetensors.torch import load_file
    from virtues.modules.multiplex_virtues import MultiplexVirtues
    from virtues.modules.segmentation.unet import VirtuesSegmentationHead

    # backbone: the marker embedding table is ours (a non-persistent buffer, not part of the weights)
    backbone = MultiplexVirtues(prior_bias_embeddings=embeddings)
    backbone.load_state_dict(load_file(hf_hub_download("bunnelab/virtues", f"{model}/model.safetensors", local_dir=model_dir)))
    backbone = backbone.cuda().eval()

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

    # standardized stack: (log1p of the clipped values - panel mean) / panel std
    channels = pd.read_csv(args.channel_stats)
    standardization = pd.read_csv(args.standardization).set_index("channel_name").loc[channels["channel_name"]]
    x = np.load(args.image_1um).astype(np.float32)
    x -= standardization["mean"].to_numpy(np.float32)[:, None, None]
    x /= (standardization["std"].to_numpy(np.float32) + 1e-8)[:, None, None]
    print(f"standardized stack {x.shape}: {', '.join(channels['channel_name'])}", flush=True)

    # channel identities: row of each channel's protein in the panel's embedding table
    embeddings = torch.load(args.embeddings)
    channel_ids = torch.tensor([embeddings["uniprot_ids"].index(u) for u in channels["uniprot_id"]], dtype=torch.long)

    # instance mask and class logits, tile by tile (VirTues stitches tissues larger than 1500 px itself)
    head = load_model(embeddings["embeddings"], args.model_dir, args.model)
    with torch.no_grad():
        instances, logits = head.segment_tissue(torch.from_numpy(x).cuda(), channel_ids.cuda(), tile_size=args.tile_size,
                                                overlap=args.overlap, batch_size=args.batch_size)
    instances, logits = instances.cpu(), logits.cpu()
    del x

    # class map: most likely class of every pixel, all head classes (None = background, Unknown)
    class_map = torch.argmax(logits, dim=0).numpy().astype(np.uint8)
    del logits
    tifffile.imwrite(args.out_class_map, class_map, tile=(1024, 1024), compression="zlib")
    classes = list(enumerate(HEAD_CLASSES))
    pd.DataFrame(classes, columns=["class_index", "cell_type"]).to_csv(args.out_classes, index=False)
    counts = np.bincount(class_map.ravel(), minlength=len(HEAD_CLASSES))
    print("pixels per class: " + ", ".join(f"{c} {counts[i] / class_map.size * 100:.1f}%" for i, c in classes), flush=True)

    # whole cells at 1 µm, small fragments dropped
    cells = remove_small_cells(instances.numpy(), min_cell_size=MIN_CELL_SIZE)

    # whole-cell mask on the image's pixel grid: each pixel takes the 1 µm pixel it falls in
    with MultiplexImage(args.image) as img:
        height, width = img.shape[1:]
    rows = np.minimum((np.arange(height) * cells.shape[0] / height).astype(np.int64), cells.shape[0] - 1)
    cols = np.minimum((np.arange(width) * cells.shape[1] / width).astype(np.int64), cells.shape[1] - 1)
    write_label_mask(args.out_cells, cells[rows][:, cols])


if __name__ == "__main__":
    main()
