"""VirTues helpers shared by its segmentation (virtues_segmentation.py) and cell embedding (virtues_embeddings.py) scripts.

The input stack: the 1 µm stack of virtues_prepare_image.py standardized with the panel's mean / std (spora
standardization, second half), then smoothed by the 3 x 3 Gaussian blur of the VirTues paper's preprocessing.
The backbone: VirTues (Hugging Face bunnelab/virtues) with our marker embedding table, on the GPU.
"""

import numpy as np
import pandas as pd
import torch
from scipy import ndimage as ndi


def standardized_stack(image_1um, channel_stats, standardization, blur_sigma=1.0, pad=0):
    """(channels, Y, X) float32 stack VirTues reads, and the channel table (channel order of the stack).

    The stack is first zero-padded by pad px on every side (VirTues' cell token demo pads the raw image before
    standardization), then (log1p of the clipped values - panel mean) / panel std, then a Gaussian blur of every
    channel with a 3 x 3 kernel of standard deviation blur_sigma and reflected borders (as torchvision's
    gaussian_blur); 0 = no blur.
    """
    # standardization: panel mean / std of every channel, in the stack's channel order
    channels = pd.read_csv(channel_stats)
    stats = pd.read_csv(standardization).set_index("channel_name").loc[channels["channel_name"]]

    # float32 stack with pad px of zeros around it (the float16 file is read lazily, without a full copy)
    raw = np.load(image_1um, mmap_mode="r")
    x = np.zeros((raw.shape[0], raw.shape[1] + 2 * pad, raw.shape[2] + 2 * pad), dtype=np.float32)
    x[:, pad:pad + raw.shape[1], pad:pad + raw.shape[2]] = raw
    del raw

    x -= stats["mean"].to_numpy(np.float32)[:, None, None]
    x /= (stats["std"].to_numpy(np.float32) + 1e-8)[:, None, None]

    # 3 x 3 Gaussian blur, separable: the same 1D kernel along rows then columns
    if blur_sigma > 0:
        kernel = np.exp(-0.5 * (np.array([-1.0, 0.0, 1.0]) / blur_sigma) ** 2)
        kernel = (kernel / kernel.sum()).astype(np.float32)
        for c in range(x.shape[0]):
            x[c] = ndi.correlate1d(ndi.correlate1d(x[c], kernel, axis=0, mode="mirror"), kernel, axis=1, mode="mirror")
    print(f"standardized stack {x.shape} (pad {pad}, blur sigma {blur_sigma}): {', '.join(channels['channel_name'])}",
          flush=True)
    return x, channels


def channel_ids(embeddings, channels):
    """Row of each channel's protein in the panel's marker embedding table (VirTues' channel identities)."""
    return torch.tensor([embeddings["uniprot_ids"].index(u) for u in channels["uniprot_id"]], dtype=torch.long)


def load_backbone(prior_embeddings, model_dir, model, device="cuda"):
    """VirTues backbone on the GPU (device), in eval mode; its weights are downloaded to model_dir on first use."""
    from huggingface_hub import hf_hub_download
    from safetensors.torch import load_file
    from virtues.modules.multiplex_virtues import MultiplexVirtues

    # the marker embedding table is ours (a non-persistent buffer, not part of the weights)
    backbone = MultiplexVirtues(prior_bias_embeddings=prior_embeddings)
    backbone.load_state_dict(load_file(hf_hub_download("bunnelab/virtues", f"{model}/model.safetensors", local_dir=model_dir)))
    return backbone.to(device).eval()
