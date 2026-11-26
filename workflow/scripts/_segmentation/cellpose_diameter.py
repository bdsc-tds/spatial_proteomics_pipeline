"""Per-sample estimate of the Cellpose cell and nucleus diameters, with Cellpose's own size models.

Cellpose's size model (models.Cellpose with diameter 0) predicts the object diameter from the network's style
vector, segments the image at that size and returns the median diameter of the objects found. SPACEc cannot be
used for this: it runs cyto3 through CellposeModel, which has no size model (it falls back to the 30 px training
diameter), and it would estimate the nuclear diameter again on every tile. The estimate is therefore made once per
sample, on a few tissue-rich crops, and the median over the crops is passed to SPACEc as a fixed diameter.

Weights: the size models (size_cyto3.npy, size_nucleitorch_0.npy) sit next to the Cellpose models and are
downloaded into CELLPOSE_LOCAL_MODELS_PATH on first use, like the models themselves.
"""

import numpy as np

# diameter value in the config that switches the estimate on
AUTO = "auto"


def parse_diameter(value):
    """Config / command-line diameter: "auto" (estimated per sample) or a number in µm."""
    if str(value).strip().lower() == AUTO:
        return AUTO
    try:
        return float(value)
    except ValueError:
        raise ValueError(f"diameter must be a number (µm) or {AUTO!r}, got {value!r}") from None


def tissue_crops(nuclear, crop_px, n_crops, downsample=8):
    """Top-left corners (row, col) of the n_crops crop_px x crop_px blocks with the most nuclear signal.

    Blocks are scored by the fraction of pixels above an Otsu threshold of the (downsampled) nuclear plane, so
    background and empty slide areas are left out. An image smaller than one crop gives a single crop at (0, 0).
    """
    from skimage.filters import threshold_otsu

    # images smaller than one crop: the whole image
    ny, nx = nuclear.shape[0] // crop_px, nuclear.shape[1] // crop_px
    if ny == 0 or nx == 0:
        return [(0, 0)]

    # foreground on a downsampled copy, enough to rank blocks
    small = nuclear[::downsample, ::downsample].astype(np.float32)
    foreground = small > threshold_otsu(small)

    # fraction of foreground pixels per block
    block = crop_px // downsample
    grid = foreground[: ny * block, : nx * block].reshape(ny, block, nx, block).mean(axis=(1, 3))

    # the blocks with the most foreground (blocks without any are never picked)
    order = np.argsort(grid, axis=None)[::-1][:n_crops]
    picks = [np.unravel_index(i, grid.shape) for i in order if grid.flat[i] > 0]
    return [(int(r) * crop_px, int(c) * crop_px) for r, c in picks] or [(0, 0)]


def estimate_diameter_px(crops, model_type, channels, use_gpu):
    """Median Cellpose size-model diameter (px) over the crops, and the per-crop diameters.

    crops: list of 2D planes (channels [0, 0]) or (Y, X, 3) RGB images (channels e.g. [1, 3]).
    A crop where the size model finds no objects gets the model's mean training diameter back from Cellpose; these
    crops are left out of the median (all crops are used if every one of them failed).
    """
    from cellpose import models

    # Cellpose (not CellposeModel) carries the size model; its SizeModel is used directly so no full masks are kept
    model = models.Cellpose(gpu=use_gpu, model_type=model_type)
    diams, _ = model.sz.eval(crops, channels=channels)
    diams = np.atleast_1d(np.asarray(diams, dtype=float))

    # drop crops where the size model fell back to its training diameter
    found = diams[~np.isclose(diams, model.sz.diam_mean)]
    median = float(np.median(found if found.size else diams))
    return median, diams.tolist()


def auto_diameters(nuclear, membrane, pixel_size, cell_model, nuclear_model, estimate_cells, estimate_nuclei,
                   cell_range_um, nuclear_range_um, use_gpu, crop_um=500, n_crops=8):
    """Estimated cell and/or nucleus diameters (µm) of one sample, with the details of the estimate.

    nuclear, membrane: full-resolution segmentation input planes (Y, X).
    estimate_cells / estimate_nuclei: which of the two diameters to estimate (the other is left as None).
    cell_range_um / nuclear_range_um: (min, max) µm; estimates outside are clamped to it, with a warning.
    Returns {"cells": {...} or None, "nuclei": {...} or None, "crops": [(row, col), ...], "crop_px": int}; each
    estimate holds diameter_um (after clamping), estimated_um, per_crop_um and clamped.
    """
    # tissue-rich crops, the same for both estimates
    crop_px = max(int(round(crop_um / pixel_size)), 224)
    corners = tissue_crops(nuclear, crop_px, n_crops)
    window = lambda plane, r, c: plane[r : r + crop_px, c : c + crop_px]
    print(f"diameter estimate on {len(corners)} crops of {crop_px} px ({crop_px * pixel_size:.0f} µm)", flush=True)

    def summarise(name, median_px, per_crop_px, value_range):
        # µm, clamped to the allowed range
        estimated = median_px * pixel_size
        low, high = value_range
        diameter = float(min(max(estimated, low), high))
        if diameter != estimated:
            print(f"WARNING: estimated {name} diameter {estimated:.1f} µm outside {low}-{high} µm, clamped to "
                  f"{diameter:.1f} µm", flush=True)
        per_crop = [round(d * pixel_size, 2) for d in per_crop_px]
        print(f"{name}: per-crop {per_crop} µm -> median {estimated:.2f} µm, using {diameter:.2f} µm", flush=True)
        return {"diameter_um": diameter, "estimated_um": estimated, "per_crop_um": per_crop,
                "clamped": diameter != estimated}

    result = {"cells": None, "nuclei": None, "crops": corners, "crop_px": crop_px}

    # nuclei: nuclear model on the nuclear signal alone, as in segmentation
    if estimate_nuclei:
        crops = [window(nuclear, r, c) for r, c in corners]
        median_px, per_crop_px = estimate_diameter_px(crops, nuclear_model, [0, 0], use_gpu)
        result["nuclei"] = summarise("nucleus", median_px, per_crop_px, nuclear_range_um)

    # cells: whole-cell model on membrane (red) + nuclear (blue), the RGB layout SPACEc builds
    if estimate_cells:
        crops = [np.stack([window(membrane, r, c), np.zeros_like(window(nuclear, r, c)), window(nuclear, r, c)],
                          axis=-1) for r, c in corners]
        median_px, per_crop_px = estimate_diameter_px(crops, cell_model, [1, 3], use_gpu)
        result["cells"] = summarise("cell", median_px, per_crop_px, cell_range_um)

    return result
