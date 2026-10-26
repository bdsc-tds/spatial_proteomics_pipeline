"""Write the cell, nuclei and cytoplasm label masks of a segmentation method, all sharing cell_id.

Shared by the segmentation scripts (mesmer.py, cellpose.py) so that every method hands the same
three masks to quantify_cells and mask_to_geojson. Masks are written as tiled, deflate-compressed
TIFFs (int32), which downstream steps read lazily tile by tile and QuPath / napari open directly.
"""

import numpy as np
import tifffile


def save_label_masks(cells, in_nucleus, out_cells, out_nuclei, out_cytoplasm):
    """Save the cell mask and the nuclei / cytoplasm masks relabelled with the id of their cell.

    cells: whole-cell label mask (int). in_nucleus: boolean mask of nuclear pixels, from a separate
    nuclear segmentation whose labels are ignored. Nuclear pixels outside any cell are dropped.
    """
    cells = cells.astype(np.int32, copy=False)

    # derive and save the label masks one at a time (a full-resolution Xenium mask is ~9 GB as int32)
    derived = {
        out_cells: lambda: cells,
        out_nuclei: lambda: np.where(in_nucleus, cells, 0),
        out_cytoplasm: lambda: np.where(in_nucleus, 0, cells),
    }
    for path, make_mask in derived.items():
        mask = make_mask()
        # 1024 px tiles, deflate (no predictor, which compresses label masks worse)
        tifffile.imwrite(path, mask, tile=(1024, 1024), compression="zlib", bigtiff=True)
        # count labels without np.unique, which would sort a copy of the whole mask
        n_labels = np.count_nonzero(np.bincount(mask.ravel())[1:])
        print(f"written: {path} ({n_labels} labels)")
        del mask
