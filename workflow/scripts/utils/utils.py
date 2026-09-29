"""Helpers to bring proteomics cells into Xenium coordinates and match them to Xenium cells.

(The former COMET-specific helpers, spatial_to_adata for Horizon/SPACEc tables, transf_centroids
and map_metadata, are in the git history; they are superseded by these table-based versions.)
"""

import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors


def apply_affine(xy, matrix):
    """Apply a 3x3 affine matrix to (N, 2) x, y coordinates."""
    xy = np.asarray(xy, dtype=float)
    homogeneous = np.column_stack([xy, np.ones(len(xy))])
    return (homogeneous @ np.asarray(matrix, dtype=float).T)[:, :2]


def mutual_nearest_neighbours(xy_a, xy_b, max_dist):
    """Pairs (i in a, j in b) that are each other's nearest neighbour and at most max_dist apart.

    Returns a DataFrame with columns index_a, index_b and dist.
    """
    # nearest neighbour of every a in b, and of every b in a
    dist_ab, nn_ab = NearestNeighbors(n_neighbors=1).fit(xy_b).kneighbors(xy_a)
    _, nn_ba = NearestNeighbors(n_neighbors=1).fit(xy_a).kneighbors(xy_b)
    nn_ab, nn_ba, dist_ab = nn_ab[:, 0], nn_ba[:, 0], dist_ab[:, 0]

    # keep pairs that point at each other and are close enough
    index_a = np.arange(len(xy_a))
    mutual = (nn_ba[nn_ab] == index_a) & (dist_ab <= max_dist)
    return pd.DataFrame({"index_a": index_a[mutual], "index_b": nn_ab[mutual], "dist": dist_ab[mutual]})
