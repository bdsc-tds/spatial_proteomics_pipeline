# Protein preprocessing (arcsinh, scaling, K-NN graph, clustering, UMAP), copied from
# analysis_repo workflow/scripts/general/preprocessing.py: preprocess_protein and its helpers.

import importlib

import numpy as np
import pandas as pd
import scanpy as sc


def _run_som_clustering(feature_mat, xdim=10, ydim=10, num_passes=10, max_k=20, seed=42,
                         sigma=1.0, learning_rate=0.5):
    
    """
    Cluster a feature matrix with a self-organizing map, then collapse SOM nodes into a
    smaller number of meta-clusters, following the Pixie/Nimbus SOM + meta-clustering strategy
    (implemented here with MiniSom + sklearn instead of ark-analysis).

    SOM training uses MiniSom's online update rule (one random sample updates the weights per
    iteration) with a fixed random_seed, matching ark's own SOM implementation
    (`pyFlowSOM`/`flowsom.c`'s `C_SOM`), which is itself a seeded online SOM, not a batch SOM.

    :param feature_mat: array-like of shape (n_obs, n_features) to cluster
    :param xdim: number of SOM nodes along the x-axis
    :param ydim: number of SOM nodes along the y-axis
    :param num_passes: number of training passes over the full dataset
    :param max_k: number of meta-clusters to collapse the SOM nodes into
    :param seed: random seed for SOM initialization/training and meta-clustering
    :param sigma: neighbourhood radius for SOM training
    :param learning_rate: initial learning rate for SOM training
    """
    if importlib.util.find_spec("minisom") is None:
        raise ImportError(
            'clustering_method="SOM" requires the `MiniSom` package. '
            "Add it to the pixi environment to use SOM clustering."
        )
    from minisom import MiniSom
    from sklearn.cluster import AgglomerativeClustering

    X = np.asarray(feature_mat, dtype=np.float64)
    n_obs, n_features = X.shape

    som = MiniSom(xdim, ydim, n_features, sigma=sigma, learning_rate=learning_rate,
                  random_seed=seed)
    som.pca_weights_init(X)
    # random_order=True: pick samples in (seeded) random order rather than cyclically,
    # matching ark's own online SOM, which draws a random sample per iteration
    som.train(X, n_obs * num_passes, random_order=True, verbose=False)

    winners = np.array([som.winner(x) for x in X])
    som_labels = winners[:, 0] * ydim + winners[:, 1]

    # average feature profile per occupied SOM node, used to meta-cluster the nodes
    node_ids = np.unique(som_labels)
    node_profiles = np.array([X[som_labels == node].mean(axis=0) for node in node_ids])

    n_meta = min(max_k, len(node_ids))
    meta_cluster = AgglomerativeClustering(n_clusters=n_meta)
    node_meta_labels = meta_cluster.fit_predict(node_profiles)

    node_to_meta = dict(zip(node_ids, node_meta_labels))
    return np.array([node_to_meta[node] for node in som_labels])



# function to preprocess protein data modality
def _parse_pcs(n_pcs, ncomps_pca):
    """
    Convert n_pcs (None, int, list/range of 1-based PC numbers, or a "start:end" string with
    inclusive ends) into a sorted list of 1-based PC numbers, checked against ncomps_pca.
    """
    # default and int: first n PCs
    if n_pcs is None:
        pcs = list(range(1, ncomps_pca + 1))
    elif isinstance(n_pcs, (int, np.integer)):
        pcs = list(range(1, int(n_pcs) + 1))
    # "start:end" string, inclusive on both ends
    elif isinstance(n_pcs, str):
        start, end = (int(s) for s in n_pcs.split(":"))
        pcs = list(range(start, end + 1))
    # list / range / array of PC numbers
    else:
        pcs = [int(p) for p in n_pcs]

    # unique, sorted, and within the computed PCs
    pcs = sorted(set(pcs))
    if not pcs:
        raise ValueError(f"n_pcs ({n_pcs}) selects no PCs")
    if pcs[0] < 1 or pcs[-1] > ncomps_pca:
        raise ValueError(f"n_pcs ({n_pcs}) must select PCs between 1 and ncomps_pca ({ncomps_pca})")

    return pcs


def preprocess_protein(adata,
                        pca=False,
                        scale=True,
                        upper_clip=None,
                        clustering_method="leiden",
                        leiden_resolution=0.8,
                        leiden_key=None,
                        nn=15,
                        raw_integers=False,
                        trans="arcsinh",
                        arcsinh_cofactor=5,
                        ncomps_pca=None,
                        n_pcs=None,
                        umap_key = "X_umap",
                        nn_key = "neighbors",
                        exclude=["DAPI", 'TRITC-autof', 'Cy5-autof'],
                        include=None,
                        som_xdim=10,
                        som_ydim=10,
                        som_num_passes=10,
                        som_max_k=20,
                        som_seed=42,
):


    """

    Function to run standard processing of protein data

    Parameters
    ------------
    adata: class, AnnData
        anndata object containing the protein data. The feature raw feature values are expecte be to in .layers["raw"]
    pca: bool
        whether to perform PCA  
    scale: bool
        whether to scale the data
    clustering_method: str
        clustering method to use, can be "leiden" or "SOM". "SOM" trains a self-organizing map
        (via MiniSom) on the same feature/PCA matrix used for the KNN graph, then collapses the
        SOM nodes into som_max_k meta-clusters via hierarchical clustering on the per-node average
        profiles (the SOM + meta-clustering strategy used in the Pixie/Nimbus pipelines).
    leiden_resolution: float
        resolution parameter for Leiden clustering. Ignored if clustering_method is "SOM"
    leiden_key: str
        key used to store the cluster labels. If not specified, f'leiden_{leiden_resolution}' is
        used for clustering_method="leiden" and f'som_meta_{som_max_k}' for clustering_method="SOM"
    nn: int
        Number of nearest neighbours to compute K-NN graph. Deafult is 15.
    raw_integers: bool  
        boolean indicating if raw features should be rounded to integers (they do this in the processing of Xenium RNA+Protein)
    trans: str
        transformation to apply to the data. Can be "arcsinh" or "log1p"
    arcsinh_cofactor: int
        cofactor to use for arcsinh transformation. Default is 5
    ncomps_pca: int
        number of PCA components to compute. Defaults to n_proteins-1.
    n_pcs: int, list, range or str
        PCA components used to compute the K-NN graph (and the SOM clustering). An int n uses the
        first n PCs; a list/range of 1-based PC numbers (e.g. [2, 3, 4, 5] or range(1, 11)) or a
        "start:end" string with inclusive 1-based ends (e.g. "1:10" or "2:5") uses exactly those
        PCs. All PCs must be <= ncomps_pca. Defaults to all ncomps_pca PCs. Only used if pca is True.
    umap_key: str
        key to use for storing UMAP coordinates in .obsm. Default is "X_umap"
    nn_key: str
        key to use for storing neighbors in .uns. Default is "neighbors"
    exclude: list
        list of markers to exclude from the pre-processing. Defaults to ["DAPI", 'TRITC-autof', 'Cy5-autof']
    include: list
        list of markers on which to run the PCA, K-NN graph and clustering (markers also in exclude
        are still dropped). All markers are still transformed and scaled. If None (default), all
        markers not in exclude are used for the K-NN graph, and all markers for the PCA.
    som_xdim: int
        number of SOM nodes along the x-axis. Only used if clustering_method is "SOM"
    som_ydim: int
        number of SOM nodes along the y-axis. Only used if clustering_method is "SOM"
    som_num_passes: int
        number of SOM training passes over the full dataset. Only used if clustering_method is "SOM"
    som_max_k: int
        number of meta-clusters to collapse the som_xdim*som_ydim SOM nodes into. Only used if
        clustering_method is "SOM"
    som_seed: int
        random seed for SOM training and meta-clustering. Only used if clustering_method is "SOM"

    """

    # make sure raw feature matrix is in .X: create the "raw" layer on first call, otherwise
    # reset .X from it so repeated calls (e.g. leiden then SOM) always start from raw data
    # instead of re-transforming whatever the previous call left in .X
    if "raw" not in adata.layers:
        adata.layers["raw"] = adata.X.copy()
    else:
        adata.X = adata.layers["raw"].copy()
    raw = ""

    # markers used for PCA, K-NN graph and clustering: the include list (if given) minus exclude
    if include is not None:
        missing = [f for f in include if f not in adata.var_names]
        if missing:
            raise ValueError(f"Markers in include not found in adata.var_names: {missing}")
    feature_mask = np.array([(include is None or f in include) and f not in exclude for f in adata.var_names])
    feature_idx = np.where(feature_mask)[0].tolist()

    if ncomps_pca is None:
        ncomps_pca = (adata.n_vars if include is None else len(feature_idx)) - 1

    # PCs used for the K-NN graph as 1-based PC numbers, defaults to all computed PCs
    pcs = _parse_pcs(n_pcs, ncomps_pca)

    # round to integers if specified
    if raw_integers:

        print("Rounding raw feature values to integers...")

        adata.layers["raw_int"] = np.round(adata.layers["raw"]).astype(np.int32)

        # move it to .X
        adata.X = adata.layers["raw_int"].copy()
        raw="_int"

    # log1p transform the data and add as a new layer
    if trans == "log1p":

        print("Applying log1p transformation...")

        sc.pp.log1p(adata)
        adata.layers[f'log1p{raw}'] = adata.X.copy()

    else:

        print(f'Applying arcsinh transformation with cofactor {arcsinh_cofactor}')

        adata.layers[f'arcsinh{raw}'] = np.arcsinh(adata.X / arcsinh_cofactor)
        adata.X = adata.layers[f'arcsinh{raw}'].copy()

    if scale:

        print("Scaling data")

        # remove previous mean and variance computed (otherwise it doesn't overwrite them...)
        for col in ["mean", "std"]:
            if col in adata.var.columns:
                adata.var.drop(columns=[col], inplace=True)

        if upper_clip is not None:

            sc.pp.scale(adata, max_value=upper_clip)

        else:

            sc.pp.scale(adata)

    # scanpy special-cases key_added="neighbors" to mean "use the default, unprefixed
    # .obsp['connectivities']/['distances'] storage", but sc.tl.leiden/sc.tl.umap's
    # neighbors_key="neighbors" always looks there too regardless of what was passed to
    # key_added. Passing key_added="neighbors" explicitly (rather than omitting it) makes
    # sc.pp.neighbors store under the prefixed .obsp['neighbors_connectivities'] instead,
    # which leiden/umap then can't find. Passing None instead of the literal "neighbors"
    # keeps storage consistent with what leiden/umap expect.
    neighbors_key_added = None if nn_key == "neighbors" else nn_key

    print(f"Computing KNN graph on PCs {pcs}" if pca else "Computing KNN graph")
    if pca:

        key= f'pca_{trans}{raw}'

        # restrict the PCA to the selected markers only when an include list is given
        sc.pp.pca(adata, n_comps=ncomps_pca, zero_center=False,
                  mask_var=feature_mask if include is not None else None)

        # first n PCs: use n_pcs directly; any other selection: build the K-NN graph on a PC subset matrix
        if pcs == list(range(1, len(pcs) + 1)):
            sc.pp.neighbors(adata, n_neighbors=nn, n_pcs=len(pcs), metric="cosine", key_added=neighbors_key_added)
        else:
            adata.obsm["X_pca_knn"] = adata.obsm["X_pca"][:, [p - 1 for p in pcs]]
            sc.pp.neighbors(adata, n_neighbors=nn, use_rep="X_pca_knn", metric="cosine", key_added=neighbors_key_added)

        # matrix used for SOM clustering, kept in sync with what the KNN graph was built on
        cluster_feature_mat = adata.obsm["X_pca"][:, [p - 1 for p in pcs]]

    else:

        key = f'scaled_{trans}{raw}'

        # add an aritifical matrix to obsm to use for dim.reduction and clustering
        adata.obsm[key] = adata.X[:,feature_idx].copy()

        # compute K-NN graph
        sc.pp.neighbors(adata, n_neighbors=nn, use_rep=key, metric="cosine", key_added=neighbors_key_added)

        cluster_feature_mat = adata.obsm[key]

    if clustering_method not in ("leiden", "SOM"):
        raise ValueError(f'Unknown clustering_method "{clustering_method}", must be "leiden" or "SOM"')

    if leiden_key is None:
        leiden_key = f'leiden_{leiden_resolution}' if clustering_method == "leiden" \
            else f'som_meta_{som_max_k}'

    print("Clustering KNN graph")
    if clustering_method == "leiden":
        sc.tl.leiden(adata,resolution=leiden_resolution, flavor="igraph", n_iterations=2,
                     key_added=leiden_key, neighbors_key=nn_key)
    else:
        som_labels = _run_som_clustering(
            cluster_feature_mat,
            xdim=som_xdim,
            ydim=som_ydim,
            num_passes=som_num_passes,
            max_k=som_max_k,
            seed=som_seed,
        )
        adata.obs[leiden_key] = pd.Categorical(som_labels.astype(str))
    print("Computing UMAP")
    sc.tl.umap(adata, neighbors_key=nn_key, key_added=umap_key, min_dist=0.3)

    return adata
