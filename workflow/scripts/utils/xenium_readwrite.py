"""Xenium readers shared with the analysis_repo repository.

read_xenium_sample() and the helpers it calls are copied verbatim from
analysis_repo/workflow/scripts/general/readwrite.py (commit 444fdc4), so that Xenium data is read
the same way in both projects. Do not edit them here: change them there and copy them over again.

Kept separate from readwrite.py because it needs spatialdata, which the spacec environment (Mesmer)
does not have.
"""

import json
import re
import warnings
from pathlib import Path
from types import MappingProxyType

import anndata as ad
import dask.dataframe as dd
import geopandas as gpd
import muon as mu
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
import spatialdata
import spatialdata_io
from sklearn.neighbors import NearestNeighbors
from spatialdata.models import (
    ShapesModel,
    PointsModel,
)


# function to add cell type labels to anndata object. Input should be an anndata object from read_xenium_sample()
def add_ct_labels(adata, 
                  xp_dir, 
                  full_sample_id,
                  seg_id, 
                  norm_id, 
                  ct_ref, 
                  ct_levels,
                  count_correction=False,
                  correction_mode = "split_fully_purified",
                  rctd_type="rctd_class_aware",
                  ct_results_file="results_df",
                  full_post_processed_info=False,

):
    """

    Add cell type information from Xenium pipeline to the .obs dataframe of anndata object.

    Parameters 
    ----------
    adata : anndata.AnnData
        Anndata object from the xenium pipeline output. 
    xp_dir : str
        Path to the processed folder of the xenium pipeline
    full_sample_id : str
        Full sample ID specified in the Xenium pipeline. Format should be: {condition}_{panel}_{donor_id}_{sample_id}
    seg_id : str
        Segmentation ID. Options are: "proseg", "10x_mm_5um", "10x_5um"
    norm_id : str
        Normalization ID used for the cell type annotation results to be added
    ct_ref : str
        Name give the to the scRNAseq reference used for annotation in the pipeline.
    ct_levels : str or list
        String or list indicating which level of granularity to be added
    count_correction : bool, optional
        Whether cell type labels before or after count correction should be added. Default is False.
    correction_mode : str, optional
        Correction mode used for SPLIT. Default="split_fully_purified"
    rctd_type : str, optional
        Mode used when running RCTD. Default="rctd_class_aware"
    ct_results_file : str, optional
        Which RCTD output table to read the labels from. Options are "results_df" (raw RCTD output) or
        "post_processed_results_df" (SPLIT post-processed output). Default is "results_df".
    full_post_processed_info : bool, optional
        Whether to also add the extra post-processing columns (second_type_class, weight_first_type,
        weight_second_type, singlet_score) for each level. Only available for
        ct_results_file="post_processed_results_df". Default is False.

    Returns
    --------

    adata : anndata.AnnData
        Anndata object with added cell type annotation in .obs dataframe. The cell type columns are named {ct_ref}_{ct_level} (e.g. "ref1_Level1") and contain the "first_type" column from the RCTD output for the specified level of granularity.
        The "spot_class" column is always added as {ct_level}_spot_class, and with full_post_processed_info the extra post-processing columns are added as {ct_level}_{column}. In both cases the level name is trimmed to the level number (e.g. "Level4_partner" -> "Level4_spot_class").

    """

    # validate the requested results table and the extra-columns option
    if ct_results_file not in ["results_df", "post_processed_results_df"]:
        raise ValueError('ct_results_file must be either "results_df" or "post_processed_results_df"')

    if full_post_processed_info and ct_results_file != "post_processed_results_df":
        raise ValueError('full_post_processed_info=True requires ct_results_file="post_processed_results_df"')

    # extra post-processing columns to carry over when full info is requested
    extra_cols = ["second_type", "weight_first_type", "weight_second_type", "singlet_score"] if full_post_processed_info else []

    # add cell type annotation to adata
    # note: adata.obs["cell_id"] (and obs_names) is already "proseg-"-prefixed for proseg
    # samples by this point (read_xenium_sample sets it right after loading), matching the
    # cell_id convention used by the RCTD/count-correction outputs read below
    meta = adata.obs.copy()

    # determine if to take cell types before or after count correction
    if count_correction:
        corr_id = "post_count_correction_cell_type_annotation"
        path_end = f'single_cell/{correction_mode}/{norm_id}/reference_based/{ct_ref}/{rctd_type}/'

    else:
        corr_id = "cell_type_annotation"


    for lvl in ct_levels:

        # the path is so different when there is count correction before...
        if count_correction:

            ct_df = pd.read_parquet(f'{xp_dir}/{corr_id}/{seg_id}/{full_sample_id}/{norm_id}/reference_based/{ct_ref}/{rctd_type}/{lvl}/{path_end}/{lvl}/single_cell/output/{ct_results_file}.parquet')

        else:

            # load cell type annotation results for the segmentation
            ct_df = pd.read_parquet(f'{xp_dir}/{corr_id}/{seg_id}/{full_sample_id}/{norm_id}/reference_based/{ct_ref}/{rctd_type}/{lvl}/single_cell/output/{ct_results_file}.parquet')

        ct_df[f'{lvl}_{ct_ref}'] = ct_df["first_type"]

        # trim the level name down to the level number (e.g. "Level4_partner" -> "Level4")
        lvl_short = re.sub(r"(\d+).*$", r"\1", lvl)

        # always carry the spot_class of the level, plus the extra post-processing columns if requested
        keep_cols = [f'{lvl}_{ct_ref}']
        for col in ["spot_class"] + extra_cols:

            ct_df[f'{lvl_short}_{col}'] = ct_df[col]
            keep_cols.append(f'{lvl_short}_{col}')

        # left join cell type annotations
        meta = meta.merge(ct_df[["cell_id"] + keep_cols], on="cell_id", how="left").set_index(adata.obs_names)

    # over-write metadata df
    adata.obs = meta

    return adata


# function to add SPLIT matrix to adata object
def add_split_matrix(adata, 
                     xp_dir,
                    seg_id, 
                    full_sample_id, 
                    norm_id,
                    ct_ref=None, 
                    ct_level=None, 
                    correction_mode = "split_fully_purified",
                    rctd_type="rctd_class_aware"): 

    """

    Add count corrected SPLIT data matrix to anndata as a layer

    Parameters 
    ----------
    adata : anndata.AnnData
        Anndata object from the xenium pipeline output. 
    xp_dir : str
        Path to the processed folder of the xenium pipeline
    full_sample_id : str
        Full sample ID specified in the Xenium pipeline. Format should be: {condition}_{panel}_{donor_id}_{sample_id}
    seg_id : str
        Segmentation ID. Options are: "proseg", "10x_mm_5um", "10x_5um"
    norm_id : str
        Normalization ID used for the cell type annotation results to be added
    ct_ref : str
        Name give the to the scRNAseq reference used for annotation in the pipeline.
    ct_level : str
        Cell type level annotation used to run SPLIT
    correction_mode : str, optional
        Correction mode used for SPLIT. Default="split_fully_purified"
    rctd_type : str, optional
        Mode used when running RCTD. Default="rctd_class_aware"

    Returns
    --------

    adata : anndata.AnnData
        Anndata object with SPLIT correct count matrix stored in .X and in layers["raw_post_split"]. The pre-split count matrix is stored in layers["raw_pre_split"].

    """

    # read split corrected counts
    split_counts = sc.read_10x_h5(f'{xp_dir}/count_correction/{seg_id}/{full_sample_id}/{norm_id}/reference_based/{ct_ref}/{rctd_type}/{ct_level}/single_cell/{correction_mode}/corrected_counts.h5')

    # make sure obs_names have the correct format
    # adata.obs_names = adata.obs["cell_id"]

    # matched cell IDs
    matched = adata.obs_names.intersection(split_counts.obs_names)

    # subset both objects to cells and genes contained in the split matrix
    adata = adata[matched, split_counts.var_names]
    split_counts = split_counts[matched, :].copy()

    # add layer with pre-split count corrected matrix
    adata.layers["raw_pre_split"] = adata.X.copy()

    # add split corrected counts as new layer
    adata.layers["raw_post_split"] = split_counts.X.copy()

    # set corrected count matrix as default assay in X
    adata.X = adata.layers["raw_post_split"].copy()

    return adata


def _get_nearest_neighbors(ref_obs, query_obs, ref_col, query_col, x_centroid="x_centroid", y_centroid="y_centroid"):
    """
    For each cell in query_obs, find its nearest neighbour (by centroid) in ref_obs.

    ref_obs, query_obs : pd.DataFrame
        obs dataframes indexed by (unique) cell obs_names, with x_centroid/y_centroid columns.

    Returns a dataframe with columns [ref_col, query_col, "dist"], keeping only the
    closest query cell for each matched reference cell.
    """
    dist, idx = (
        NearestNeighbors(n_neighbors=1)
        .fit(ref_obs[[x_centroid, y_centroid]])
        .kneighbors(query_obs[[x_centroid, y_centroid]])
    )
    ids = [i[0] for i in idx]
    dists = [d[0] for d in dist]

    matches = pd.DataFrame({
        ref_col: ref_obs.index[ids].values,
        query_col: query_obs.index.values,
        "dist": dists,
    })

    # keep only the closest query cell for each matched reference cell (ties broken
    # by distance). Avoid groupby().apply() here: pandas >=2.2 changed its default
    # `include_groups` behaviour and drops the grouping column from what's passed
    # to the applied function, which silently breaks this kind of "keep top-1 per
    # group" pattern.
    return (
        matches.sort_values("dist")
        .drop_duplicates(subset=ref_col, keep="first")
        .reset_index(drop=True)
    )


def match_cells_spatial_mnn(proseg_obs, xr_obs, x_centroid="x_centroid", y_centroid="y_centroid", max_dist=10):
    """
    Match proseg and xenium ranger cells by spatial mutual nearest neighbours (MNN).

    A proseg cell and a xenium ranger cell are considered matched if each is the
    other's nearest neighbour by centroid distance (adapted from the xenium_analysis_pipeline's spatial
    MNN cell-mapping script).

    Parameters
    ----------
    proseg_obs : pd.DataFrame
        adata.obs of the proseg sample, indexed by (unique) proseg obs_names, with
        x_centroid/y_centroid columns.
    xr_obs : pd.DataFrame
        adata.obs of the xenium ranger sample, indexed by (unique) xenium ranger
        obs_names, with x_centroid/y_centroid columns.
    max_dist : float, optional
        Maximum centroid distance allowed for a mutually-matched pair to be kept.
        Matched pairs further apart than this are dropped. Default 10.

    Returns
    -------
    pd.DataFrame with columns "proseg_cell", "xr_cell" and "dist", one row per
    mutually matched cell pair within max_dist (unmatched or too-distant cells on
    either side are dropped).
    """
    proseg2xr = _get_nearest_neighbors(xr_obs, proseg_obs, "xr_cell", "proseg_cell", x_centroid, y_centroid)
    xr2proseg = _get_nearest_neighbors(proseg_obs, xr_obs, "proseg_cell", "xr_cell", x_centroid, y_centroid)

    matches = pd.merge(
        proseg2xr[["proseg_cell", "xr_cell", "dist"]],
        xr2proseg[["proseg_cell", "xr_cell"]],
        how="inner",
        on=["proseg_cell", "xr_cell"],
    )

    return matches.loc[matches["dist"] <= max_dist].reset_index(drop=True)


def read_xenium_specs(xenium_specs_file):
    xenium_specs_file = Path(xenium_specs_file)
    with open(xenium_specs_file) as f:
        specs = json.load(f)
    return specs


def xenium_proseg(
    path,
    cells_boundaries=True,
    cells_boundaries_layers=True,
    nucleus_boundaries=False,
    cells_as_circles=False,
    cells_labels=False,
    nucleus_labels=False,
    transcripts=True,
    morphology_mip=False,
    morphology_focus=False,
    aligned_images=False,
    cells_table=True,
    n_jobs=1,
    imread_kwargs=MappingProxyType({}),
    image_models_kwargs=MappingProxyType({}),
    labels_models_kwargs=MappingProxyType({}),
    cells_metadata=True,
    xeniumranger_dir=None,
    xenium_specs=True,
    pandas_engine="pyarrow",
    verbose=False,
):
    """
    Reads a Xenium segmentation run and returns a `SpatialData` object.

    Parameters
    ----------
    path : str
        The directory path of the segmentation run.
    cells_boundaries : bool or str, optional
        The file path of the cell boundaries GeoJSON file.
    cells_boundaries_layers : bool or str, optional
        The file path of the cell boundaries layers GeoJSON file.
    nucleus_boundaries : bool or str, optional
        Whether to read nucleus boundaries. Not implemented.
    cells_as_circles : bool or str, optional
        Whether to read cells as circles.
    cells_labels : bool or str, optional
        Whether to read cells labels. Not implemented.
    nucleus_labels : bool or str, optional
        Whether to read nucleus labels. Not implemented.
    transcripts : bool or str, optional
        The file path of the transcripts CSV file.
    morphology_mip : bool, optional
        Whether to read morphology MIP images.
    morphology_focus : bool, optional
        Whether to read morphology focus images.
    aligned_images : bool, optional
        Whether to read aligned images.
    cells_table : bool or str, optional
        The file path of the cells table CSV file.
    n_jobs : int, optional
        The number of jobs to use. Not implemented.
    imread_kwargs : dict, optional
        Keyword arguments to pass to `imread`.
    image_models_kwargs : dict, optional
        Keyword arguments to pass to `ImageModel`.
    labels_models_kwargs : dict, optional
        Keyword arguments to pass to `LabelsModel`.
    cells_metadata : bool or str, optional
        The file path of the cells metadata CSV file.
    xeniumranger_dir : str, optional
        The directory path of the XeniumRanger run.
    xenium_specs : bool or str, optional
        The file path of the Xenium specs file.
    pandas_engine : str, optional
        The pandas engine to use when reading CSV files.

    Returns
    -------
    sdata : SpatialData
        The `SpatialData` object.
    """
    path = Path(path)

    # unsupported options compared to spatialdata_io.xenium
    if nucleus_boundaries:
        raise ValueError("reading nucleus_boundaries not implemented for proseg")
    if cells_labels:
        raise ValueError("reading cells_labels not implemented for proseg")
    if nucleus_labels:
        raise ValueError("reading nucleus_labels not implemented for proseg")
    if n_jobs > 1:
        raise ValueError("n_jobs>1 not supported")

    # default expected file paths
    def parse_arg(arg, default):
        if arg:
            return default
        elif isinstance(arg, str):
            return Path(arg)
        else:
            return arg

    cells_metadata = parse_arg(cells_metadata, path / "cell-metadata.csv.gz")
    cells_boundaries = parse_arg(cells_boundaries, path / "cell-polygons.geojson.gz")
    cells_boundaries_layers = parse_arg(
        cells_boundaries_layers, path / "cell-polygons-layers.geojson.gz"
    )
    transcripts = parse_arg(transcripts, path / "transcript-metadata.csv.gz")
    cells_table = parse_arg(cells_table, path / "expected-counts.csv.gz")
    xeniumranger_dir = parse_arg(xeniumranger_dir, None)
    if xeniumranger_dir is not None or isinstance(xenium_specs, str):
        xenium_specs = parse_arg(xenium_specs, xeniumranger_dir / "experiment.xenium")

    ### images
    if morphology_mip or morphology_focus or aligned_images:
        if verbose:
            print("Reading images...")
        sdata_images = spatialdata_io.xenium(
            xeniumranger_dir,
            cells_table=False,
            cells_as_circles=False,
            cells_boundaries=False,
            nucleus_boundaries=False,
            cells_labels=False,
            nucleus_labels=False,
            transcripts=False,
            morphology_mip=morphology_mip,
            morphology_focus=morphology_focus,
            aligned_images=aligned_images,
            imread_kwargs=imread_kwargs,
            image_models_kwargs=image_models_kwargs,
            labels_models_kwargs=labels_models_kwargs,
        )

        images = sdata_images.images
    else:
        images = {}

    ### tables
    region = "cell_polygons"
    region_key = "region"
    instance_key = "cell_id"

    # flag columns not corresponding to genes
    if isinstance(cells_table, Path):
        if verbose:
            print("Reading cells table...")
        df_table = pd.read_csv(cells_table, engine=pandas_engine, dtype=np.float32) # otherwise the float64 matrix is huge

        control_columns = df_table.columns.str.contains(
            "|".join(["BLANK_", "Codeword", "NegControl"])
        )

        table = ad.AnnData(
            df_table.iloc[:, ~control_columns],
            uns={
                "spatialdata_attrs": {
                    "region": region,
                    "region_key": region_key,
                    "instance_key": instance_key,
                }
            },
        )

        if isinstance(cells_metadata, Path):
            if verbose:
                print("Reading cells metadata...")

            df_cells_metadata = pd.read_csv(
                cells_metadata, engine=pandas_engine
            ).rename(columns={"cell": "cell_id"})
            table.obs = pd.concat(
                (df_cells_metadata, df_table.iloc[:, control_columns]), axis=1
            )
            table.obsm["spatial"] = table.obs[["centroid_x", "centroid_y"]].values
            table.obs = table.obs
        table.obs[region_key] = pd.Categorical([region] * table.n_obs)

        # sparsify .X
        table.X = scipy.sparse.csr_matrix(table.X)
        tables = {"table": table}
    else:
        tables = {}

    ### labels
    # not implemented
    labels = {}

    ### points
    if isinstance(transcripts, Path):
        if verbose:
            print("Reading transcripts...")

        df_transcripts = dd.read_csv(transcripts, blocksize=None).rename(
            columns={"gene": "feature_name", "assignment": "cell_id"}
        )
        points = {"transcripts": PointsModel.parse(df_transcripts)}
    else:
        points = {}

    ### shapes
    shapes = {}

    # read specs
    if isinstance(xenium_specs, Path):
        if verbose:
            print("Reading specs...")
        specs = read_xenium_specs(xenium_specs)

        # get xenium pixel size
        scale = spatialdata.transformations.Scale(
            [1.0 / specs["pixel_size"], 1.0 / specs["pixel_size"]], axes=("x", "y")
        )
        transformations = {"global": scale}
    else:
        transformations = None
        if isinstance(cells_boundaries, Path) or isinstance(
            cells_boundaries_layers, Path
        ):
            warnings.warn(
                """
                Couldn't load xenium specs file with pixel size. 
                Not applying scale transformations to shapes.
                Please specify xeniumranger_dir or xenium_specs
                """
            )

    # read cells boundaries
    if isinstance(cells_boundaries, Path):
        if verbose:
            print("Reading cells boundaries...")

        df_cells_boundaries = gpd.read_file(
            "gzip://" + cells_boundaries.as_posix()
        ).rename(columns={"cell": "cell_id"})
        shapes["cells_boundaries"] = ShapesModel.parse(
            df_cells_boundaries, transformations=transformations
        )

    # read cells boundaries layers
    if isinstance(cells_boundaries_layers, Path):
        if verbose:
            print("Reading cells boundaries layers...")

        df_cells_boundaries_layers = gpd.read_file(
            "gzip://" + cells_boundaries_layers.as_posix()
        ).rename(columns={"cell": "cell_id"})
        shapes["cells_boundaries_layers"] = ShapesModel.parse(
            df_cells_boundaries_layers, transformations=transformations
        )

    # convert cells boundaries to circles
    if cells_as_circles:
        if verbose:
            print("Converting cells boundaries to circle...")

        shapes["cells_boundaries_circles"] = spatialdata.to_circles(
            shapes["cells_boundaries"]
        )

    ### sdata
    sdata = spatialdata.SpatialData(
        images=images, labels=labels, points=points, shapes=shapes, tables=tables
    )

    return sdata


def _sync_shapes_to_table(sdata, table_name, adata, is_proseg=False, instance_key="cell_id"):
    """
    Keep per-cell shape elements in sync with a (possibly already subsetted/reordered) table.

    None of the row-reducing steps upstream in read_xenium_sample (cell type merge, SPLIT
    correction, proseg<->xenium ranger MNN matching, RNA/protein cell intersection) touch the
    shapes stored in sdata -- they stay at their original, full extent. Also, table.uns["spatialdata_attrs"]
    ends up pointing at the wrong (or a nonexistent) region element: spatialdata_io.xenium() always
    stamps it as "cell_labels" regardless of whether cell_boundaries was what got loaded, and
    xenium_proseg() stamps it as "cell_polygons", which isn't the name of any shapes element it
    ever creates. This subsets each loaded shape element down to the cells retained in `adata`
    and re-points the table's region metadata at it.

    The two readers disagree on how shapes are keyed, so `is_proseg` controls which convention
    is used (pass the same check read_xenium_sample uses to decide the obs_names prefix):

    - xeniumranger (spatialdata_io.xenium): shapes are "cell_boundaries"/"cell_circles", indexed
      directly by cell_id -- filtering is a plain `.loc[]` on the shapes' index.
    - proseg (xenium_proseg): shapes are "cells_boundaries"/"cells_boundaries_circles"/
      "cells_boundaries_layers", indexed by row position instead, with the real cell_id kept as
      a *column*. read_xenium_sample prefixes both adata's obs_names and its obs["cell_id"]
      column with "proseg-" right after loading, but the shapes' cell_id column is never
      touched and stays raw -- so the prefix is stripped to recover the id the shapes were
      built with, used to look up rows in the shapes' cell_id column, and the matched shapes
      are then reindexed to obs[instance_key] (prefixed) so they line up with the table.

    Only covers per-cell shapes. Does not cover cell_labels/nucleus_labels/nucleus_boundaries,
    which use different indexing schemes and aren't currently loaded by this pipeline's calls
    to read_xenium_sample.
    """
    table_cell_ids = adata.obs[instance_key]

    if is_proseg:
        shape_keys = ("cells_boundaries", "cells_boundaries_circles", "cells_boundaries_layers")
        raw_cell_ids = table_cell_ids.astype(str).str.removeprefix("proseg-")
    else:
        shape_keys = ("cell_boundaries", "cell_circles")

    for key in shape_keys:
        if key not in sdata.shapes:
            continue

        if is_proseg:
            shp = sdata.shapes[key]
            matched = shp.set_index(shp["cell_id"].astype(str), drop=False).loc[raw_cell_ids]
            # realign to obs[instance_key] as-is (prefixed or not) rather than the raw id used
            # for the lookup, so the shapes' index matches the table's instance_key column
            matched.index = table_cell_ids.values
            sdata.shapes[key] = matched
        else:
            sdata.shapes[key] = sdata.shapes[key].loc[table_cell_ids]

        # set_table_annotates_spatialelement checks that obs[region_key] already agrees
        # with the target region before it will update the table's region metadata
        adata.obs["region"] = key
        sdata.set_table_annotates_spatialelement(
            table_name, region=key, region_key="region", instance_key=instance_key
        )
    return sdata


def read_xenium_sample(
    path,
    cells_as_circles=False,
    cells_boundaries=False,
    cells_boundaries_layers=False,
    nucleus_boundaries=False,
    cells_labels=False,
    nucleus_labels=False,
    transcripts=False,
    morphology_mip=False,
    morphology_focus=False,
    aligned_images=False,
    cells_table=True,
    gex_only=True,
    xr_segmentation="10x_mm_5um",
    mnn_max_dist=10,
    anndata=False,
    add_meta=False,
    celltypes=False,
    norm_id=None,
    ct_ref=None,
    ct_levels=None,
    count_correction=False,
    split_matrix=False,
    split_level=None,
    correction_mode = "split_fully_purified",
    rctd_type="rctd_class_aware",
    ct_results_file="results_df",
    full_post_processed_info=False,
    xeniumranger_dir=None,
    sample_name=None,
):
    """
    Reads a xenium sample from a directory path.

    Parameters
    ----------
    path (str): The directory path of the segmentation run. For xeniumranger segmentation is the normalised_results/outs folder, for proseg is the raw_results folder.
    cells_as_circles (bool): Whether to include cell polygons as circles or not.
    cells_boundaries (bool): Whether to include cell boundaries or not.
    nucleus_boundaries (bool): Whether to include nucleus boundaries or not.
    cells_labels (bool): Whether to include cell labels or not.
    nucleus_labels (bool): Whether to include nucleus labels or not.
    transcripts (bool): Whether to include transcript locations or not.
    morphology_mip (bool): Whether to include morphology MIP or not.
    morphology_focus (bool): Whether to include morphology focus or not.
    aligned_images (bool): Whether to include aligned images or not.
    cells_table (bool): Whether to include cells table or not.
    anndata (bool): Whether to return only the anndata object or the full spatialdata object.
    gex_only (bool): Whether to load only transcript information. Set to False to load protein data when available (Default = True).
    xr_segmentation (str): Segmentation ID for xenium ranger to get protein values. Default: "10x_mm_5um"
    mnn_max_dist (float): For proseg samples with gex_only=False, maximum centroid distance allowed between a proseg cell and its matched xenium ranger (protein) cell. Pairs matched further apart than this are dropped. Default: 10
    add_meta (bool): Whether to add sample information to .obs (sample_name, segmentation_id, condition, gene_panel, donor, sample).
    celltypes (bool): Whether to add cell type annotation to .obs dataframe. Requires norm_id, ct_ref, and ct_levels to be specified.
    norm_id (str): Normalization ID used for the cell type annotation results to be added. Required if celltypes is True.
    ct_ref (str): Name of the scRNAseq reference used for annotation in the pipeline. Required if celltypes is True.
    ct_levels (list): List indicating which cell type levels of granularity to add
    count_correction (bool): Whether cell type labels before or after count correction should be added. Default is False.
    split_matrix (bool): Whether to add SPLIT corrected count matrix to the anndata object. Default is False.
    split_level (str): Cell type level used to run SPLIT. Required if split_matrix is True.
    correction_mode (str): Correction mode used for SPLIT. Default="split_fully_purified"
    rctd_type (str, optional): Mode used when running RCTD. Default="rctd_class_aware"
    ct_results_file (str): Which RCTD output table to read cell type labels from. Options are "results_df" (raw RCTD output) or "post_processed_results_df" (SPLIT post-processed output). Default="results_df"
    full_post_processed_info (bool): Whether to also add the extra post-processing columns (second_type_class, weight_first_type, weight_second_type, singlet_score) per level. Requires ct_results_file="post_processed_results_df". Default is False.
    xeniumranger_dir (str): Path to xeniumranger output directory (for proseg raw only)
    sample_name (str): The sample name.
     

    Returns
    -------
    If anndata, returns a tuple of the sample name and anndata object.
    Otherwise, returns a tuple of the sample name and spatialdata object.

    If sample_name is None, sample_name is not returned
    """
    path = Path(path)
    kwargs = dict(
        cells_as_circles=cells_as_circles,
        cells_boundaries=cells_boundaries,
        nucleus_boundaries=nucleus_boundaries,
        cells_labels=cells_labels,
        nucleus_labels=nucleus_labels,
        transcripts=transcripts,
        morphology_mip=morphology_mip,
        morphology_focus=morphology_focus,
        aligned_images=aligned_images,
        cells_table=cells_table,
    )

    # automatically check whether path is a folder with proseg raw outputs or in xeniumranger format
    if (path / "expected-counts.csv.gz").exists():
        reader = xenium_proseg
        kwargs["cells_boundaries_layers"] = cells_boundaries_layers
        kwargs["xeniumranger_dir"] = xeniumranger_dir

    else:
        reader = spatialdata_io.xenium

        # include the possibility of loading protein data if available if xenium ranger segmentation is selected
        kwargs["gex_only"] = gex_only

    sdata = reader(path, **kwargs)

    adata = sdata["table"]

    # add spatial coordinates to .obs
    adata.obs["x_centroid"] = adata.obsm["spatial"][:,0]
    adata.obs["y_centroid"] = adata.obsm["spatial"][:,1]

    # add obs_names to anndata
    is_proseg = "proseg" in str(path)
    if is_proseg:
        # prefix both the cell_id column and obs_names so they stay consistent for the rest
        # of the pipeline (shapes in sdata keep the raw, unprefixed proseg cell id though)
        adata.obs["cell_id"] = "proseg-" + adata.obs["cell_id"].astype(str)
        adata.obs_names = adata.obs["cell_id"].values

    else:
        adata.obs_names = adata.obs["cell_id"].values

    # extract sample information from path
    trim_path = re.sub(".+segmentation/", "", str(path))

    # extract metadata information from path
    segmentation = trim_path.split("/")[0]

    # sometimes is written proseg and sometimes proseg_expected...
    if segmentation == "proseg":

        seg_id="proseg_expected"

    else:
        seg_id = segmentation

    condition = trim_path.split("/")[1]
    panel = trim_path.split("/")[2]
    donor = trim_path.split("/")[3]
    sample = trim_path.split("/")[4]

    # add metadata information to adata.obs
    if add_meta:

        adata = add_metadata2ad(adata, 
                                sample_id=sample,
                                segmentation_id=seg_id,
                                condition=condition,
                                gene_panel=panel,
                                donor=donor,
                                sample=sample,
                                segmentation_method=segmentation)
        
    # cell type information if specified
    if celltypes:

        full_sample_id = f"{condition}/{panel}/{donor}/{sample}"

        xp_dir = re.sub("/segmentation.+$", "", str(path))

        adata = add_ct_labels(adata=adata, 
                              xp_dir=xp_dir, 
                              seg_id=seg_id, 
                              full_sample_id=full_sample_id,
                              norm_id = norm_id, 
                              ct_ref = ct_ref,
                              ct_levels = ct_levels,
                              count_correction = count_correction,
                              correction_mode = correction_mode,
                              rctd_type = rctd_type,
                              ct_results_file = ct_results_file,
                              full_post_processed_info = full_post_processed_info)

    # split anndata into RNA and protein if specified and xenium ranger output is loaded
    matched_by_mnn = False
    if not gex_only and seg_id != "proseg_expected":

        adata_prot = adata[:, adata.var["feature_types"] == "Protein Expression"].copy()

        adata = adata[:, adata.var["feature_types"] == "Gene Expression"].copy()

    # in case we want proseg corrected counts and protein data
    # since protein quantification is done on xenium ranger segmentation, we need to map proseg and xr cells
    # using spatial mutual nearest neighbours
    elif not gex_only and seg_id == "proseg_expected":

        print(f"Loading data from xenium ranger {xr_segmentation} to get protein intensities")

        # get the sibling xenium ranger (protein) path for this proseg sample
        xr_path = re.sub(r"/segmentation/proseg/", f"/segmentation/{xr_segmentation}/", str(path))
        xr_path = Path(re.sub(r"raw_results$", "normalised_results/outs", xr_path))

        # load anndata with protein information (minimal read, only need the table)
        sdata_xr = spatialdata_io.xenium(xr_path, cells_table=True, gex_only=False)
        adata_xr = sdata_xr["table"]

        # add spatial coordinates to .obs
        adata_xr.obs["x_centroid"] = adata_xr.obsm["spatial"][:, 0]
        adata_xr.obs["y_centroid"] = adata_xr.obsm["spatial"][:, 1]
        adata_xr.obs_names = adata_xr.obs["cell_id"].values

        # subset for the protein features
        adata_prot = adata_xr[:, adata_xr.var["feature_types"] == "Protein Expression"].copy()

        print("Spatial MNN matching between Proseg and XR cells")

        # match proseg (adata) and xenium ranger (adata_prot) cells by spatial mutual nearest neighbours
        matches = match_cells_spatial_mnn(adata.obs, adata_prot.obs, max_dist=mnn_max_dist)

        # subset to the mutually-matched cells (aligned by row order)
        adata = adata[matches["proseg_cell"].values, :].copy()
        adata_prot = adata_prot[matches["xr_cell"].values, :].copy()

        # keep each side's original cell ID for traceability, before aligning obs_names
        adata.obs["xr_cell_id"] = matches["xr_cell"].values
        adata_prot.obs["proseg_cell_id"] = matches["proseg_cell"].values

        # align obs_names so adata and adata_prot share the same cell index
        adata_prot.obs_names = adata.obs_names

        matched_by_mnn = True

    # add cell type labels to .obs of protein anndata
    if celltypes and not gex_only:

        # make dataframe with shared cell IDs and cell type levels (including spot_class and any extra post-processing columns)
        ct_cols = [c for c in adata.obs.columns if any(c == f'{lvl}_{ct_ref}' or c.startswith(re.sub(r"(\d+).*$", r"\1", lvl) + "_") for lvl in ct_levels)]
        ct_df = adata.obs[ct_cols].copy()
        ct_df["cellID"] = adata.obs_names

        # left join cell type annotations to protein anndata
        adata_prot.obs = adata_prot.obs.merge(ct_df, left_index=True, right_on="cellID", how="left").set_index(adata_prot.obs_names)

        # drop cellID column
        adata_prot.obs = adata_prot.obs.drop(columns=["cellID"])

    # make sure names are unique if protein data is loaded
    if not gex_only:
        adata_prot.var_names = [f"prot_{i}" for i in adata_prot.var_names]


    # add post-SPLIT count matrix
    if split_matrix:

        full_sample_id = f"{condition}/{panel}/{donor}/{sample}"

        xp_dir = re.sub("/segmentation.+$", "", str(path))

        adata = add_split_matrix(adata=adata, 
                                 xp_dir=xp_dir, 
                                 seg_id=seg_id, 
                                 full_sample_id=full_sample_id,
                                 norm_id = norm_id, 
                                 ct_ref = ct_ref,
                                 ct_level = split_level,
                                 correction_mode = correction_mode,
                                 rctd_type = rctd_type)

        # for some reason it makes a copy of the adata inside the add_split_matrix function...
        sdata["table"] = adata
        

    # add sample level summary stats to .uns
    metrics_path = path / "metrics_summary.csv"
    if metrics_path.exists():
        adata.uns["metrics_summary"] = pd.read_csv(metrics_path)
    else:
        print("metrics_summary.csv not found at:", metrics_path)

    # compute basic qc metrics for the sample of interest (e.g.transcripts and number of genes detected per cell)
    # create a boolean for all the genes
    adata.var["all_genes"] = np.ones(adata.shape[1], dtype=bool)

    # mitochondrial genes, "MT-" for human, "Mt-" for mouse
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))

    # compute total number of genes detected per cell and total estimated transcripts/cell
    sc.pp.calculate_qc_metrics(adata, qc_vars=["all_genes", "mt", "ribo"], inplace=True, percent_top=[50])

    if sample_name is None:

        if gex_only:

            if anndata:
                return adata
            else:
                sdata = _sync_shapes_to_table(sdata, "table", adata, is_proseg=is_proseg)
                return sdata
        # return muon object with RNA and protein adata objects
        else:

            # make sure that RNA and protein objects have the same cells.
            cells = adata.obs_names.intersection(adata_prot.obs_names)
            adata = adata[cells, :].copy()
            adata_prot = adata_prot[cells, :].copy()

            if anndata:
                mdata = mu.MuData({"rna": adata, "prot": adata_prot})
                return mdata
            # add the anndata objects as separate entries in sdata["tables"]
            else:
                sdata["table_rna"] = adata
                sdata["table_prot"] = adata_prot
                del sdata["table"]
                sdata = _sync_shapes_to_table(sdata, "table_rna", adata, is_proseg=is_proseg)
                return sdata
    else:
        if anndata:
            return sample_name, adata
        else:
            sdata = _sync_shapes_to_table(sdata, "table", adata, is_proseg=is_proseg)
            return sample_name, sdata


######### Add metadata to anndata
def add_metadata2ad(
    adata: ad.AnnData,
    sample_id: str,
    segmentation_id: str,
    condition: str,
    gene_panel: str,
    donor: str,
    sample: str,
    segmentation_method: str,
) -> ad.AnnData:
    """
    Adds metadata to an AnnData object.

    Parameters
    ----------
    adata : ad.AnnData
        The AnnData object to add metadata to.
    sample_id : str
        The sample ID.
    segmentation_id : str
        The segmentation ID.
    condition : str
        The condition.
    gene_panel : str
        The gene panel.
    donor : str
        The donor.
    sample : str
        The sample.
    segmentation_method : str
        The segmentation method.

    Returns
    -------
    ad.AnnData
        The AnnData object with metadata added.
    """
    adata.obs["sample_id"] = sample_id
    adata.obs["segmentation_id"] = segmentation_id
    adata.obs["condition"] = condition
    adata.obs["gene_panel"] = gene_panel
    adata.obs["donor"] = donor
    adata.obs["sample"] = sample
    adata.obs["segmentation_method"] = segmentation_method

    return adata
