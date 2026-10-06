"""Proteomics + Xenium readers shared with the analysis_repo repository.

Copied verbatim from analysis_repo/workflow/scripts/general/multimodal_readwrite.py (uncommitted
there as of 2026-10-06), so that proteomics cells are aligned and matched to Xenium cells the same
way in both projects. Do not edit it here: change it there and copy it over again. The only change is
the import below, from xenium_readwrite (this repo's copy of readwrite.py) instead of readwrite.
"""

# read/write helpers to load a Xenium sample together with the spatial proteomics (COMET, CODEX)
# measured on the same section, as processed by the spatial_proteomics_pipeline:
#   processed/segmentation/{seg_method}/{sample_id}/feature_table.csv + *_boundaries.geojson
#   processed/xenium_alignment/{align_method}/max_width_{N}/{sample_id}/affine_matrix.csv
# The affine matrix maps proteomics image pixels onto Xenium morphology pixels.
# readwrite.py (same folder) must be importable, i.e. on sys.path, before importing this module.

import os
import re

import numpy as np
import pandas as pd
import geopandas as gpd
import anndata as ad
import muon as mu
from spatialdata.models import ShapesModel, TableModel
from spatialdata.transformations import Affine, Scale, Sequence, get_transformation

from xenium_readwrite import read_xenium_sample, read_xenium_specs, match_cells_spatial_mnn, _sync_shapes_to_table


# feature table columns that are cell metadata rather than marker intensities
MESMER_META_COLS = ["cell_id", "x_px", "y_px", "area_px", "x_um", "y_um", "area_um2"]

# segmentation boundary files written by the spatial_proteomics_pipeline, keyed by element suffix
MESMER_BOUNDARY_FILES = {"cell_boundaries": "cell_boundaries.geojson",
                         "nucleus_boundaries": "nuclei_boundaries.geojson",
                         "cytoplasm_boundaries": "cytoplasm_boundaries.geojson"}

# read_xenium_sample flags that add spatial elements (shapes, points, images, labels)
XENIUM_ELEMENT_FLAGS = ["cells_as_circles", "cells_boundaries", "cells_boundaries_layers", "nucleus_boundaries",
                        "cells_labels", "nucleus_labels", "transcripts", "morphology_mip", "morphology_focus",
                        "aligned_images"]


def proteomics_sample_paths(processed_dir,
                            sample_id,
                            seg_method="mesmer",
                            align_method="insitupy",
                            max_width=8000):
    """
    Paths of the segmentation folder and of the affine matrix of one proteomics sample processed with
    the spatial_proteomics_pipeline.

    Arguments
    ---------
    processed_dir: str
        processed/ folder of the pipeline output (e.g. cfg["codex_cx_dir"] + "processed/")
    sample_id: str
        sample ID as {tissue}/{panel}/{donor}/{sample}
    seg_method: str
        segmentation method. Default: "mesmer"
    align_method: str
        method used to align the proteomics image onto the Xenium image. Default: "insitupy"
    max_width: int
        max_width the insitupy alignment was run with. Default: 8000

    Returns
    -------
    dict with keys "seg_dir" and "affine"
    """
    return dict(seg_dir=os.path.join(processed_dir, "segmentation", seg_method, sample_id),
                affine=os.path.join(processed_dir, "xenium_alignment", align_method, f"max_width_{max_width}",
                                    sample_id, "affine_matrix.csv"))


# function to find the Xenium pixel size (µm/px) of a sample
def xenium_pixel_size(xenium_path, xeniumranger_dir=None):
    """
    Look for experiment.xenium in the Xenium segmentation folder, in xeniumranger_dir, or, for proseg
    runs, in the matching raw Xenium output folder, and return its pixel size in µm/px.
    """
    # candidate locations of experiment.xenium
    xenium_path = str(xenium_path)
    candidates = [os.path.join(xenium_path, "experiment.xenium")]
    if xeniumranger_dir is not None:
        candidates.append(os.path.join(str(xeniumranger_dir), "experiment.xenium"))
    raw_run = re.sub(r"/raw_results/?$", "", re.sub(r"/processed/segmentation/[^/]+/", "/raw/", xenium_path))
    candidates.append(os.path.join(raw_run, "experiment.xenium"))

    for c in candidates:
        if os.path.exists(c):
            return read_xenium_specs(c)["pixel_size"]

    raise FileNotFoundError(f"experiment.xenium not found in {candidates}; pass xenium_pixel_size explicitly")


# function to apply a 3x3 affine matrix to (N, 2) x, y coordinates
def apply_affine(xy, matrix):
    xy = np.asarray(xy, dtype=float)
    homogeneous = np.column_stack([xy, np.ones(len(xy))])
    return (homogeneous @ np.asarray(matrix, dtype=float).T)[:, :2]


def mesmer_feature_table_to_adata(feature_table, dapi_to_obs=True):
    """
    Load the per-cell mean marker intensities of a Mesmer segmentation (spatial_proteomics_pipeline
    quantify_cells output) into an AnnData object.

    Arguments
    ---------
    feature_table: str
        path to feature_table.csv (cell_id, x_px, y_px, area_px, x_um, y_um, area_um2, one column per marker)
    dapi_to_obs: bool
        whether to move the DAPI channel from the marker matrix to .obs. Default: True

    Returns
    -------
    AnnData with marker intensities in .X (and layers["raw"]), cell metadata in .obs indexed by the
    Mesmer cell label (as str), and the proteomics image coordinates (µm) in .obsm["spatial"]
    """
    df = pd.read_csv(feature_table)

    # split cell metadata and marker intensities
    obs = df[MESMER_META_COLS].copy()
    markers = [c for c in df.columns if c not in MESMER_META_COLS]

    # move DAPI to .obs, it is not a marker
    if dapi_to_obs and "DAPI" in markers:
        obs["DAPI"] = df["DAPI"].to_numpy()
        markers.remove("DAPI")

    # build AnnData indexed by the Mesmer cell label
    obs.index = obs["cell_id"].astype(str).values
    adata = ad.AnnData(X=df[markers].to_numpy(dtype=np.float32),
                       obs=obs,
                       var=pd.DataFrame(index=markers))
    adata.layers["raw"] = adata.X.copy()
    adata.obsm["spatial"] = obs[["x_um", "y_um"]].to_numpy()

    return adata


def read_mesmer_boundaries(geojson, cell_ids=None):
    """
    Read one Mesmer boundary geojson (polygons in proteomics image pixels) into a GeoDataFrame indexed
    by the Mesmer cell label. If cell_ids is given, only those cells are kept.
    """
    # read polygons and index them by the Mesmer cell label
    gdf = gpd.read_file(geojson)
    gdf.index = gdf["label"].astype(int).values
    gdf = gdf[["geometry"]]

    # keep only the requested cells
    if cell_ids is not None:
        gdf = gdf.loc[gdf.index.isin(np.asarray(cell_ids, dtype=int))]

    return gdf


# function to get the transformation from Xenium µm to the "global" coordinate system of a Xenium sdata
def _xenium_um_to_global(sdata, pixel_size):
    # take it from a Xenium shapes or points element, which are stored in µm
    for elements in (sdata.shapes, sdata.points):
        for name in elements:
            return get_transformation(elements[name], "global")

    # spatialdata_io convention otherwise: "global" is in Xenium morphology pixels
    return Scale([1.0 / pixel_size, 1.0 / pixel_size], axes=("x", "y"))


def read_xenium_proteomics_sample(xenium_path,
                                  seg_dir,
                                  affine,
                                  seg_method="mesmer",
                                  max_dist=5,
                                  protein_cell_boundaries=False,
                                  protein_nucleus_boundaries=False,
                                  protein_cytoplasm_boundaries=False,
                                  dapi_to_obs=True,
                                  xenium_pixel_size_um=None,
                                  **xenium_kwargs):
    """
    Load a Xenium sample together with the spatial proteomics (COMET, CODEX) segmentation of the same
    section, bring the proteomics cells into Xenium coordinates with the image alignment affine matrix,
    and match them to the Xenium cells by spatial mutual nearest neighbours.

    Returns a MuData object ("rna", "prot") when no spatial element is requested, or a SpatialData
    object when any Xenium element flag (cells_boundaries, nucleus_boundaries, transcripts,
    morphology_focus, ...) or any protein_*_boundaries flag is True.

    Only cells matched in both modalities are kept, and both modalities share the Xenium cell IDs as
    obs_names (as for the proseg/xenium ranger matching in read_xenium_sample).

    Arguments
    ---------
    xenium_path: str
        Xenium segmentation run passed to read_xenium_sample (normalised_results/outs for xenium ranger,
        raw_results for proseg)
    seg_dir: str
        folder of the proteomics segmentation with feature_table.csv and the *_boundaries.geojson files
        (see proteomics_sample_paths)
    affine: str
        3x3 affine matrix (.csv) mapping proteomics image pixels onto Xenium morphology pixels
    seg_method: str
        proteomics segmentation method. Only "mesmer" is supported for now. Default: "mesmer"
    max_dist: float
        maximum distance (µm) between a Xenium cell and its matched proteomics cell. Default: 5
    protein_cell_boundaries, protein_nucleus_boundaries, protein_cytoplasm_boundaries: bool
        whether to load the proteomics cell / nucleus / cytoplasm polygons as shapes in Xenium
        coordinates (SpatialData output). Default: False
    dapi_to_obs: bool
        whether to move the proteomics DAPI channel to .obs. Default: True
    xenium_pixel_size_um: float, optional
        Xenium pixel size (µm/px). By default read from experiment.xenium (see xenium_pixel_size)
    **xenium_kwargs:
        any other read_xenium_sample argument (cells_boundaries, transcripts, add_meta, celltypes,
        ct_ref, ct_levels, xeniumranger_dir, ...). anndata, gex_only and sample_name are set here.

    Returns
    -------
    MuData or SpatialData. In SpatialData, "table_rna" annotates the Xenium cell shapes and
    "table_prot" annotates the {seg_method}_cell_boundaries shapes (when loaded) through the
    protein_cell_id column.
    """
    if seg_method != "mesmer":
        raise ValueError(f"segmentation method {seg_method} not supported, only 'mesmer'")
    for k in ("anndata", "gex_only", "sample_name"):
        if k in xenium_kwargs:
            raise ValueError(f"{k} is set by read_xenium_proteomics_sample, do not pass it")

    # decide the output type from the requested spatial elements
    protein_shapes = {"cell_boundaries": protein_cell_boundaries,
                      "nucleus_boundaries": protein_nucleus_boundaries,
                      "cytoplasm_boundaries": protein_cytoplasm_boundaries}
    as_sdata = any(xenium_kwargs.get(f, False) for f in XENIUM_ELEMENT_FLAGS) or any(protein_shapes.values())

    # Xenium RNA (AnnData, or SpatialData with the table in "table")
    xe = read_xenium_sample(xenium_path, anndata=not as_sdata, gex_only=True, **xenium_kwargs)
    adata_rna = xe["table"] if as_sdata else xe

    # Xenium pixel size to convert aligned pixels to µm
    if xenium_pixel_size_um is None:
        xenium_pixel_size_um = xenium_pixel_size(xenium_path, xenium_kwargs.get("xeniumranger_dir"))

    # proteomics cells with their centroids aligned to Xenium µm
    adata_prot = mesmer_feature_table_to_adata(os.path.join(seg_dir, "feature_table.csv"), dapi_to_obs=dapi_to_obs)
    matrix = np.loadtxt(affine, delimiter=",")
    xy_aligned = apply_affine(adata_prot.obs[["x_px", "y_px"]].to_numpy(), matrix) * xenium_pixel_size_um
    adata_prot.obs["x_centroid_aligned"] = xy_aligned[:, 0]
    adata_prot.obs["y_centroid_aligned"] = xy_aligned[:, 1]
    print(f"{adata_rna.n_obs} Xenium cells, {adata_prot.n_obs} {seg_method} cells")

    # mutual nearest neighbours between Xenium and aligned proteomics cells
    prot_xy = adata_prot.obs[["x_centroid_aligned", "y_centroid_aligned"]].set_axis(["x_centroid", "y_centroid"], axis=1)
    matches = match_cells_spatial_mnn(adata_rna.obs, prot_xy, max_dist=max_dist)
    matches = matches.rename(columns={"proseg_cell": "xenium_cell_id", "xr_cell": "protein_cell_id"})
    print(f"{len(matches)} cells matched within {max_dist} µm")

    # subset both modalities to the matched pairs (aligned by row order)
    adata_rna = adata_rna[matches["xenium_cell_id"].values, :].copy()
    adata_prot = adata_prot[matches["protein_cell_id"].values, :].copy()

    # keep each side's original cell ID and the match distance
    adata_rna.obs["protein_cell_id"] = adata_prot.obs["cell_id"].to_numpy()
    adata_rna.obs["match_dist"] = matches["dist"].to_numpy()
    adata_prot.obs = adata_prot.obs.rename(columns={"cell_id": "protein_cell_id"})
    adata_prot.obs["xenium_cell_id"] = matches["xenium_cell_id"].to_numpy()
    adata_prot.obs["match_dist"] = matches["dist"].to_numpy()

    # proteomics image coordinates kept apart, aligned Xenium µm become the spatial coordinates
    adata_prot.obsm["spatial_px"] = adata_prot.obs[["x_px", "y_px"]].to_numpy()
    adata_prot.obsm["spatial"] = adata_prot.obs[["x_centroid_aligned", "y_centroid_aligned"]].to_numpy()

    # share the Xenium cell IDs as obs_names and keep protein names unique
    adata_prot.obs_names = adata_rna.obs_names
    adata_prot.var_names = [f"prot_{i}" for i in adata_prot.var_names]

    if not as_sdata:
        return mu.MuData({"rna": adata_rna, "prot": adata_prot})

    # --- SpatialData output ---
    sdata = xe
    is_proseg = "proseg" in str(xenium_path)  # same rule as read_xenium_sample

    # Xenium table subset to the matched cells, with its shapes synced to it
    del sdata["table"]
    sdata["table_rna"] = adata_rna
    sdata = _sync_shapes_to_table(sdata, "table_rna", adata_rna, is_proseg=is_proseg)

    # proteomics pixels -> Xenium pixels -> Xenium µm -> "global" of the Xenium sdata
    to_global = Sequence([Affine(matrix, input_axes=("x", "y"), output_axes=("x", "y")),
                          Scale([xenium_pixel_size_um, xenium_pixel_size_um], axes=("x", "y")),
                          _xenium_um_to_global(sdata, xenium_pixel_size_um)])

    # proteomics boundaries of the matched cells
    prot_ids = adata_prot.obs["protein_cell_id"].to_numpy()
    for key, load in protein_shapes.items():
        if not load:
            continue
        gdf = read_mesmer_boundaries(os.path.join(seg_dir, MESMER_BOUNDARY_FILES[key]), cell_ids=prot_ids)
        sdata.shapes[f"{seg_method}_{key}"] = ShapesModel.parse(gdf, transformations={"global": to_global})
        print(f"loaded {len(gdf)} {seg_method} {key}")

    # proteomics table, annotating the proteomics cell shapes if they were loaded
    region = f"{seg_method}_cell_boundaries"
    if region in sdata.shapes:
        adata_prot.obs["region"] = pd.Categorical([region] * adata_prot.n_obs)
        sdata["table_prot"] = TableModel.parse(adata_prot, region=region, region_key="region",
                                               instance_key="protein_cell_id")
    else:
        sdata["table_prot"] = TableModel.parse(adata_prot)

    return sdata
