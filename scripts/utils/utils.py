from affine import Affine
import anndata as ad
import geopandas as gpd
import glob
import json
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
from shapely.geometry import shape
from sklearn.neighbors import NearestNeighbors
import tifffile
import xml.etree.ElementTree as ET

def get_channels(image_path):
    '''
    Extracts the channel names from the OME-XML metadata of a
    multi-channel TIFF image.
    
    Parameters:
    - image_path: str, path to the input multi-channel TIF image.
    '''
    with tifffile.TiffFile(image_path) as tif:
        print("File opened successfully.")
        # OME metadata is stored in the description tag of the first page
        xml_data = tif.pages[0].description
        # Parse the XML
        root = ET.fromstring(xml_data)
        # OME-XML uses namespaces, we need to account for them
        namespace = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06'}
        # Find all Channel elements and get their 'Name' attribute
        # it is crucial that the channel names order is equal to how the channels appear in the image 
        # the reference channel is taken not by name but by index
        channels = root.findall('.//ome:Channel', namespace)
        channel_names = [c.attrib.get('Name') for c in channels]
    return channel_names

def spatial_to_adata(csv_file, hid, rid, method):
    df = pd.read_csv(csv_file)
    id_info = pd.DataFrame()
    id_info["rsamp_id"] = rid
    id_info["samp_id"] = hid
    if method == "mesmer":
        df.columns = ["Cells/Mean Intensity " + col if "autof" in col else col for col in df.columns]
        area_col = "area"
        x_col = "x"
        y_col = "y"
    elif method == "horizon":
        area_col = "Cells/Area in μm²"
        x_col = "Cells/X Position in μm"
        y_col = "Cells/Y Position in μm"
    else:
        raise ValueError("Method must be either 'mesmer' or 'horizon'")
    
    obs_df = pd.DataFrame.from_dict(dict(area = df[area_col],
                                        samp_id = hid,
                                        cell_id = [f'{method}-{i}' for i in range(len(df[area_col]))],
                                        x_centroid = df[x_col],
                                        y_centroid = df[y_col]))
    
    obs_df = obs_df.merge(id_info, left_on="samp_id", how="left", right_on="samp_id")
    obs_df["library_id"] = hid
    feat_df = df[[col for col in df.columns if "Intensity" in col and "elu" not in col]]
    cols = [re.sub(".+sity \(", "", i) for i in feat_df.columns]
    cols = [re.sub(" -.+$", "", i) for i in cols]
    cols = [re.sub(" .+|_.+", "", i) for i in cols]
    cols = [re.sub("(α)", "alpha-", i) for i in cols]
    cols = [re.sub("\)", "", i) for i in cols]
    feat_df.columns = cols

    adata = ad.AnnData(X=feat_df, obs=obs_df, obsm=dict(spatial=obs_df[["x_centroid", "y_centroid"]].to_numpy()))
    adata.uns["spatial"] = {"hsamp_id": {}}
    adata.layers["raw"] = adata.X.copy()
    return adata

# function to apply affine transformation matrix to cell centroids
def transf_centroids(adata, 
                     aff_mtx,
                     cell_id_col,
                     mpp_input=None,
                     mpp_aligned=None,
                     centroids_microns=True,
                     x_col="x_centroid", 
                     y_col="y_centroid",
                     crop_gpd=None,
                     add_to_obs=True):
    
    """

    Function to apply an affine transformation the cell centroids. It requires an anndata object with cell centroids stored in .obs and an affine transformation matrix.
    So far I only tested with the automated alignment using Insitupy.

    Arguments
    ------------

    adata: class, AnnData
        anndata object with cell centroids stored in the .obs dataframe
    aff_mtx: str
        path to .csv file containing the affine matrix computed after image alignment
    cell_id_col: str
        column name which contains the unique cell IDs in .obs
    mpp_input: float,optional
        Microns per pixel of the input image. Needed in case cell centroids given in pixels
    mpp_aligned: float,optional
        Microns per pixel of the input image. Needed in case cell centroids given in pixels 
    centroid_microns: boolean
        Boolean indicating if cell centroids are stored in microns. If false it is assumed to be pixels (Default=True)
    x_col: str, optional
        Column name in .obs which stores the x coordinate of cell centroid (Default='x_centroid')
    y_col: str, optional
        Column name in .obs which stores the y coordinate of cell centroids (Default='y_centroid')
    crop_gpd: str, optional
        Path to geojson file which stores the polygon used to crop the input image before alignment. Only needed if cell segmentation was done in the uncropped image.
    add_to_obs: boolean
        Boolean indicating if centered centroids should be added to the .obs slot. False will return a dataframe

    Returns 
    -------------

    anndata: class, AnnData
        Returns the anndata object with affine transformed centroids as columns x_centroid_aligned in .obs. If add_to_obs=False, it returns a dataframe

    """

    # subset .obs for columns of interest
    centr_df = adata.obs.loc[:,[cell_id_col, x_col,y_col]]
    # if cropping shape file is specified loaded and put the coordinates as dataframe
    # for now I assume that the shape is a rectangle
    if crop_gpd is not None:
        crop_gdf=gpd.read_file(crop_gpd)
        crop_coords=pd.DataFrame([pt for geom in crop_gdf.geometry for pt in geom.exterior.coords],
                                columns=["x", "y"])
        # if centroids are given in microns, I need to transform the rectangle coordinates from pixels to um
        if centroids_microns:
            crop_coords = crop_coords*mpp_input
        # recenter the cell centroids taking the bottom left corner of the cropping rectangle as origin
        centr_df[x_col] = centr_df[x_col]-crop_coords["x"].min()
        centr_df[y_col] = centr_df[y_col]-crop_coords["y"].min()
    # Load matrix for affine tranformation
    # same approach Bernd use to align transcripts

    aff_mtx_np = np.loadtxt(aff_mtx, delimiter=",")
    aff_mtx = Affine(*aff_mtx_np.flatten())
    # extract cell centroids
    xy = centr_df[[x_col, y_col]].values
    # tranform them to pixels if specified
    if centroids_microns:
        xy = xy / mpp_input
    # apply affine transformation
    xy_aligned = np.apply_along_axis(lambda x: aff_mtx * x, 1, xy)
    # transform them back to microns if specified
    if centroids_microns:
        xy_aligned = xy_aligned*mpp_aligned
    # add align coordinates to data frame
    centr_df["x_centroid_aligned"] = xy_aligned[:,0]
    centr_df["y_centroid_aligned"] = xy_aligned[:,1]
    # add to algined coordinates to .obs if specified

    if add_to_obs:
        adata.obs = adata.obs.merge(centr_df, how="left", left_on=cell_id_col, right_on=cell_id_col)
        return adata
    else:
        return centr_df

def get_coords(geojson_file):
    with open(geojson_file) as f:
        geojson = json.load(f)
    bounds = shape(geojson["features"][0]["geometry"]).bounds
    return bounds

def uniquify(df: pd.DataFrame) -> pd.DataFrame:
    if df.shape[0] == 1:
        return df
    
    ret = df.sort_values("dist")
    return ret.iloc[[0]]

def get_nearest_neighbors(
    ref_df: pd.DataFrame,
    inquery_df: pd.DataFrame,
    ref_cell_col_name: str,
    inquery_cell_col_name: str,
    *,
    x_centroid: str = "x_centroid",
    y_centroid: str = "y_centroid",
    cell_id: str = "cell_id",
) -> pd.DataFrame:
    # find closest nearest neighbour in spatial space
    dist, idx = (
        NearestNeighbors(n_neighbors=1)
        .fit(ref_df[[x_centroid, y_centroid]])
        .kneighbors(inquery_df[[x_centroid, y_centroid]])
    )

    ids = [x[0] for x in idx]
    dists = [x[0] for x in dist]
    return pd.DataFrame(
            {
                ref_cell_col_name: ref_df[cell_id].iloc[ids].values,
                inquery_cell_col_name: inquery_df[cell_id].values,
                "dist": dists,
            }).groupby(ref_cell_col_name,).apply(uniquify,).reset_index(drop=True,) # type: ignore

def map_metadata(adata_comet, adata_xe, ct_cols=None, max_dist=5):

    """
    Function to add cell type labels to comet.obs file and also to add comet cell IDs to xenium anndata object

    Arguments
    ------------

    adata_comet: class, AnnData
        anndata object containing the comet data. It should contain the columns "x_centroid_aligned" and "y_centroid_aligned" in .obs obtained after running transf_centroids function
    adata_xe: class, AnnData
        anndata object containing the xenium data. Cell centroids should be stored in .obs as "x_centroid" and "y_centroid" columns
    ct_cols: list
        list of column names in .obs Xenium anndata object which contain cell type labels to be mapped to the comet anndata object
    max_dist: float, optional
        Maximum distance in microns between comet and xenium cell centroids to be considered a match. Default: 5

    Returns
    ------------
    tuple: (adata_comet, adata_xe)
        return a tuple with the two input anndata objects with cell type labels in comet and comet cell IDs in xenium .obs.       

    """
    # make copies of anndata objects to avoid inplace modifications
    adata_comet = adata_comet.copy()
    adata_xe = adata_xe.copy()
    comet_obs = adata_comet.obs.copy()
    xe_obs = adata_xe.obs.copy()

    # the cell centroid columns must have the same name
    comet_obs["x_centroid"] = comet_obs["x_centroid_aligned"]
    comet_obs["y_centroid"] = comet_obs["y_centroid_aligned"]
    xe_obs["x_centroid"] = xe_obs["centroid_x"]
    xe_obs["y_centroid"] = xe_obs["centroid_y"]

    # NN using the xenium data as reference and comet as query
    comet2xr = get_nearest_neighbors(
        ref_df=xe_obs,
        inquery_df=comet_obs,
        ref_cell_col_name= "xe_cell_id",
        inquery_cell_col_name= "comet_cell_id",
    )
    # do the opposite
    xr2comet = get_nearest_neighbors(
        ref_df=comet_obs,
        inquery_df=xe_obs,
        ref_cell_col_name= "comet_cell_id",
        inquery_cell_col_name= "xe_cell_id",
    )
    # df with mapped cells
    mapped = pd.merge(
        comet2xr,
        xr2comet,
        how="inner",
        on=["xe_cell_id", "comet_cell_id"],
        suffixes=("_comet2xe", "_xe2comet"),
    )
    # merge Xenium cell type labels to mapped dataframe
    mrgd_df = pd.merge(mapped, adata_xe.obs.loc[:,["cell_id"] + ct_cols].reset_index(drop=True), 
                    how="left", left_on="xe_cell_id", right_on="cell_id").drop(columns=["cell_id"])
    
    # add Xenium cell type labels to adata_comet.obs 
    mrgd_df = pd.merge(comet_obs, mrgd_df[mrgd_df["dist_comet2xe"] <= max_dist], 
                        how="left", left_on="cell_id", right_on="comet_cell_id").drop(columns=["cell_id"])
    
    # left join the dataframes
    mrgd_df = pd.merge(adata_comet.obs, mrgd_df.loc[:,(["comet_cell_id", "xe_cell_id"] + ct_cols)], 
                    how="left", left_on="cell_id", right_on="comet_cell_id").drop(columns=["comet_cell_id"])
    
    # add information to the .obs dataframe of the comet anndata object
    adata_comet.obs = mrgd_df
    # add column indicating if the cell was matched or not
    adata_comet.obs["matched"] = ~adata_comet.obs["xe_cell_id"].isna()
    # add comet cell IDs to adata_xe
    # add columns to the Xenium anndata .obs with the comet cell ID and the marker intensities
    adata_xe.obs = pd.merge(adata_xe.obs, 
                    mrgd_df.loc[:,["xe_cell_id", "cell_id"]].rename(columns={"cell_id": "comet_cell_id"}), 
                    how="left", left_on="cell_id", right_on="xe_cell_id").drop(columns=["xe_cell_id"])
    adata_xe.obs["matched"] = ~adata_xe.obs["comet_cell_id"].isna()
    return adata_comet, adata_xe