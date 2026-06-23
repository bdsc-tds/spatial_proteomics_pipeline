import argparse
import geopandas as gpd
import glob
from insitupy.io import read_xenium
from insitupy.tools import register_images
from ome_types import from_xml, to_xml
import os
import pandas as pd
import tifffile
from shapely.geometry import box
import shutil
import sys
import xml.etree.ElementTree as ET

def crop_image(image_path, output_dir, crop_coords):
    """
    Crops a region from a multi-channel TIFF image and saves it as a
    new TIF file along with a GeoJSON file containing the cropping box.

    Parameters:
    - image_path: str, path to the input multi-channel TIFF image.
    - output_dir: str, directory to save the cropped image and GeoJSON file.
    - crop_coords: tuple, (xmin, ymin, xmax, ymax) coordinates for the cropping region.
    """
    xmin, ymin, xmax, ymax = crop_coords
    # Read the multi-channel TIFF image
    with tifffile.TiffFile(image_path) as tif:
        ome_xml = tif.ome_metadata
        data = tif.asarray()
    print("Loaded file")
    # Crop the image using the specified coordinates
    cropped_img = data[:, ymin:ymax, xmin:xmax]
    '''
    new_y = ymax - ymin
    new_x = xmax - xmin

    if ome_xml:
        try:
            ome = from_xml(ome_xml)
            ome.images[0].pixels.size_x = new_x
            ome.images[0].pixels.size_y = new_y
            ome_xml = to_xml(ome)
            greek_map = {
                'α': 'alpha', 'β': 'beta', 'γ': 'gamma', 'δ': 'delta',
                'ε': 'epsilon', 'κ': 'kappa', 'λ': 'lambda', 'μ': 'mu',
                'µ': 'mu', 'π': 'pi', 'σ': 'sigma', 'τ': 'tau',
            }
            for greek, latin in greek_map.items():
                ome_xml = ome_xml.replace(greek, latin)
        except Exception as e:
            print(f"Error updating OME-XML metadata: {e}")
            ome_xml = None
    print(ome_xml) '''
    # Save the cropped image as a new TIFF file
    cropped_path = os.path.join(output_dir, "image_cropped.ome.tif")
    os.makedirs(os.path.dirname(cropped_path), exist_ok=True)
    tifffile.imwrite(cropped_path, cropped_img, description=ome_xml)
    print(f"Cropped image (ome.tif) saved to {cropped_path}")

    # Create a GeoDataFrame with the cropping box and save it as a GeoJSON file
    geo_path = os.path.join(output_dir, "image_cropped.geojson")
    corners_gdf = gpd.GeoDataFrame(geometry=[box(xmin, ymin, xmax, ymax)], crs="EPSG:4326")
    corners_gdf.to_file(geo_path, driver="GeoJSON")
    print(f"Cropped image (geojson) saved to {geo_path}")

def parse_crop_coordinates(csv_file):
    """
    Parses cropping coordinates from a CSV file.

    Parameters:
    - csv_file: str, path to the CSV file containing cropping coordinates.

    Returns:
    - tuple: (xmin, ymin, xmax, ymax) coordinates for cropping.
    """
    print(f"Parsing crop coordinates from {csv_file}")
    df = pd.read_csv(csv_file)
    print(df.head)
    if not {'xmin', 'ymin', 'xmax', 'ymax'}.issubset(df.columns):
        raise ValueError("CSV file must contain 'xmin', 'ymin', 'xmax', and 'ymax' columns.")
    return int(df.iloc[0]['xmin']), int(df.iloc[0]['ymin']), int(df.iloc[0]['xmax']), int(df.iloc[0]['ymax'])

def main():
    print("Starting image cropping")
    parser = argparse.ArgumentParser(
        description='Crop a region from a multi-channel TIFF image, and save the cropped image as a TIF file along with a GeoJSON file containing the cropping box.'
    )
    parser.add_argument('image_path', type=str,
        help='Path to the Comet .ome.tiff image file.')
    parser.add_argument('csv_file', type=str,
        help='Path to the CSV file containing the cropping coordinates.')
    parser.add_argument('output_dir', type=str,
        help='Directory to save the cropped image and GeoJSON file; defaults to the Comet sample directory.')

    args = parser.parse_args()

    crop_coords = parse_crop_coordinates(args.csv_file)
    print(f"Crop coordinates: {crop_coords}")
    print("Cropping image")
    crop_image(args.image_path, args.output_dir, crop_coords)


if __name__ == '__main__':
    main()