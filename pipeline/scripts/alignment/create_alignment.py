import argparse
import glob
from insitupy.io import read_xenium
from insitupy.tools import register_images
import numpy as np
import os
import tifffile
import shutil
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

def remove_excess_files(root_dir):
    '''
    Removes the excess files created during the alignment process, and moves
    the alignment matrix to the root directory.
        
    Parameters:
    - root_dir: str, path to the root directory containing the alignment results.
    '''
    print("Removing excess files.")
    normalized_dir = os.path.join(root_dir, 'normalised_results')

    # Find and move the alignment matrix CSV to the root directory
    csv_files = glob.glob(
        os.path.join(normalized_dir, '**', '*_T.csv'),
        recursive=True
    )
    if len(csv_files) == 0:
        print("Warning: No alignment matrix CSV found in normalised_results.")
    else:
        dest = os.path.join(root_dir, 'affine_matrix.csv')
        shutil.move(csv_files[0], dest)
        print(f"Moved alignment matrix: {os.path.basename(csv_files[0])} → {root_dir}")
        print(f"Alignment matrix saved to affine_matrix.csv")

    # Remove all subdirectories inside normalised_results

    shutil.rmtree(normalized_dir)
    print(f"Removed folder: {normalized_dir}")

    print("Cleanup complete.")

def align_matrix(xd_path, image_path, output_dir, min_good_matches_per_area=5):
    """
    Creates the aligned matrix between Xenium and Comet and saves it in the
    save directory as the comet sample.
    If the samples aligned are not from the same slice, try to decrease the
    'min_good_matches_per_area' parameter, to allow the alignment to be done
    with fewer matches. If 1 is not enough, the samples are likely too
    different to be aligned, and the alignment should be done manually.

    Parameters:
    - xd_path: str, path to read the Xenium data, should be a folder
    - image_path: str, path to import the Comet data, should be the exact file
    - min_good_matches_per_area: int, minimum matches needed per 5mm² areas, to complete the alignment
    """

    xd = read_xenium(xd_path)

    channel_names = get_channels(image_path)
    print(f"Channel names extracted from OME-XML: {channel_names}")
        # First channel whatever its name
    first_channel = channel_names[0]

    # First channel that is different from the first one
    second_channel = next(c for c in channel_names if c != first_channel)

    selected_channels = [first_channel, second_channel]
    print(selected_channels)
    # create folder for aligned images
    output_dir = os.path.join(output_dir, 'normalised_results')
    os.makedirs(output_dir, exist_ok=True)
    register_images(data=xd, 
        image_to_be_registered=image_path,
        axes_image="CYX",
        channel_names=selected_channels,
        channel_name_for_registration=channel_names[0],
        template_image_name="nuclei",
        axes_template="YX",
        output_dir=output_dir,
        save_registered_images=True,
        min_good_matches_per_area=min_good_matches_per_area
        )
    print("Image registration done!")
    print(f"Output directory: {output_dir}")

def main():
    parser = argparse.ArgumentParser(description="Align Comet on Xenium images.")
    parser.add_argument("xd_path", type=str, help="Path to the Xenium data folder.")
    parser.add_argument("image_path", type=str, help="Path to the Comet image file; has to be a .tif file.")
    parser.add_argument("output_dir", type=str, help="Directory to save the aligned matrix.")
    parser.add_argument("--min_good_matches_per_area", type=int, default=5, help="Minimum matches needed per 5mm² areas for alignment.")
    
    args = parser.parse_args()
    
    align_matrix(args.xd_path, args.image_path, args.output_dir, args.min_good_matches_per_area)
    remove_excess_files(args.output_dir)


if __name__ == "__main__":
    main()