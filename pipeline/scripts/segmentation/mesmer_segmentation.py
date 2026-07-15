import argparse
import numpy as np
import os
import pandas as pd
from pathlib import Path
import spacec as sp
import tifffile
import warnings
import xml.etree.ElementTree as ET

warnings.filterwarnings('ignore')

# Same as in utils, but without all the specific imports not used for segmentation.
# Cleaner conda environments
def get_channels(markerout, channel_path):
    '''
    Extracts the channel names from the OME-XML metadata of a
    multi-channel TIFF image.
    
    Parameters:
    - image_path: str, path to the input multi-channel TIF image.
    '''
    df = pd.read_csv(markerout)
    channel_names = df['marker_name'].tolist()
    with open(channel_path, 'w') as f:
        f.write('\n'.join(channel_names))

    print(f"✅ {len(channel_names)} channels written to {channel_path}:")
    for i, n in enumerate(channel_names):
        print(f"  {i}: {n}")
    return channel_names

def mesmer_segmentation(img, markerout, output_dir):
    print("Segmenting image: ", img)
    os.makedirs(output_dir, exist_ok=True)
    channel_path = os.path.join(output_dir, "channels.txt")
    if not os.path.exists(channel_path):
        channel_names = get_channels(markerout, channel_path)
    else:
        channel_names = []
        with open(channel_path, 'r') as f:
            channel_names = [line.strip() for line in f]
        print(channel_names)

    seg_output1 = sp.tl.cell_segmentation(
        file_name = img,
        channel_file = channel_path,
        output_dir = output_dir,
        seg_method ='mesmer',
        nuclei_channel = channel_names[0],
        output_fname = 'feature_table',
        membrane_channel_list = channel_names[1:], #default is None; if provide more than one channel, then they will be combined
        # compartment = 'whole-cell', # mesmer # segment whole cells or nuclei only
        differentiate_nucleus_cytoplasm = True,
        input_format ='Multichannel', # Phenocycler or codex
        resize_factor=1, # default is 1; if the image is too large, lower the value. Lower values will speed up the segmentation but may reduce the accuracy.
        size_cutoff = 0)

    print(seg_output1.keys())

    # rename combined csv in feature_table.csv
    src = Path(os.path.join(output_dir, "feature_table_combined_features.csv"))
    dst = Path(os.path.join(output_dir, "feature_table.csv"))

    if os.path.exists(src):
        src.rename(dst)
    else:
        print("CSV Files were not created.")

    # Save segmentation masks output as npy
    np.save(os.path.join(output_dir, 'cell_boundaries_mask.npy'), seg_output1["masks"])
    np.save(os.path.join(output_dir, 'nuclei_boundaries_mask.npy'), seg_output1["masks_nuclei"])
    np.save(os.path.join(output_dir, 'cytoplasm_boundaries_mask.npy'), seg_output1["masks_cytoplasm"])

def main():
    parser = argparse.ArgumentParser(description="Segment the input image with the mesmer algorithm")
    parser.add_argument("img", help="The .ome.tiff image to segment.")
    parser.add_argument("markerout", help="The marker_out file created after the background substraction.")
    parser.add_argument("output_dir", help="Output folder.")
    args = parser.parse_args()
    
    mesmer_segmentation(args.img, args.markerout, args.output_dir)

if __name__ == "__main__":
    main()