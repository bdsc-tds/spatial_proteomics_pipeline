# not named backsub.py: the script's folder comes first on sys.path, so that name would shadow the
# backsub package imported below
import os
import argparse
import subprocess
import pandas as pd
from backsub import metadata2markers
from pathlib import Path

def generate_markers_csv(
    input_img: str,
    output_dir: str,
    output_file_name: str = "markers.csv",
    registration_filter: str = "DAPI",
    remove_background_channels: bool = False,
    remove_dapi: bool = False,
    keep_registration_markername: str = "DAPI",
    channel_names: str = None,
):
    """
    Generate a markers CSV file from an OME-TIFF image.
 
    Args:
        input_img:                    Path to the input image stack (.tif).
        output_dir:                   Directory where the CSV will be written.
        output_file_name:             Name of the output CSV file (default: markers.csv).
        registration_filter:          Filter used for image registration (default: DAPI).
                                      Must match a FluorescenceChannel value in the image metadata.
        remove_background_channels:   If True, marks background channels with remove=TRUE.
        remove_dapi:                  If True, marks extra DAPI reference channels with remove=TRUE.
        keep_registration_markername: Marker name to keep when remove_dapi=True (default: DAPI).
        channel_names:                CSV (raw_name, channel_name) renaming the kept channels; raw_name
                                      is the marker name without its " - <filter>-autof" suffix.
                                      None keeps the raw names.
    """
    img_path = Path(input_img)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
 
    df = metadata2markers.meta_from_file(img_path, registration_filter=registration_filter)
    df = metadata2markers.assign_background(
        df,
        rmv_ref=remove_dapi,
        rmv_back=remove_background_channels,
        ref_marker=keep_registration_markername,
        registration_filter=registration_filter,
    )

    # if marker_name contains TRICT-autof set background to TRITC-autof and if marker_name contains Cy5-autof set background to Cy5-autof
    df.loc[df["marker_name"].str.contains("TRITC-autof"), "background"] = "TRITC-autof"
    df.loc[df["marker_name"].str.contains("Cy5-autof"), "background"] = "Cy5-autof"

    # set remove to true on channels where names = "TRITC-autof" and "Cy5-autof" and background to the None
    df.loc[df["marker_name"].isin(["TRITC-autof", "Cy5-autof"]), "remove"] = "TRUE"
    df.loc[df["marker_name"].isin(["TRITC-autof", "Cy5-autof"]), "background"] = None

    # remove markers with names containing "elu" and set remove to true
    df.loc[df["marker_name"].str.contains("elu"), "remove"] = "TRUE"

    # rename the kept channels: strip the " - TRITC-autof" / " - Cy5-autof" suffix of instrument
    # AF-subtracted exports, then map to the clean name; backsub writes these names into the output.
    # The background channels keep their names, backsub finds them by exact name
    if channel_names is not None:
        mapping = pd.read_csv(channel_names).set_index("raw_name")["channel_name"]
        base = df["marker_name"].str.replace(r" - (TRITC|Cy5)-autof$", "", regex=True)
        kept = df["remove"] != "TRUE"
        unmapped = sorted(set(base[kept]) - set(mapping.index))
        if unmapped:
            raise ValueError(f"channels missing from {channel_names}: {unmapped}")
        df.loc[kept, "marker_name"] = base[kept].map(mapping)

    output_path = out_dir / output_file_name
    df.to_csv(output_path, index=False)
    print(f"Markers CSV saved to: {output_path}")
    return df


def main(img_path: str, output_dir: str, channel_names: str = None):
    df = generate_markers_csv(
        input_img=img_path,
        output_dir=output_dir,
        registration_filter="DAPI",
        remove_background_channels=False,
        remove_dapi=True,
        keep_registration_markername="DAPI",
        channel_names=channel_names,
    )
    # run the backsub CLI; check=True so that a failed subtraction fails the rule
    subprocess.run(
        ["backsub",
         "-r", img_path,
         "-o", os.path.join(output_dir, "image.ome.tiff"),
         "-m", os.path.join(output_dir, "markers.csv"),
         "-mo", os.path.join(output_dir, "markerout.csv")],
        check=True,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run background subtraction pipeline.")
    parser.add_argument("--img_path", required=True, help="Path to the input OME-TIFF image file.")
    parser.add_argument("--output_dir", required=True, help="Directory to save the output files")
    parser.add_argument("--channel_names", help="CSV (raw_name, channel_name) renaming the kept channels")
    args = parser.parse_args()

    main(img_path=args.img_path, output_dir=args.output_dir, channel_names=args.channel_names)