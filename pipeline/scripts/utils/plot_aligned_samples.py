import argparse
import matplotlib.pyplot as plt
from plotnine import *
import os
import scanpy as sc

from utils import *

def plot_aligned_samples(main_dir, xenium_file, data_type, hid, rid, segmentation_method, align_method):
    raw_folder = os.path.join(main_dir, "raw", data_type, "cropped", hid, rid)
    processed_folder = os.path.join(main_dir, "processed", data_type, hid, rid)
    matrix = os.path.join(processed_folder, 'image_alignment', align_method, 'affine_matrix.csv')
    if os.path.exists(os.path.join(os.path.dirname(matrix), "error.txt")):
        print("Affine matrix does not exists")
        exit(0)
    
    print(f"Processing {data_type} with method {segmentation_method}...")
    csv_file = os.path.join(processed_folder, "segmentation", segmentation_method, "feature_table.csv")
    adata_comet = spatial_to_adata(csv_file, hid, rid, segmentation_method)
    xmin, ymin, xmax, ymax = get_coords(os.path.join(raw_folder, "image_cropped.geojson"))
    print(f"Cropping adata_comet to coordinates: xmin={xmin}, ymin={ymin}, xmax={xmax}, ymax={ymax}...")

    if segmentation_method == "horizon":
        # the cropping coordinates are in microns, so we need to convert them to pixels by multiplying by 0.28 (the mpp of the comet images)
        xmin *= 0.28
        xmax *= 0.28
        ymin *= 0.28
        ymax *= 0.28
    
    adata_subs = adata_comet[((adata_comet.obs["y_centroid"] <= ymax) &
                                (adata_comet.obs["y_centroid"] >= ymin) &
                                (adata_comet.obs["x_centroid"] <= xmax) &
                                (adata_comet.obs["x_centroid"] >= xmin)),:].copy()
    
    if segmentation_method == "mesmer":
        adata_subs.obs["x_centroid"] *= 0.28
        adata_subs.obs["y_centroid"] *= 0.28
        adata_subs.obs["area"] *= 0.28 ** 2
    
    print(adata_subs)
    print(f"Applying affine transformation for {align_method}...")
    print(f"Using affine matrix from: {matrix}")
    
    mapped_adata = transf_centroids(
        adata_subs,
        cell_id_col="cell_id",
        aff_mtx=matrix,
        crop_gpd=os.path.join(raw_folder, "image_cropped.geojson"),
        mpp_input=0.28,
        mpp_aligned=0.2125,
        centroids_microns=True,
        x_col="x_centroid",
        y_col="y_centroid",
        add_to_obs=True,
    )

    xenium_adata = sc.read_h5ad(xenium_file)
    print(xenium_adata)
    
    plt.figure(figsize=(8, 8))
    plt.scatter(xenium_adata.obs["centroid_x"], xenium_adata.obs["centroid_y"], s=0.5, alpha=0.4, label=f"xenium sample")
    plt.scatter(mapped_adata.obs["x_centroid_aligned"], mapped_adata.obs["y_centroid_aligned"], s=0.5, alpha=0.4, label=f"{data_type} - {segmentation_method} - {align_method}")
    plt.title(f"Aligned Centroids for {hid} : {data_type} with method {segmentation_method} and alignment {align_method}")
    plt.xlabel("X Centroid Aligned")
    plt.ylabel("Y Centroid Aligned")
    plt.legend()
    plt.gca().set_aspect("equal")
    
    out_path = f"{processed_folder}/image_alignment/{align_method}/plots"
    os.makedirs(out_path, exist_ok=True)
    plt.savefig(f"{out_path}/{segmentation_method}_aligned_samples.png")
    
    mapped_adata, xenium_adata = map_metadata(mapped_adata, xenium_adata, ct_cols=[], max_dist=5)
    mapped = mapped_adata.obs["matched"].sum()
    unmapped = mapped_adata.n_obs-mapped
    # put in a dataframe ncells for aligment
    df_align = pd.DataFrame({"ncells": [xenium_adata.n_obs, mapped_adata.n_obs,
                                        mapped, unmapped],
                            "modality": ["Xenium", "comet", "mapped", "unmapped"]})
    df_align["modality"] = pd.Categorical(df_align["modality"], categories=["Xenium", "comet", "mapped", "unmapped"], ordered=True)
    print(f"{hid}|{rid}|{data_type}|{segmentation_method}|{align_method}")
    print(df_align)

    cols_mapped = {
        'Xenium': '#19618A', 
        'comet': '#09B39C', # (Your existing teal/green)
        'mapped': '#FFB300',     # Bright Amber/Yellow
        'unmapped': '#D32F2F'   # Strong, Clear Red
    }
    
    # make barplot
    bar = (ggplot(df_align, aes(x="modality", y="ncells", fill="modality")) +
        geom_bar(stat="identity") +
        geom_text(aes(label="ncells"), va="bottom", nudge_y=50, size=10, format_string="{:,}") +
        theme_classic() +
        theme(panel_border=element_rect(color="black", fill=None),
                legend_text=element_text(size=15),
                legend_title=element_text(size=16),
                axis_text_x=element_text(size=14, color="black", angle=45, hjust=1, vjust=1),
                axis_text_y=element_text(size=12, color="black"),
                axis_title_y = element_text(size = 15),
                figure_size = (6,5)) +
            ylab("Number of cells") +
            xlab("") +
            scale_y_continuous(expand = [0,0,0.08,0]) +
            scale_fill_manual(values=cols_mapped)
            )
    
    bar.save(f"{out_path}/{segmentation_method}_bar_plot.png", dpi=300, height=5, width=6)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("main_dir", help="Main directory")
    parser.add_argument("xenium_file", help="Xenium file")
    parser.add_argument("data_type", help="Data type")
    parser.add_argument("hid", help="HID")
    parser.add_argument("rid", help="RID")
    parser.add_argument("segment", help="Segmentation method")
    parser.add_argument("align", help="Alignment method")
    args = parser.parse_args()

    plot_aligned_samples(args.main_dir, args.xenium_file, args.data_type, args.hid, args.rid, args.segment, args.align)

if __name__ == "__main__":
    main()