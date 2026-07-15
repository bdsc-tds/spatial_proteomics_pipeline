configfile: "config.yaml"

import glob
from itertools import product
import numpy as np
import os
import pandas as pd
from pathlib import Path
from pprint import pprint

main_dir = config['main_dir']
ALIGN = config['align_methods']
SEGMENT = config['segmentation_methods']

SAMPLES = glob_wildcards(main_dir + "lunaphore_protein/raw/{data_type}/uncropped/{hid}/{rid}/{sample_name}.ome.tiff")

filtered_tuples = [
    (dt, hid, rid) 
    for dt, hid, rid in zip(SAMPLES.data_type, SAMPLES.hid, SAMPLES.rid) 
    if dt == 'postXenium'
]

# 3. Reconstruct your SAMPLES object (or a custom object) for your expand functions
class SampleConfig:
    def __init__(self, tuples):
        self.data_type = [t[0] for t in tuples]
        self.hid = [t[1] for t in tuples]
        self.rid = [t[2] for t in tuples]

SAMPLES = SampleConfig(filtered_tuples)

def get_path(key):
    return os.path.join(config["main_dir"],config["paths"][key])

# function to compute the memory RAM dynamically based on whether the job fails
def get_mem_mb(wildcards, attempt):
    return attempt * 30000

def get_mem_mb_exp(wildcards, attempt):
    return attempt**2 * 25000

rule all:
    input:
        # background substraction
        expand(get_path("backsub_img"), zip, data_type=SAMPLES.data_type, hid=SAMPLES.hid, rid=SAMPLES.rid),
        expand(get_path("backsub_markers"), zip, data_type=SAMPLES.data_type, hid=SAMPLES.hid, rid=SAMPLES.rid),
        expand(get_path("backsub_markerout"), zip, data_type=SAMPLES.data_type, hid=SAMPLES.hid, rid=SAMPLES.rid),
        
        # cropping image
        expand(get_path("cropped_img"), zip, data_type=SAMPLES.data_type, hid=SAMPLES.hid, rid=SAMPLES.rid),
        expand(get_path("cropped_geojson"), zip, data_type=SAMPLES.data_type, hid=SAMPLES.hid, rid=SAMPLES.rid),
        
        # align image
        [
            get_path("align_matrix").format(**locals())
            for data_type, hid, rid in zip(SAMPLES.data_type, SAMPLES.hid, SAMPLES.rid)
            for align in ALIGN
        ],

        # plot aligned images
        [
            get_path("plot_aligned_samples").format(**locals())
            for data_type, hid, rid in zip(SAMPLES.data_type, SAMPLES.hid, SAMPLES.rid)
            for align in ALIGN
            for segment in SEGMENT
        ],
        [
            get_path("plot_aligned_bar").format(**locals())
            for data_type, hid, rid in zip(SAMPLES.data_type, SAMPLES.hid, SAMPLES.rid)
            for align in ALIGN
            for segment in SEGMENT
        ],

        # cell segmentation
        [
            get_path("seg_csv").format(**locals())
            for data_type, hid, rid in zip(SAMPLES.data_type, SAMPLES.hid, SAMPLES.rid)
            for segment in SEGMENT
        ],
        [
            get_path("seg_npy_cells").format(**locals())
            for data_type, hid, rid in zip(SAMPLES.data_type, SAMPLES.hid, SAMPLES.rid)
            for segment in ["mesmer"]
        ],
        [
            get_path("seg_npy_nuclei").format(**locals())
            for data_type, hid, rid in zip(SAMPLES.data_type, SAMPLES.hid, SAMPLES.rid)
            for segment in ["mesmer"]
        ],
        [
            get_path("seg_npy_cytoplasm").format(**locals())
            for data_type, hid, rid in zip(SAMPLES.data_type, SAMPLES.hid, SAMPLES.rid)
            for segment in ["mesmer"]
        ],
        [
            get_path("seg_txt").format(**locals())
            for data_type, hid, rid in zip(SAMPLES.data_type, SAMPLES.hid, SAMPLES.rid)
            for segment in ["mesmer"]
        ],

        # convert mask to geojson
        [
            get_path("conversion_geojson_cells_file").format(**locals())
            for data_type, hid, rid in zip(SAMPLES.data_type, SAMPLES.hid, SAMPLES.rid)
            for segment in ["mesmer"]
        ],
        [
            get_path("conversion_geojson_nuclei_file").format(**locals())
            for data_type, hid, rid in zip(SAMPLES.data_type, SAMPLES.hid, SAMPLES.rid)
            for segment in ["mesmer"]
        ],
        [
            get_path("conversion_geojson_cytoplasm_file").format(**locals())
            for data_type, hid, rid in zip(SAMPLES.data_type, SAMPLES.hid, SAMPLES.rid)
            for segment in ["mesmer"]
        ],

rule background_substraction:
    input:
        img = lambda wildcards: glob.glob(get_path("img").format(data_type=wildcards.data_type, hid=wildcards.hid, rid=wildcards.rid))
    output:
        backsubed_img = get_path("backsub_img"),
        markers_csv = get_path("backsub_markers"),
        markers_out_csv = get_path("backsub_markerout")
    resources:
        mem_mb = get_mem_mb_exp
    log:
        get_path("backsub_log")
    conda:
        "envs/backsub.yaml"
    shell:
        "python scripts/background_substraction/background_substraction.py "
        "--img_path '{input.img}' "
        "--output_dir $(dirname '{output.backsubed_img}') > {log} 2>&1"

rule crop_image:
    input:
        img = lambda wildcards: glob.glob(get_path("img").format(data_type=wildcards.data_type, hid=wildcards.hid, rid=wildcards.rid)),
        csv_file = lambda wildcards: glob.glob(get_path("cropped_coordinates").format(data_type=wildcards.data_type, hid=wildcards.hid, rid=wildcards.rid)),
    output:
        cropped_img = get_path("cropped_img"),
        cropped_geojson = get_path("cropped_geojson")
    resources:
        mem_mb = get_mem_mb_exp
    log:
        get_path("cropped_log")
    conda:
        "envs/insitupy.yaml"
    shell:
        "python scripts/alignment/image_cropping.py "
        "'{input.img}' "
        "'{input.csv_file}' "
        "$(dirname '{output.cropped_img}') > {log} 2>&1"

rule align_image:
    input:
        xenium_path = lambda wildcards: glob.glob(get_path("xenium_outs").format(hid=wildcards.hid)),
        cropped_img = rules.crop_image.output.cropped_img
    output:
        affine_matrix = get_path("align_matrix")
    resources:
        mem_mb = get_mem_mb_exp,
        current_attempt = lambda wildcards, attempt: attempt
    log:
        get_path("align_log")
    conda:
        "envs/insitupy.yaml"
    run:
        error_path =  os.path.join(os.path.dirname(output.affine_matrix), "error.txt")
        if not os.path.exists(error_path):
            try:
                shell("python scripts/alignment/create_alignment.py '{input.xenium_path}' '{input.cropped_img}' $(dirname '{output[0]}') > {log} 2>&1")
            except Exception as e:
                # Try to automate and take the retries number from the profile config file.
                if resources.current_attempt == 4:
                    with open(error_path, "w") as f:
                        f.write(f"Alignment failed")
                raise
        else:
            print("Alignment already failed.")

rule plot_alignment_and_data:
    # output a png of aligned plots, barplots for cell informations 
    input:
        #  xenium adata for comparison and plots
        xenium_path = lambda wildcards: glob.glob(get_path("xenium_adata").format(hid=wildcards.hid)),
        
        # segmentation csv file for centroid positions
        csv_file = get_path("seg_csv"),
        
        # matrix used to align images
        affine_matrix = rules.align_image.output.affine_matrix
    output:
        # plot aligned samples and barplots
        aligned = get_path("plot_aligned_samples"),
        bar = get_path("plot_aligned_bar")
    resources:
        # if the job fails, retry with more memory
        mem_mb = get_mem_mb_exp
    params:
        # to avoid trailing slash issues
        main_dir=lambda wildcards: config["main_dir"].rstrip("/")
    log:
        get_path("plot_aligned_log")
    conda:
        "envs/squid_lazy_expl.yaml"
    shell:
        "python scripts/utils/plot_aligned_samples.py "
        "{params.main_dir}/lunaphore_protein/ "
        "{input.xenium_path} "
        "{wildcards.data_type} "
        "{wildcards.hid} "
        "{wildcards.rid} "
        "{wildcards.segment} "
        "{wildcards.align} > {log} 2>&1"

rule segmentation_mesmer:
    input:
        backsubed_img = rules.background_substraction.output.backsubed_img,
        markers_out_csv = rules.background_substraction.output.markers_out_csv
    output:
        csv_file = get_path("seg_csv"),
        npy_cells_file = get_path("seg_npy_cells"),
        npy_nuclei_file = get_path("seg_npy_nuclei"),
        npy_cytoplasm_file = get_path("seg_npy_cytoplasm"),
        channels = get_path("seg_txt")
    wildcard_constraints:
        segment = "mesmer"
    resources:
        mem_mb = get_mem_mb_exp
    log:
        get_path("seg_log")
    conda:
        "envs/spacec.yaml"
    shell:
        "python scripts/segmentation/mesmer_segmentation.py "
        "{input.backsubed_img} "
        "{input.markers_out_csv} "
        "$(dirname '{output.csv_file}') > {log} 2>&1"

rule npy_to_geojson:
    input: 
        npy_cells_file = rules.segmentation_mesmer.output.npy_cells_file,
        npy_nuclei_file = rules.segmentation_mesmer.output.npy_nuclei_file,
        npy_cytoplasm_file = rules.segmentation_mesmer.output.npy_cytoplasm_file
    output:
        geojson_cells_file = get_path("conversion_geojson_cells_file"),
        geojson_nuclei_file = get_path("conversion_geojson_nuclei_file"),
        geojson_cytoplasm_file = get_path("conversion_geojson_cytoplasm_file")
    resources:
        mem_mb = get_mem_mb_exp
    log:
        get_path("conversion_log")
    conda:
        "envs/squid_lazy_expl.yaml"
    shell:
        "python scripts/utils/convert_mask.py "
        "{input.npy_cells_file} "
        "{output.geojson_cells_file} > {log} 2>&1 && "
        "python scripts/utils/convert_mask.py "
        "{input.npy_nuclei_file} "
        "{output.geojson_nuclei_file} >> {log} 2>&1 && "
        "python scripts/utils/convert_mask.py "
        "{input.npy_cytoplasm_file} "
        "{output.geojson_cytoplasm_file} >> {log} 2>&1"