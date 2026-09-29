# rules to segment cells with Mesmer (deepcell, through SPACEc) on the nuclear + membrane channels,
# measure every channel per cell and export the masks as GeoJSON.
# quantify_cells and mask_to_geojson work on {segmentation_dir} with a {seg_method} wildcard, so a
# new segmentation method only needs its own mask rule writing *_boundaries_mask.npy there

segmentation_dir = OUT + "/segmentation/{seg_method}/{sample_id}/"
mesmer_dir = OUT + "/segmentation/mesmer/{sample_id}/"


def gpu_threads(gpu=1, cpu=8):
    """Threads for a GPU rule: few when the GPU does the work, more when it falls back to CPU."""
    return gpu if GPU else cpu


rule segmentation_mesmer:
    input:
        script=f"{workflow.basedir}/scripts/_segmentation/mesmer.py",
        image=rules.prepare_image.output.segmentation_input,
    output:
        cells=mesmer_dir + "cell_boundaries_mask.npy",
        nuclei=mesmer_dir + "nuclei_boundaries_mask.npy",
        cytoplasm=mesmer_dir + "cytoplasm_boundaries_mask.npy",
    params:
        nuclear_channel=lambda wc: shlex.quote(panel_meta(wc)["nuclear_channel"]),
        membrane_channels=lambda wc: shell_list(panel_meta(wc)["membrane_channels"]),
        model_dir=config["segmentation"]["mesmer_model_dir"],
    threads: gpu_threads(gpu=1, cpu=8)
    resources:
        # Mesmer holds the nuclear and combined membrane planes plus its outputs in memory
        runtime=360,
        mem='64G',
        partition=GPU_PARTITION,
        slurm_extra=GPU_EXTRA,
    log:
        mesmer_dir + "logs/segmentation_mesmer.log"
    container:
        config["containers"]["pixi"]
    shell:
        """
        {PX} -e spacec python {input.script} \
        --image {input.image} \
        --nuclear_channel {params.nuclear_channel} \
        --membrane_channels {params.membrane_channels} \
        --model_dir {params.model_dir} \
        --out_cells {output.cells} \
        --out_nuclei {output.nuclei} \
        --out_cytoplasm {output.cytoplasm} \
        &> {log}
        """

# area, centroid and mean intensity of every channel per cell, from the cell mask and the full
# (raw or background-subtracted) image
rule quantify_cells:
    input:
        script=f"{workflow.basedir}/scripts/_segmentation/quantify_cells.py",
        image=source_image,
        mask=segmentation_dir + "cell_boundaries_mask.npy",
    output:
        feature_table=segmentation_dir + "feature_table.csv",
    params:
        pixel_size=lambda wc: panel_meta(wc)["pixel_size_um"],
    threads: 4
    resources:
        # one channel plane plus one row block of the memory-mapped mask at a time
        runtime=480,
        mem='32G',
    log:
        segmentation_dir + "logs/quantify_cells.log"
    container:
        config["containers"]["pixi"]
    shell:
        """
        {PX} python {input.script} \
        --image {input.image} \
        --mask {input.mask} \
        --pixel_size_fallback {params.pixel_size} \
        --out_csv {output.feature_table} \
        --threads {threads} \
        &> {log}
        """

# cell, nuclei and cytoplasm masks as GeoJSON polygons (e.g. for QuPath)
rule mask_to_geojson:
    input:
        script=f"{workflow.basedir}/scripts/_segmentation/mask_to_geojson.py",
        cells=segmentation_dir + "cell_boundaries_mask.npy",
        nuclei=segmentation_dir + "nuclei_boundaries_mask.npy",
        cytoplasm=segmentation_dir + "cytoplasm_boundaries_mask.npy",
    output:
        cells=segmentation_dir + "cell_boundaries.geojson",
        nuclei=segmentation_dir + "nuclei_boundaries.geojson",
        cytoplasm=segmentation_dir + "cytoplasm_boundaries.geojson",
    resources:
        runtime=480,
        mem='64G',
    log:
        segmentation_dir + "logs/mask_to_geojson.log"
    container:
        config["containers"]["pixi"]
    shell:
        """
        ( {PX} python {input.script} {input.cells} {output.cells} && \
        {PX} python {input.script} {input.nuclei} {output.nuclei} && \
        {PX} python {input.script} {input.cytoplasm} {output.cytoplasm} ) \
        &> {log}
        """

rule segmentation_all:
    input:
        expand(rules.quantify_cells.output.feature_table, seg_method=SEG_METHODS, sample_id=SAMPLE_IDS),
        expand(rules.mask_to_geojson.output, seg_method=SEG_METHODS, sample_id=SAMPLE_IDS),
