# Cellpose 3 (via SPACEc) masks: nuclear model on the nuclear signal, cell model on membrane + nuclear.
# Diameters from the config, or estimated per sample ("auto", Cellpose size models); the values used go to diameters.json.
# Quantified in quantify_cells.smk

# variables defined in other .smk files:
# OUT:           defined in Snakefile
# GPU_PARTITION: defined in Snakefile
# GPU_EXTRA:     defined in Snakefile

cellpose_dir = OUT + "/segmentation/cellpose/{sample_id}/"

rule segmentation_cellpose:
    input:
        script=f"{workflow.basedir}/scripts/_segmentation/cellpose_segmentation.py",
        # from prepare_segmentation_input
        image=OUT + "/image_preparation/{sample_id}/segmentation_input.ome.tif",
    output:
        cells=cellpose_dir + "cell_boundaries_mask.tif",
        nuclei=cellpose_dir + "nuclei_boundaries_mask.tif",
        cytoplasm=cellpose_dir + "cytoplasm_boundaries_mask.tif",
        diameters=cellpose_dir + "diameters.json",
    params:
        model_dir=config["segmentation"]["cellpose"]["model_dir"],
        model=config["segmentation"]["cellpose"]["model"],
        nuclear_model=config["segmentation"]["cellpose"]["nuclear_model"],
        diameter_um=config["segmentation"]["cellpose"]["diameter_um"],
        nuclear_diameter_um=config["segmentation"]["cellpose"]["nuclear_diameter_um"],
        # clamping range of the estimated diameters (diameter "auto")
        auto_range_um=shell_list(config["segmentation"]["cellpose"].get("auto_diameter_range_um", [5, 25])),
        auto_nuclear_range_um=shell_list(config["segmentation"]["cellpose"].get("auto_nuclear_diameter_range_um", [4, 15])),
    threads: 1
    resources:
        # both input planes, SPACEc's copies and the masks in memory; a few more minutes for the diameter estimate
        runtime=200,
        mem='100G',
        partition=GPU_PARTITION,
        slurm_extra=GPU_EXTRA,
    log:
        cellpose_dir + "logs/segmentation_cellpose.log"
    container:
        config["containers"]["pixi"]
    shell:
        """
        {PX} -e spacec python {input.script} \
        --image {input.image} \
        --model_dir {params.model_dir} \
        --model {params.model} \
        --nuclear_model {params.nuclear_model} \
        --diameter_um {params.diameter_um} \
        --nuclear_diameter_um {params.nuclear_diameter_um} \
        --auto_diameter_range_um {params.auto_range_um} \
        --auto_nuclear_diameter_range_um {params.auto_nuclear_range_um} \
        --out_cells {output.cells} \
        --out_nuclei {output.nuclei} \
        --out_cytoplasm {output.cytoplasm} \
        --out_diameters {output.diameters} \
        &> {log}
        """
