# Cellpose 3 (via SPACEc) masks: nuclear model on the nuclear signal, cell model on membrane + nuclear.
# Quantified in quantify_cells.smk

cellpose_dir = OUT + "/segmentation/cellpose/{sample_id}/"

rule segmentation_cellpose:
    input:
        script=f"{workflow.basedir}/scripts/_segmentation/cellpose_segmentation.py",
        # from prepare_image
        image=OUT + "/image_preparation/{sample_id}/segmentation_input.ome.tif",
    output:
        cells=temp(cellpose_dir + "cell_boundaries_mask.npy"),
        nuclei=temp(cellpose_dir + "nuclei_boundaries_mask.npy"),
        cytoplasm=temp(cellpose_dir + "cytoplasm_boundaries_mask.npy"),
    params:
        model_dir=config["segmentation"]["cellpose"]["model_dir"],
        model=config["segmentation"]["cellpose"]["model"],
        nuclear_model=config["segmentation"]["cellpose"]["nuclear_model"],
        diameter_um=config["segmentation"]["cellpose"]["diameter_um"],
        nuclear_diameter_um=config["segmentation"]["cellpose"]["nuclear_diameter_um"],
    threads: 1
    resources:
        # both input planes, SPACEc's copies and the masks in memory
        runtime=360,
        mem='80G',
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
        --out_cells {output.cells} \
        --out_nuclei {output.nuclei} \
        --out_cytoplasm {output.cytoplasm} \
        &> {log}
        """
