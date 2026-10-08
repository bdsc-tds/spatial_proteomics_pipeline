# Mesmer (deepcell, via SPACEc) masks from the nuclear + membrane signals.
# Quantified in quantify_cells.smk

mesmer_dir = OUT + "/segmentation/mesmer/{sample_id}/"

rule segmentation_mesmer:
    input:
        script=f"{workflow.basedir}/scripts/_segmentation/mesmer.py",
        # from prepare_image
        image=OUT + "/image_preparation/{sample_id}/segmentation_input.ome.tif",
    output:
        cells=temp(mesmer_dir + "cell_boundaries_mask.npy"),
        nuclei=temp(mesmer_dir + "nuclei_boundaries_mask.npy"),
        cytoplasm=temp(mesmer_dir + "cytoplasm_boundaries_mask.npy"),
    params:
        model_dir=config["segmentation"]["mesmer"]["model_dir"],
    #threads: gpu_threads(gpu=1, cpu=8)
    threads: 1
    resources:
        # both input planes and the masks in memory
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
        --model_dir {params.model_dir} \
        --out_cells {output.cells} \
        --out_nuclei {output.nuclei} \
        --out_cytoplasm {output.cytoplasm} \
        &> {log}
        """
