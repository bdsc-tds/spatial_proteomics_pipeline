# VirTues cell embeddings of every mask: the overlap-weighted mean of the backbone's patch summary tokens around
# each cell (virtues_embeddings.py). The 1 µm stack, standardization and marker embeddings come from the VirTues
# segmentation rules (_segmentation/virtues.smk). One rule on {seg_method} for every mask. Samples imaged on Xenium
# also embed 10x's own cells (xenium_cells.smk)

# variables defined in other .smk files:
# SEG_METHODS:            defined in Snakefile
# step_seg_methods:       defined in Snakefile
# seg_targets:            defined in Snakefile
# results:                defined in Snakefile
# script:                 defined in Snakefile
# SAMPLES:                defined in Snakefile
# source_image:           defined in Snakefile
# VIRTUES:                defined in _segmentation/virtues.smk
# VIRTUES_GPU_PARTITION:  defined in Snakefile
# GPU_EXTRA:              defined in Snakefile
# XENIUM_CT_SEG:          defined in Snakefile
# SAMPLE_IDS:             defined in Snakefile
# XENIUM_MASK_SAMPLE_IDS: defined in _dimensionality_reduction/xenium_cells.smk

if "virtues" not in SEG_METHODS:
    raise WorkflowError("VirTues embeddings read the inputs prepared by the VirTues segmentation rules: "
                        "add virtues to segmentation: methods")

# segmentation methods embedded (dimensionality_reduction: virtues: seg_methods, default every method)
VIRTUES_DR_SEG_METHODS = step_seg_methods(config["dimensionality_reduction"].get("virtues") or {})

# embedding of every cell of one sample's mask (cell_id, 512-d embedding)
rule virtues_embeddings:
    input:
        script=script("_dimensionality_reduction/virtues_embeddings.py"),
        # from virtues_prepare_image
        image_1um=results("segmentation/virtues/{sample_id}/image_1um.npy"),
        channel_stats=results("segmentation/virtues/{sample_id}/channel_stats.csv"),
        # from virtues_standardization_stats and virtues_marker_embeddings (the sample's panel)
        standardization=lambda wc: results(f"segmentation/virtues/panels/{SAMPLES[wc.sample_id]['panel']}/standardization.csv"),
        embeddings=lambda wc: results(f"segmentation/virtues/panels/{SAMPLES[wc.sample_id]['panel']}/marker_embeddings.pt"),
        # raw or background-subtracted image (source_image, Snakefile): its size only
        image=source_image,
        # from the method's segmentation rule
        mask=results("segmentation/{seg_method}/{sample_id}/cell_boundaries_mask.tif"),
    output:
        embeddings=results("dimensionality_reduction/virtues/{seg_method}/{sample_id}/cell_embeddings.h5"),
    params:
        model_dir=VIRTUES["model_dir"],
        model=VIRTUES["model"],
        blur_sigma=VIRTUES["blur_sigma"],
        cell_ids="",  # cell_id = mask label
    threads: 1
    resources:
        runtime=120,
        mem='96G',
        partition=VIRTUES_GPU_PARTITION,
        slurm_extra=GPU_EXTRA,
    log:
        results("dimensionality_reduction/virtues/{seg_method}/{sample_id}/logs/virtues_embeddings.log")
    container:
        config["containers"]["pixi"]
    shell:
        """
        {PX} -e gpu python -u {input.script} \
        --image_1um {input.image_1um} \
        --channel_stats {input.channel_stats} \
        --standardization {input.standardization} \
        --embeddings {input.embeddings} \
        --image {input.image} \
        --mask {input.mask} \
        --cell_ids "{params.cell_ids}" \
        --model_dir {params.model_dir} \
        --model {params.model} \
        --blur_sigma {params.blur_sigma} \
        --out_h5 {output.embeddings} \
        &> {log}
        """

use rule virtues_embeddings as virtues_embeddings_xenium with:
    input:
        script=script("_dimensionality_reduction/virtues_embeddings.py"),
        # from virtues_prepare_image
        image_1um=results("segmentation/virtues/{sample_id}/image_1um.npy"),
        channel_stats=results("segmentation/virtues/{sample_id}/channel_stats.csv"),
        # from virtues_standardization_stats and virtues_marker_embeddings (the sample's panel)
        standardization=lambda wc: results(f"segmentation/virtues/panels/{SAMPLES[wc.sample_id]['panel']}/standardization.csv"),
        embeddings=lambda wc: results(f"segmentation/virtues/panels/{SAMPLES[wc.sample_id]['panel']}/marker_embeddings.pt"),
        # raw or background-subtracted image (source_image, Snakefile): its size only
        image=source_image,
        # from xenium_cell_mask: 10x cells of the XENIUM_CT_SEG segmentation and their cell_ids
        mask=results("dimensionality_reduction/xenium_" + XENIUM_CT_SEG + "/{sample_id}/cell_boundaries_mask.tif"),
        cell_ids=results("dimensionality_reduction/xenium_" + XENIUM_CT_SEG + "/{sample_id}/cell_ids.csv"),
    output:
        embeddings=results("dimensionality_reduction/virtues/xenium_" + XENIUM_CT_SEG + "/{sample_id}/cell_embeddings.h5"),
    params:
        model_dir=VIRTUES["model_dir"],
        model=VIRTUES["model"],
        blur_sigma=VIRTUES["blur_sigma"],
        cell_ids=lambda wc, input: input.cell_ids,  # cell_id = 10x cell_id
    log:
        results("dimensionality_reduction/virtues/xenium_" + XENIUM_CT_SEG + "/{sample_id}/logs/virtues_embeddings.log")

# VirTues embeddings of every method in VIRTUES_DR_SEG_METHODS (on its own samples), and of the Xenium cells
# (samples imaged on Xenium)
rule virtues_embeddings_all:
    input:
        seg_targets(results("dimensionality_reduction/virtues/{seg_method}/{sample_id}/cell_embeddings.h5"),
                    methods=VIRTUES_DR_SEG_METHODS),
        expand(results("dimensionality_reduction/virtues/xenium_" + XENIUM_CT_SEG + "/{sample_id}/cell_embeddings.h5"),
               sample_id=XENIUM_MASK_SAMPLE_IDS),
