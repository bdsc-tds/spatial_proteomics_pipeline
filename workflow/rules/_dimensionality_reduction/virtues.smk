# VirTues cell embeddings of every mask: the overlap-weighted mean of the backbone's patch summary tokens around
# each cell (virtues_embeddings.py). The 1 µm stack, standardization and marker embeddings come from the VirTues
# segmentation rules (_segmentation/virtues.smk). The mesmer rule is the base, other masks derive from it with
# `use rule` (all inputs restated). Samples imaged on Xenium also embed 10x's own cells (xenium_cells.smk)

# variables defined in other .smk files:
# SEG_METHODS:            defined in Snakefile
# OUT:                    defined in Snakefile
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

# embedding of every cell of one sample's mesmer mask (cell_id, 512-d embedding)
rule virtues_embeddings_mesmer:
    input:
        script=f"{workflow.basedir}/scripts/_dimensionality_reduction/virtues_embeddings.py",
        # from virtues_prepare_image
        image_1um=OUT + "/segmentation/virtues/{sample_id}/image_1um.npy",
        channel_stats=OUT + "/segmentation/virtues/{sample_id}/channel_stats.csv",
        # from virtues_standardization_stats and virtues_marker_embeddings (the sample's panel)
        standardization=lambda wc: OUT + f"/segmentation/virtues/panels/{SAMPLES[wc.sample_id]['panel']}/standardization.csv",
        embeddings=lambda wc: OUT + f"/segmentation/virtues/panels/{SAMPLES[wc.sample_id]['panel']}/marker_embeddings.pt",
        # raw or background-subtracted image (source_image, Snakefile): its size only
        image=source_image,
        # from segmentation_mesmer
        mask=OUT + "/segmentation/mesmer/{sample_id}/cell_boundaries_mask.tif",
    output:
        embeddings=OUT + "/dimensionality_reduction/virtues/mesmer/{sample_id}/cell_embeddings.npz",
    params:
        model_dir=VIRTUES["model_dir"],
        model=VIRTUES["model"],
        blur_sigma=VIRTUES["blur_sigma"],
        cell_ids="",  # cell_id = mask label
    threads: 1
    resources:
        runtime=240,
        mem='64G',
        partition=VIRTUES_GPU_PARTITION,
        slurm_extra=GPU_EXTRA,
    log:
        OUT + "/dimensionality_reduction/virtues/mesmer/{sample_id}/logs/virtues_embeddings.log"
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
        --out_npz {output.embeddings} \
        &> {log}
        """

use rule virtues_embeddings_mesmer as virtues_embeddings_cellpose with:
    input:
        script=f"{workflow.basedir}/scripts/_dimensionality_reduction/virtues_embeddings.py",
        # from virtues_prepare_image
        image_1um=OUT + "/segmentation/virtues/{sample_id}/image_1um.npy",
        channel_stats=OUT + "/segmentation/virtues/{sample_id}/channel_stats.csv",
        # from virtues_standardization_stats and virtues_marker_embeddings (the sample's panel)
        standardization=lambda wc: OUT + f"/segmentation/virtues/panels/{SAMPLES[wc.sample_id]['panel']}/standardization.csv",
        embeddings=lambda wc: OUT + f"/segmentation/virtues/panels/{SAMPLES[wc.sample_id]['panel']}/marker_embeddings.pt",
        # raw or background-subtracted image (source_image, Snakefile): its size only
        image=source_image,
        # from segmentation_cellpose
        mask=OUT + "/segmentation/cellpose/{sample_id}/cell_boundaries_mask.tif",
    output:
        embeddings=OUT + "/dimensionality_reduction/virtues/cellpose/{sample_id}/cell_embeddings.npz",
    log:
        OUT + "/dimensionality_reduction/virtues/cellpose/{sample_id}/logs/virtues_embeddings.log"

use rule virtues_embeddings_mesmer as virtues_embeddings_virtues with:
    input:
        script=f"{workflow.basedir}/scripts/_dimensionality_reduction/virtues_embeddings.py",
        # from virtues_prepare_image
        image_1um=OUT + "/segmentation/virtues/{sample_id}/image_1um.npy",
        channel_stats=OUT + "/segmentation/virtues/{sample_id}/channel_stats.csv",
        # from virtues_standardization_stats and virtues_marker_embeddings (the sample's panel)
        standardization=lambda wc: OUT + f"/segmentation/virtues/panels/{SAMPLES[wc.sample_id]['panel']}/standardization.csv",
        embeddings=lambda wc: OUT + f"/segmentation/virtues/panels/{SAMPLES[wc.sample_id]['panel']}/marker_embeddings.pt",
        # raw or background-subtracted image (source_image, Snakefile): its size only
        image=source_image,
        # from segmentation_virtues
        mask=OUT + "/segmentation/virtues/{sample_id}/cell_boundaries_mask.tif",
    output:
        embeddings=OUT + "/dimensionality_reduction/virtues/virtues/{sample_id}/cell_embeddings.npz",
    log:
        OUT + "/dimensionality_reduction/virtues/virtues/{sample_id}/logs/virtues_embeddings.log"

use rule virtues_embeddings_mesmer as virtues_embeddings_xenium with:
    input:
        script=f"{workflow.basedir}/scripts/_dimensionality_reduction/virtues_embeddings.py",
        # from virtues_prepare_image
        image_1um=OUT + "/segmentation/virtues/{sample_id}/image_1um.npy",
        channel_stats=OUT + "/segmentation/virtues/{sample_id}/channel_stats.csv",
        # from virtues_standardization_stats and virtues_marker_embeddings (the sample's panel)
        standardization=lambda wc: OUT + f"/segmentation/virtues/panels/{SAMPLES[wc.sample_id]['panel']}/standardization.csv",
        embeddings=lambda wc: OUT + f"/segmentation/virtues/panels/{SAMPLES[wc.sample_id]['panel']}/marker_embeddings.pt",
        # raw or background-subtracted image (source_image, Snakefile): its size only
        image=source_image,
        # from xenium_cell_mask: 10x cells of the XENIUM_CT_SEG segmentation and their cell_ids
        mask=OUT + "/dimensionality_reduction/xenium_" + XENIUM_CT_SEG + "/{sample_id}/cell_boundaries_mask.tif",
        cell_ids=OUT + "/dimensionality_reduction/xenium_" + XENIUM_CT_SEG + "/{sample_id}/cell_ids.csv",
    output:
        embeddings=OUT + "/dimensionality_reduction/virtues/xenium_" + XENIUM_CT_SEG + "/{sample_id}/cell_embeddings.npz",
    params:
        model_dir=VIRTUES["model_dir"],
        model=VIRTUES["model"],
        blur_sigma=VIRTUES["blur_sigma"],
        cell_ids=lambda wc, input: input.cell_ids,  # cell_id = 10x cell_id
    log:
        OUT + "/dimensionality_reduction/virtues/xenium_" + XENIUM_CT_SEG + "/{sample_id}/logs/virtues_embeddings.log"

# VirTues embeddings of every method in segmentation: methods, and of the Xenium cells (samples imaged on Xenium)
rule virtues_embeddings_all:
    input:
        expand(OUT + "/dimensionality_reduction/virtues/{seg_method}/{sample_id}/cell_embeddings.npz",
               seg_method=SEG_METHODS, sample_id=SAMPLE_IDS),
        expand(OUT + "/dimensionality_reduction/virtues/xenium_" + XENIUM_CT_SEG + "/{sample_id}/cell_embeddings.npz",
               sample_id=XENIUM_MASK_SAMPLE_IDS),
