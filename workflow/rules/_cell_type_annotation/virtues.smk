# VirTues cell types per segmentation method: every cell gets the majority class of its pixels in the class map of
# segmentation_virtues. The mesmer rule is the base, other masks derive from it with `use rule` (all inputs
# restated). Samples imaged on Xenium also get cell types for 10x's own cells (xenium_{xenium_segmentation} folder)

# samples imaged on Xenium (technology xenium) whose run has the XENIUM_CT_SEG segmentation (xenium_ct_seg_dir: Snakefile)
VIRTUES_XENIUM_SAMPLE_IDS = [s for s in SAMPLE_IDS
                             if SAMPLES[s]["panel_meta"]["technology"] == "xenium" and xenium_ct_seg_dir(s)]

# cell type of every cell of one sample (cell_id, cell_type, confidence, pixel fraction of every class)
rule virtues_cell_types_mesmer:
    input:
        script=f"{workflow.basedir}/scripts/_cell_type_annotation/virtues_cell_types.py",
        # from segmentation_virtues
        class_map=OUT + "/segmentation/virtues/{sample_id}/class_map.tif",
        classes=OUT + "/segmentation/virtues/{sample_id}/classes.csv",
        # from segmentation_mesmer
        mask=OUT + "/segmentation/mesmer/{sample_id}/cell_boundaries_mask.tif",
    output:
        cell_types=OUT + "/cell_type_annotation/virtues/mesmer/{sample_id}/cell_types.csv",
    threads: 1
    resources:
        runtime=60,
        mem='16G',
    log:
        OUT + "/cell_type_annotation/virtues/mesmer/{sample_id}/logs/virtues_cell_types.log"
    container:
        config["containers"]["pixi"]
    shell:
        """
        {PX} python -u {input.script} \
        --class_map {input.class_map} \
        --classes {input.classes} \
        --mask {input.mask} \
        --out_csv {output.cell_types} \
        &> {log}
        """

use rule virtues_cell_types_mesmer as virtues_cell_types_cellpose with:
    input:
        script=f"{workflow.basedir}/scripts/_cell_type_annotation/virtues_cell_types.py",
        # from segmentation_virtues
        class_map=OUT + "/segmentation/virtues/{sample_id}/class_map.tif",
        classes=OUT + "/segmentation/virtues/{sample_id}/classes.csv",
        # from segmentation_cellpose
        mask=OUT + "/segmentation/cellpose/{sample_id}/cell_boundaries_mask.tif",
    output:
        cell_types=OUT + "/cell_type_annotation/virtues/cellpose/{sample_id}/cell_types.csv",
    log:
        OUT + "/cell_type_annotation/virtues/cellpose/{sample_id}/logs/virtues_cell_types.log"

use rule virtues_cell_types_mesmer as virtues_cell_types_virtues with:
    input:
        script=f"{workflow.basedir}/scripts/_cell_type_annotation/virtues_cell_types.py",
        # from segmentation_virtues (class map and its own cell mask)
        class_map=OUT + "/segmentation/virtues/{sample_id}/class_map.tif",
        classes=OUT + "/segmentation/virtues/{sample_id}/classes.csv",
        mask=OUT + "/segmentation/virtues/{sample_id}/cell_boundaries_mask.tif",
    output:
        cell_types=OUT + "/cell_type_annotation/virtues/virtues/{sample_id}/cell_types.csv",
    log:
        OUT + "/cell_type_annotation/virtues/virtues/{sample_id}/logs/virtues_cell_types.log"

use rule virtues_cell_types_mesmer as virtues_cell_types_xenium with:
    input:
        script=f"{workflow.basedir}/scripts/_cell_type_annotation/virtues_cell_types.py",
        # from segmentation_virtues
        class_map=OUT + "/segmentation/virtues/{sample_id}/class_map.tif",
        classes=OUT + "/segmentation/virtues/{sample_id}/classes.csv",
        # 10x cell labels of the XENIUM_CT_SEG segmentation, read with read_xenium_sample
        mask=lambda wc: os.path.join(xenium_ct_seg_dir(wc.sample_id), "cells.zarr.zip"),
    output:
        cell_types=OUT + "/cell_type_annotation/virtues/xenium_" + XENIUM_CT_SEG + "/{sample_id}/cell_types.csv",
    log:
        OUT + "/cell_type_annotation/virtues/xenium_" + XENIUM_CT_SEG + "/{sample_id}/logs/virtues_cell_types.log"

# VirTues cell types of every method in segmentation: methods, and of the Xenium cells (samples imaged on Xenium)
rule virtues_all:
    input:
        expand(OUT + "/cell_type_annotation/virtues/{seg_method}/{sample_id}/cell_types.csv",
               seg_method=SEG_METHODS, sample_id=SAMPLE_IDS),
        expand(OUT + "/cell_type_annotation/virtues/xenium_" + XENIUM_CT_SEG + "/{sample_id}/cell_types.csv",
               sample_id=VIRTUES_XENIUM_SAMPLE_IDS),
