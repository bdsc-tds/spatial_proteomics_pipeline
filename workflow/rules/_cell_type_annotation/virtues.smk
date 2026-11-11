# VirTues cell types of every mask: the majority class of every cell's 1 µm pixels in the VirTues segmentation
# head's cell type logits, with the None (background) logit suppressed (Unknown kept, unlike VirTues' segmentation demo).
# The virtues rule is the base, the mesmer and cellpose masks derive from it with `use rule` (all inputs restated).
# Samples imaged on Xenium also label 10x's own cells, from the mask exported by xenium_cell_mask
# (_dimensionality_reduction/xenium_cells.smk)

# variables defined in other .smk files:
# os:                 imported in Snakefile
# OUT:                defined in Snakefile
# SEG_METHODS:        defined in Snakefile
# SAMPLE_IDS:         defined in Snakefile
# SAMPLES:            defined in Snakefile
# xenium_ct_seg_dir:  defined in Snakefile
# XENIUM_CT_SEG:      defined in Snakefile

if "virtues" not in SEG_METHODS:
    raise WorkflowError("VirTues cell types read the logits of the VirTues segmentation: "
                        "add virtues to segmentation: methods")

# samples imaged on Xenium (technology xenium) whose run has the XENIUM_CT_SEG segmentation (xenium_ct_seg_dir: Snakefile)
VIRTUES_XENIUM_SAMPLE_IDS = [s for s in SAMPLE_IDS
                             if SAMPLES[s]["panel_meta"]["technology"] == "xenium" and xenium_ct_seg_dir(s)]

# cell type of every cell of one sample (cell_id, cell_type, confidence, pixel fraction of every class)
rule virtues_cell_types_virtues:
    input:
        script=f"{workflow.basedir}/scripts/_cell_type_annotation/virtues_cell_types.py",
        # from segmentation_virtues (logits, class names and its own cell mask)
        logits=OUT + "/segmentation/virtues/{sample_id}/class_logits.tif",
        classes=OUT + "/segmentation/virtues/{sample_id}/classes.csv",
        mask=OUT + "/segmentation/virtues/{sample_id}/cell_boundaries_mask.tif",
    output:
        cell_types=OUT + "/cell_type_annotation/virtues/virtues/{sample_id}/cell_types.csv",
    params:
        cell_ids="",  # cell_id = mask label
    threads: 1
    resources:
        runtime=30,
        mem='32G',
    log:
        OUT + "/cell_type_annotation/virtues/virtues/{sample_id}/logs/virtues_cell_types.log"
    container:
        config["containers"]["pixi"]
    shell:
        """
        {PX} python -u {input.script} \
        --logits {input.logits} \
        --classes {input.classes} \
        --mask {input.mask} \
        --cell_ids "{params.cell_ids}" \
        --out_csv {output.cell_types} \
        &> {log}
        """

use rule virtues_cell_types_virtues as virtues_cell_types_xenium with:
    input:
        script=f"{workflow.basedir}/scripts/_cell_type_annotation/virtues_cell_types.py",
        # from segmentation_virtues (logits and class names)
        logits=OUT + "/segmentation/virtues/{sample_id}/class_logits.tif",
        classes=OUT + "/segmentation/virtues/{sample_id}/classes.csv",
        # from xenium_cell_mask: 10x cells of the XENIUM_CT_SEG segmentation and their cell_ids
        mask=OUT + "/dimensionality_reduction/xenium_" + XENIUM_CT_SEG + "/{sample_id}/cell_boundaries_mask.tif",
        cell_ids=OUT + "/dimensionality_reduction/xenium_" + XENIUM_CT_SEG + "/{sample_id}/cell_ids.csv",
    output:
        cell_types=OUT + "/cell_type_annotation/virtues/xenium_" + XENIUM_CT_SEG + "/{sample_id}/cell_types.csv",
    params:
        cell_ids=lambda wc, input: input.cell_ids,  # cell_id = 10x cell_id
    log:
        OUT + "/cell_type_annotation/virtues/xenium_" + XENIUM_CT_SEG + "/{sample_id}/logs/virtues_cell_types.log"

use rule virtues_cell_types_virtues as virtues_cell_types_mesmer with:
    input:
        script=f"{workflow.basedir}/scripts/_cell_type_annotation/virtues_cell_types.py",
        # from segmentation_virtues (logits and class names)
        logits=OUT + "/segmentation/virtues/{sample_id}/class_logits.tif",
        classes=OUT + "/segmentation/virtues/{sample_id}/classes.csv",
        # from segmentation_mesmer
        mask=OUT + "/segmentation/mesmer/{sample_id}/cell_boundaries_mask.tif",
    output:
        cell_types=OUT + "/cell_type_annotation/virtues/mesmer/{sample_id}/cell_types.csv",
    log:
        OUT + "/cell_type_annotation/virtues/mesmer/{sample_id}/logs/virtues_cell_types.log"

use rule virtues_cell_types_virtues as virtues_cell_types_cellpose with:
    input:
        script=f"{workflow.basedir}/scripts/_cell_type_annotation/virtues_cell_types.py",
        # from segmentation_virtues (logits and class names)
        logits=OUT + "/segmentation/virtues/{sample_id}/class_logits.tif",
        classes=OUT + "/segmentation/virtues/{sample_id}/classes.csv",
        # from segmentation_cellpose
        mask=OUT + "/segmentation/cellpose/{sample_id}/cell_boundaries_mask.tif",
    output:
        cell_types=OUT + "/cell_type_annotation/virtues/cellpose/{sample_id}/cell_types.csv",
    log:
        OUT + "/cell_type_annotation/virtues/cellpose/{sample_id}/logs/virtues_cell_types.log"

# VirTues cell types of every method in segmentation: methods, and of the Xenium cells (samples imaged on Xenium)
rule virtues_all:
    input:
        expand(OUT + "/cell_type_annotation/virtues/{seg_method}/{sample_id}/cell_types.csv",
               seg_method=SEG_METHODS, sample_id=SAMPLE_IDS),
        expand(OUT + "/cell_type_annotation/virtues/xenium_" + XENIUM_CT_SEG + "/{sample_id}/cell_types.csv",
               sample_id=VIRTUES_XENIUM_SAMPLE_IDS),
