# VirTues cell types of every mask: the majority class of every cell's 1 µm pixels in the VirTues segmentation
# head's cell type logits, with the None (background) logit suppressed (Unknown kept, unlike VirTues' segmentation demo).
# One rule on {seg_method} for every mask (its logits always come from the VirTues segmentation).
# Samples imaged on Xenium also label 10x's own cells, from the mask exported by xenium_cell_mask
# (_dimensionality_reduction/xenium_cells.smk)

# variables defined in other .smk files:
# os:                 imported in Snakefile
# results:            defined in Snakefile
# script:             defined in Snakefile
# step_seg_methods:   defined in Snakefile
# seg_targets:        defined in Snakefile
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

# segmentation methods labelled (cell_type_annotation: virtues: seg_methods, default every method)
VIRTUES_CT_SEG_METHODS = step_seg_methods(config["cell_type_annotation"].get("virtues") or {})

# cell type of every cell of one sample (cell_id, cell_type, confidence, pixel fraction of every class)
rule virtues_cell_types:
    input:
        script=script("_cell_type_annotation/virtues_cell_types.py"),
        # from segmentation_virtues (logits and class names)
        logits=results("segmentation/virtues/{sample_id}/class_logits.tif"),
        classes=results("segmentation/virtues/{sample_id}/classes.csv"),
        # from the method's segmentation rule
        mask=results("segmentation/{seg_method}/{sample_id}/cell_boundaries_mask.tif"),
    output:
        cell_types=results("cell_type_annotation/virtues/{seg_method}/{sample_id}/cell_types.csv"),
    params:
        cell_ids="",  # cell_id = mask label
    threads: 1
    resources:
        runtime=30,
        mem='32G',
    log:
        results("cell_type_annotation/virtues/{seg_method}/{sample_id}/logs/virtues_cell_types.log")
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

use rule virtues_cell_types as virtues_cell_types_xenium with:
    input:
        script=script("_cell_type_annotation/virtues_cell_types.py"),
        # from segmentation_virtues (logits and class names)
        logits=results("segmentation/virtues/{sample_id}/class_logits.tif"),
        classes=results("segmentation/virtues/{sample_id}/classes.csv"),
        # from xenium_cell_mask: 10x cells of the XENIUM_CT_SEG segmentation and their cell_ids
        mask=results("dimensionality_reduction/xenium_" + XENIUM_CT_SEG + "/{sample_id}/cell_boundaries_mask.tif"),
        cell_ids=results("dimensionality_reduction/xenium_" + XENIUM_CT_SEG + "/{sample_id}/cell_ids.csv"),
    output:
        cell_types=results("cell_type_annotation/virtues/xenium_" + XENIUM_CT_SEG + "/{sample_id}/cell_types.csv"),
    params:
        cell_ids=lambda wc, input: input.cell_ids,  # cell_id = 10x cell_id
    log:
        results("cell_type_annotation/virtues/xenium_" + XENIUM_CT_SEG + "/{sample_id}/logs/virtues_cell_types.log")

# VirTues cell types of every method in VIRTUES_CT_SEG_METHODS (on its own samples), and of the Xenium cells
# (samples imaged on Xenium)
rule virtues_all:
    input:
        seg_targets(results("cell_type_annotation/virtues/{seg_method}/{sample_id}/cell_types.csv"),
                    methods=VIRTUES_CT_SEG_METHODS),
        expand(results("cell_type_annotation/virtues/xenium_" + XENIUM_CT_SEG + "/{sample_id}/cell_types.csv"),
               sample_id=VIRTUES_XENIUM_SAMPLE_IDS),
