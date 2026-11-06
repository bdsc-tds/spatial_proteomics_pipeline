# VirTues cell types of the VirTues mask: the majority class of every cell's 1 µm pixels in the segmentation
# head's cell type logits (None / Unknown suppressed, as VirTues' segmentation demo)

# variables defined in other .smk files:
# OUT:                defined in Snakefile
# SEG_METHODS:        defined in Snakefile
# SAMPLE_IDS:         defined in Snakefile

if "virtues" not in SEG_METHODS:
    raise WorkflowError("VirTues cell types read the logits of the VirTues segmentation: "
                        "add virtues to segmentation: methods")

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
        --out_csv {output.cell_types} \
        &> {log}
        """

# VirTues cell types of the VirTues mask of every sample it segments
rule virtues_all:
    input:
        expand(OUT + "/cell_type_annotation/virtues/virtues/{sample_id}/cell_types.csv", sample_id=SAMPLE_IDS),
