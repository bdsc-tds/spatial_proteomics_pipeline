# KRONOS2 cell embeddings of every mask: the CLS token of a box around each cell, at native resolution
# (kronos2_embeddings.py). The mesmer rule is the base, other masks derive from it with `use rule` (all inputs
# restated). Samples imaged on Xenium also embed 10x's own cells (xenium_cells.smk)

# variables defined in other .smk files:
# os:                     imported in Snakefile
# CONFIG_PATH:            defined in Snakefile
# SAMPLES:                defined in Snakefile
# SAMPLE_IDS:             defined in Snakefile
# source_image:           defined in Snakefile
# xenium_qc_mask:         defined in Snakefile
# OUT:                    defined in Snakefile
# GPU_PARTITION:          defined in Snakefile
# GPU_EXTRA:              defined in Snakefile
# XENIUM_CT_SEG:          defined in Snakefile
# SEG_METHODS:            defined in Snakefile
# XENIUM_MASK_SAMPLE_IDS: defined in _dimensionality_reduction/xenium_cells.smk

KRONOS2 = config["dimensionality_reduction"]["kronos2"]


def kronos2_markers(panel):
    """Marker table of a panel (channel_name, kronos2_marker): kronos2/{panel}.csv next to the dataset config."""
    return os.path.join(os.path.dirname(CONFIG_PATH), "kronos2", panel + ".csv")


# every panel needs a marker table: KRONOS2 identifies channels by its own marker names
_missing_kronos2_markers = [kronos2_markers(p) for p in sorted({SAMPLES[s]["panel"] for s in SAMPLE_IDS})
                            if not os.path.exists(kronos2_markers(p))]
if _missing_kronos2_markers:
    raise WorkflowError(f"KRONOS2 needs a marker table per panel, missing: {_missing_kronos2_markers}")

# embedding of every cell of one sample's mesmer mask (cell_id, 768-d embedding)
rule kronos2_embeddings_mesmer:
    input:
        script=f"{workflow.basedir}/scripts/_dimensionality_reduction/kronos2_embeddings.py",
        markers=lambda wc: kronos2_markers(SAMPLES[wc.sample_id]["panel"]),
        # raw or background-subtracted image (source_image, Snakefile)
        image=source_image,
        # from segmentation_mesmer
        mask=OUT + "/segmentation/mesmer/{sample_id}/cell_boundaries_mask.tif",
        # Xenium only: saturation QC masks of the Xenium run (xenium_qc_mask)
        qc_mask=xenium_qc_mask,
    output:
        embeddings=OUT + "/dimensionality_reduction/kronos2/mesmer/{sample_id}/cell_embeddings.h5",
    params:
        model_dir=KRONOS2["model_dir"],
        box_size=KRONOS2["box_size"],
        cell_ids="",  # cell_id = mask label
        # Xenium only: QC-flagged pixels are set to 0, like pixels outside the cell
        qc_mask=lambda wc: f"--qc_mask {xenium_qc_mask(wc)}" if xenium_qc_mask(wc) else "",
    threads: 1
    resources:
        runtime=240,
        mem='64G',
        partition=GPU_PARTITION,
        slurm_extra=GPU_EXTRA,
    log:
        OUT + "/dimensionality_reduction/kronos2/mesmer/{sample_id}/logs/kronos2_embeddings.log"
    container:
        config["containers"]["pixi"]
    shell:
        """
        {PX} -e gpu python -u {input.script} \
        --image {input.image} \
        --markers {input.markers} \
        --mask {input.mask} \
        --cell_ids "{params.cell_ids}" \
        --model_dir {params.model_dir} \
        --box_size {params.box_size} \
        {params.qc_mask} \
        --out_h5 {output.embeddings} \
        &> {log}
        """

use rule kronos2_embeddings_mesmer as kronos2_embeddings_cellpose with:
    input:
        script=f"{workflow.basedir}/scripts/_dimensionality_reduction/kronos2_embeddings.py",
        markers=lambda wc: kronos2_markers(SAMPLES[wc.sample_id]["panel"]),
        # raw or background-subtracted image (source_image, Snakefile)
        image=source_image,
        # from segmentation_cellpose
        mask=OUT + "/segmentation/cellpose/{sample_id}/cell_boundaries_mask.tif",
        # Xenium only: saturation QC masks of the Xenium run (xenium_qc_mask)
        qc_mask=xenium_qc_mask,
    output:
        embeddings=OUT + "/dimensionality_reduction/kronos2/cellpose/{sample_id}/cell_embeddings.h5",
    log:
        OUT + "/dimensionality_reduction/kronos2/cellpose/{sample_id}/logs/kronos2_embeddings.log"

# VirTues mask disabled for now
# use rule kronos2_embeddings_mesmer as kronos2_embeddings_virtues with:
#     input:
#         script=f"{workflow.basedir}/scripts/_dimensionality_reduction/kronos2_embeddings.py",
#         markers=lambda wc: kronos2_markers(SAMPLES[wc.sample_id]["panel"]),
#         # raw or background-subtracted image (source_image, Snakefile)
#         image=source_image,
#         # from segmentation_virtues
#         mask=OUT + "/segmentation/virtues/{sample_id}/cell_boundaries_mask.tif",
#         # Xenium only: saturation QC masks of the Xenium run (xenium_qc_mask)
#         qc_mask=xenium_qc_mask,
#     output:
#         embeddings=OUT + "/dimensionality_reduction/kronos2/virtues/{sample_id}/cell_embeddings.h5",
#     log:
#         OUT + "/dimensionality_reduction/kronos2/virtues/{sample_id}/logs/kronos2_embeddings.log"

use rule kronos2_embeddings_mesmer as kronos2_embeddings_xenium with:
    input:
        script=f"{workflow.basedir}/scripts/_dimensionality_reduction/kronos2_embeddings.py",
        markers=lambda wc: kronos2_markers(SAMPLES[wc.sample_id]["panel"]),
        # raw or background-subtracted image (source_image, Snakefile)
        image=source_image,
        # from xenium_cell_mask: 10x cells of the XENIUM_CT_SEG segmentation and their cell_ids
        mask=OUT + "/dimensionality_reduction/xenium_" + XENIUM_CT_SEG + "/{sample_id}/cell_boundaries_mask.tif",
        cell_ids=OUT + "/dimensionality_reduction/xenium_" + XENIUM_CT_SEG + "/{sample_id}/cell_ids.csv",
        # Xenium only: saturation QC masks of the Xenium run (xenium_qc_mask)
        qc_mask=xenium_qc_mask,
    output:
        embeddings=OUT + "/dimensionality_reduction/kronos2/xenium_" + XENIUM_CT_SEG + "/{sample_id}/cell_embeddings.h5",
    params:
        model_dir=KRONOS2["model_dir"],
        box_size=KRONOS2["box_size"],
        cell_ids=lambda wc, input: input.cell_ids,  # cell_id = 10x cell_id
        qc_mask=lambda wc: f"--qc_mask {xenium_qc_mask(wc)}" if xenium_qc_mask(wc) else "",
    log:
        OUT + "/dimensionality_reduction/kronos2/xenium_" + XENIUM_CT_SEG + "/{sample_id}/logs/kronos2_embeddings.log"

# KRONOS2 embeddings of every method in segmentation: methods (VirTues disabled for now), and of the Xenium cells
# (samples imaged on Xenium)
rule kronos2_embeddings_all:
    input:
        expand(OUT + "/dimensionality_reduction/kronos2/{seg_method}/{sample_id}/cell_embeddings.h5",
               seg_method=[m for m in SEG_METHODS if m != "virtues"], sample_id=SAMPLE_IDS),
        expand(OUT + "/dimensionality_reduction/kronos2/xenium_" + XENIUM_CT_SEG + "/{sample_id}/cell_embeddings.h5",
               sample_id=XENIUM_MASK_SAMPLE_IDS),
