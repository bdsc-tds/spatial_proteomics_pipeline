# 10x's cells of samples imaged on Xenium as a label mask in the segmentation methods' format, with their 10x
# cell_ids, for the embedding methods (their GPU environment cannot read SpatialData)

# variables defined in other .smk files:
# os:                imported in Snakefile
# DR_SAMPLE_IDS:     defined in Snakefile
# SAMPLES:           defined in Snakefile
# xenium_ct_seg_dir: defined in Snakefile
# OUT:               defined in Snakefile
# XENIUM_CT_SEG:     defined in Snakefile

# samples imaged on Xenium (technology xenium) whose run has the XENIUM_CT_SEG segmentation (xenium_ct_seg_dir: Snakefile)
XENIUM_MASK_SAMPLE_IDS = [s for s in DR_SAMPLE_IDS
                          if SAMPLES[s]["panel_meta"]["technology"] == "xenium" and xenium_ct_seg_dir(s)]

# cell label mask (tiled TIFF) and label -> 10x cell_id table of one sample
rule xenium_cell_mask:
    input:
        script=f"{workflow.basedir}/scripts/_dimensionality_reduction/xenium_cell_mask.py",
        # 10x cell labels and cell table of the XENIUM_CT_SEG segmentation, read with read_xenium_sample
        cells=lambda wc: os.path.join(xenium_ct_seg_dir(wc.sample_id), "cells.zarr.zip"),
    output:
        mask=OUT + "/dimensionality_reduction/xenium_" + XENIUM_CT_SEG + "/{sample_id}/cell_boundaries_mask.tif",
        cell_ids=OUT + "/dimensionality_reduction/xenium_" + XENIUM_CT_SEG + "/{sample_id}/cell_ids.csv",
    params:
        seg_dir=lambda wc: xenium_ct_seg_dir(wc.sample_id),
    threads: 1
    resources:
        runtime=60,
        mem='48G',
    log:
        OUT + "/dimensionality_reduction/xenium_" + XENIUM_CT_SEG + "/{sample_id}/logs/xenium_cell_mask.log"
    container:
        config["containers"]["pixi"]
    shell:
        """
        {PX} python -u {input.script} \
        --seg_dir {params.seg_dir} \
        --out_mask {output.mask} \
        --out_cell_ids {output.cell_ids} \
        &> {log}
        """
