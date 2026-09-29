# rules to register the proteomics DAPI image onto the Xenium DAPI (nuclei) image with insitupy, then
# overlay the aligned proteomics cells on the Xenium cells and count matched cells.
# Only samples with a Xenium run (exact sample-name match, see experiments.py) get targets

xenium_alignment_dir = OUT + "/xenium_alignment/{align_method}/{sample_id}/"
insitupy_dir = OUT + "/xenium_alignment/insitupy/{sample_id}/"


def xenium_dir(wildcards):
    return SAMPLES[wildcards.sample_id]["xenium"]


rule align_to_xenium_insitupy:
    input:
        script=f"{workflow.basedir}/scripts/_xenium_alignment/insitupy_alignment.py",
        dapi=rules.prepare_image.output.dapi,
        xenium=lambda wc: os.path.join(xenium_dir(wc), "experiment.xenium"),
    output:
        affine=insitupy_dir + "affine_matrix.csv",
    params:
        xenium_dir=xenium_dir,
        min_good_matches=config["xenium_alignment"]["min_good_matches_per_area"],
    threads: 4
    resources:
        runtime=480,
        mem='64G',
    log:
        insitupy_dir + "logs/align_to_xenium_insitupy.log"
    container:
        config["containers"]["pixi"]
    shell:
        """
        {PX} python {input.script} \
        {params.xenium_dir} \
        {input.dapi} \
        $(dirname {output.affine}) \
        --min_good_matches_per_area {params.min_good_matches} \
        &> {log}
        """

# overlay of aligned proteomics and Xenium cells, per segmentation method
rule plot_alignment:
    input:
        script=f"{workflow.basedir}/scripts/_xenium_alignment/plot_alignment.py",
        feature_table=segmentation_dir + "feature_table.csv",
        affine=xenium_alignment_dir + "affine_matrix.csv",
        xenium_cells=lambda wc: os.path.join(xenium_dir(wc), "cells.parquet"),
    output:
        overlay=xenium_alignment_dir + "figures/{seg_method}_overlay.png",
        counts=xenium_alignment_dir + "figures/{seg_method}_counts.png",
    params:
        pixel_size=lambda wc: panel_meta(wc)["pixel_size_um"],
        xenium_pixel_size=config["xenium_alignment"]["xenium_pixel_size_um"],
        max_dist=config["xenium_alignment"]["max_match_dist_um"],
    resources:
        runtime=120,
        mem='32G',
    log:
        xenium_alignment_dir + "logs/{seg_method}_plot_alignment.log"
    container:
        config["containers"]["pixi"]
    shell:
        """
        {PX} python {input.script} \
        --feature_table {input.feature_table} \
        --affine {input.affine} \
        --xenium_cells {input.xenium_cells} \
        --pixel_size {params.pixel_size} \
        --xenium_pixel_size {params.xenium_pixel_size} \
        --max_dist {params.max_dist} \
        --out_overlay {output.overlay} \
        --out_counts {output.counts} \
        &> {log}
        """

rule xenium_alignment_all:
    input:
        expand(xenium_alignment_dir + "affine_matrix.csv", align_method=ALIGN_METHODS, sample_id=XENIUM_SAMPLE_IDS),
        expand(rules.plot_alignment.output, align_method=ALIGN_METHODS, seg_method=SEG_METHODS, sample_id=XENIUM_SAMPLE_IDS),
