# rules to register the proteomics DAPI image onto the Xenium DAPI (nuclei) image with insitupy, then
# overlay the aligned proteomics cells on the Xenium cells and count matched cells.
# Only samples with a Xenium run (donor + sample name match, see experiments.match_xenium_runs) get targets

# {max_width}: size (px) of the shrunk copies insitupy matches features on, one output folder per value
xenium_alignment_dir = OUT + "/xenium_alignment/{align_method}/max_width_{max_width}/{sample_id}/"
insitupy_dir = OUT + "/xenium_alignment/insitupy/max_width_{max_width}/{sample_id}/"


def xenium_dir(wildcards):
    return SAMPLES[wildcards.sample_id]["xenium"]


rule align_to_xenium_insitupy:
    input:
        script=f"{workflow.basedir}/scripts/_xenium_alignment/insitupy_alignment.py",
        dapi=rules.prepare_image.output.dapi,
        xenium=lambda wc: os.path.join(xenium_dir(wc), "experiment.xenium"),
    output:
        affine=insitupy_dir + "affine_matrix.csv",
        matches=insitupy_dir + "figures/feature_matches.png",
    params:
        xenium_dir=xenium_dir,
        min_good_matches=config["xenium_alignment"]["insitupy"]["min_good_matches_per_area"],
    resources:
        # proteomics DAPI and Xenium nuclei images in memory, features matched on max_width px copies
        runtime=360,
        mem='50G',
    log:
        insitupy_dir + "logs/align_to_xenium_insitupy.log"
    container:
        config["containers"]["pixi"]
    shell:
        """
        {PX} python {input.script} \
        --xenium_dir {params.xenium_dir} \
        --dapi {input.dapi} \
        --out_affine {output.affine} \
        --out_matches {output.matches} \
        --min_good_matches_per_area {params.min_good_matches} \
        --max_width {wildcards.max_width} \
        &> {log}
        """

# aligned proteomics cells matched to Xenium cells, with overlay and count plots, per segmentation method
rule plot_alignment:
    input:
        script=f"{workflow.basedir}/scripts/_xenium_alignment/plot_alignment.py",
        feature_table=rules.quantify_cells.output.feature_table, # output from mesmer rule
        affine=xenium_alignment_dir + "affine_matrix.csv",
        xenium_cells=lambda wc: os.path.join(xenium_dir(wc), "cells.parquet"),
        xenium_experiment=lambda wc: os.path.join(xenium_dir(wc), "experiment.xenium"),
    output:
        matches=xenium_alignment_dir + "{seg_method}_cell_matches.csv",
        overlay=xenium_alignment_dir + "figures/{seg_method}_overlay.png",
        counts=xenium_alignment_dir + "figures/{seg_method}_counts.png",
        counts_table=xenium_alignment_dir + "figures/{seg_method}_counts.csv", # written by savefig_plotnine
    params:
        xenium_dir=xenium_dir,
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
        --xenium_dir {params.xenium_dir} \
        --max_dist {params.max_dist} \
        --sample_id {wildcards.sample_id} \
        --out_matches {output.matches} \
        --out_overlay {output.overlay} \
        --out_counts {output.counts} \
        &> {log}
        """

# fraction of Xenium cells matched to a proteomics cell, per alignment method, over all samples,
# segmentation methods and max_width values
xenium_alignment_summary_dir = OUT + "/xenium_alignment/{align_method}/summary/"

rule plot_alignment_summary:
    input:
        script=f"{workflow.basedir}/scripts/_xenium_alignment/plot_alignment_summary.py",
        counts=expand(rules.plot_alignment.output.counts_table, max_width=ALIGN_MAX_WIDTHS, seg_method=SEG_METHODS,
                      sample_id=XENIUM_SAMPLE_IDS, allow_missing=True),
    output:
        table=xenium_alignment_summary_dir + "fraction_matched.csv",
        barplot=xenium_alignment_summary_dir + "figures/fraction_matched_per_sample.png",
        boxplot=xenium_alignment_summary_dir + "figures/fraction_matched_per_max_width.png",
    params:
        root=OUT + "/xenium_alignment",
    resources:
        runtime=30,
        mem='4G',
    log:
        xenium_alignment_summary_dir + "logs/plot_alignment_summary.log"
    container:
        config["containers"]["pixi"]
    shell:
        """
        {PX} python {input.script} \
        --counts {input.counts:q} \
        --root {params.root} \
        --out_table {output.table} \
        --out_barplot {output.barplot} \
        --out_boxplot {output.boxplot} \
        &> {log}
        """

rule xenium_alignment_all:
    input:
        expand(xenium_alignment_dir + "affine_matrix.csv", align_method=ALIGN_METHODS, max_width=ALIGN_MAX_WIDTHS,
               sample_id=XENIUM_SAMPLE_IDS),
        expand(rules.plot_alignment.output, align_method=ALIGN_METHODS, max_width=ALIGN_MAX_WIDTHS,
               seg_method=SEG_METHODS, sample_id=XENIUM_SAMPLE_IDS),
        expand(rules.plot_alignment_summary.output, align_method=ALIGN_METHODS),
