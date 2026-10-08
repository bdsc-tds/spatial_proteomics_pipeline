# Register the proteomics DAPI onto the Xenium DAPI with insitupy, then match cells and plot per
# segmentation method. Samples with a Xenium run only. SEG_METHODS: segmentation.smk

# {max_width}: px size of the shrunk copies used for feature matching, one folder per value
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
        # Xenium DAPI + registered proteomics DAPI (Xenium px grid)
        registered=insitupy_dir + "dapi_registered.ome.tif",
        dapi_overlay=insitupy_dir + "figures/dapi_overlay.png",
        # each modality alone
        dapi_xenium=insitupy_dir + "figures/dapi_xenium.png",
        dapi_proteomics=insitupy_dir + "figures/dapi_proteomics_registered.png",
    params:
        xenium_dir=xenium_dir,
        min_good_matches=config["xenium_alignment"]["insitupy"]["min_good_matches_per_area"],
    resources:
        # codex: SIFT on the full-size image (~23k² px) runs out of memory at 50G
        runtime=lambda wildcards: 360 if wildcards.max_width == 'codex' else 120,
        mem=lambda wildcards: '200G' if wildcards.max_width == 'codex' else '50G' # adjust for minimal rescaling option
    log:
        insitupy_dir + "logs/align_to_xenium_insitupy.log"
    container:
        config["containers"]["pixi"]
    shell:
        """
        {PX} python -u {input.script} \
        --xenium_dir {params.xenium_dir} \
        --dapi {input.dapi} \
        --out_affine {output.affine} \
        --out_matches {output.matches} \
        --out_registered {output.registered} \
        --out_overlay {output.dapi_overlay} \
        --out_xenium_dapi {output.dapi_xenium} \
        --out_proteomics_dapi {output.dapi_proteomics} \
        --min_good_matches_per_area {params.min_good_matches} \
        --max_width {wildcards.max_width} \
        &> {log}
        """

# match aligned proteomics cells to Xenium cells; overlay and count plots per segmentation method
rule plot_alignment:
    input:
        script=f"{workflow.basedir}/scripts/_xenium_alignment/plot_alignment.py",
        # from quantify_cells_{seg_method}
        feature_table=OUT + "/segmentation/{seg_method}/{sample_id}/feature_table.csv",
        affine=xenium_alignment_dir + "affine_matrix.csv",
        xenium_cells=lambda wc: os.path.join(xenium_dir(wc), "cells.parquet"),
        xenium_experiment=lambda wc: os.path.join(xenium_dir(wc), "experiment.xenium"),
    output:
        matches=xenium_alignment_dir + "{seg_method}_cell_matches.csv",
        overlay=xenium_alignment_dir + "figures/{seg_method}_overlay.png",
        matched_xenium=xenium_alignment_dir + "figures/{seg_method}_matched_spatial_xenium.png",
        matched_proteomics=xenium_alignment_dir + "figures/{seg_method}_matched_spatial_proteomics.png",
        counts=xenium_alignment_dir + "figures/{seg_method}_counts.png",
        counts_table=xenium_alignment_dir + "figures/{seg_method}_counts.csv", # written by savefig_plotnine
    params:
        xenium_dir=xenium_dir,
        max_dist=config["xenium_alignment"]["max_match_dist_um"],
    resources:
        runtime=60,
        mem='16G',
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
        --out_matched_xenium {output.matched_xenium} \
        --out_matched_proteomics {output.matched_proteomics} \
        --out_counts {output.counts} \
        &> {log}
        """

# fraction of Xenium cells matched to a proteomics cell over all samples, segmentation methods and
# max_width values, for the default Xenium segmentation and Proseg (when the run has one)
xenium_alignment_summary_dir = OUT + "/xenium_alignment/{align_method}/summary/"

rule plot_alignment_summary:
    input:
        script=f"{workflow.basedir}/scripts/_xenium_alignment/plot_alignment_summary.py",
        counts=expand(rules.plot_alignment.output.counts_table, max_width=ALIGN_MAX_WIDTHS, seg_method=SEG_METHODS,
                      sample_id=XENIUM_SAMPLE_IDS, allow_missing=True),
    output:
        table=xenium_alignment_summary_dir + "fraction_matched.csv",
        barplot_multimodal=xenium_alignment_summary_dir + "figures/fraction_matched_per_sample_multimodal.png",
        barplot_proseg=xenium_alignment_summary_dir + "figures/fraction_matched_per_sample_proseg.png",
        boxplot=xenium_alignment_summary_dir + "figures/fraction_matched_per_max_width.png",
    params:
        root=OUT + "/xenium_alignment",
        seg_root=OUT + "/segmentation", # feature tables, for the Proseg matching
        experiments=EXPERIMENTS_PATH, # Xenium run of each sample, for the Proseg folder
        max_dist=config["xenium_alignment"]["max_match_dist_um"],
    resources:
        # Proseg cells of every Xenium run read in turn
        runtime=240,
        mem='32G',
    log:
        xenium_alignment_summary_dir + "logs/plot_alignment_summary.log"
    container:
        config["containers"]["pixi"]
    shell:
        """
        {PX} python {input.script} \
        --counts {input.counts:q} \
        --root {params.root} \
        --seg_root {params.seg_root} \
        --experiments {params.experiments} \
        --max_dist {params.max_dist} \
        --out_table {output.table} \
        --out_barplot_multimodal {output.barplot_multimodal} \
        --out_barplot_proseg {output.barplot_proseg} \
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
