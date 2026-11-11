# Register the proteomics DAPI onto the Xenium DAPI with insitupy (affine, optional local refinement), then match
# cells and plot per segmentation method. Samples with a Xenium run only. SEG_METHODS: Snakefile

# variables defined in other .smk files:
# os:                   imported in Snakefile
# OUT:                  defined in Snakefile
# SAMPLES:              defined in Snakefile
# xenium_segmentations: imported in Snakefile (experiments.py)
# XENIUM_SAMPLE_IDS:    defined in Snakefile
# ALIGN_MAX_WIDTHS:     defined in xenium_alignment.smk
# SEG_METHODS:          defined in Snakefile
# EXPERIMENTS_PATH:     defined in Snakefile
# ALIGN_METHODS:        defined in xenium_alignment.smk

# {max_width}: px size of the shrunk copies used for feature matching, one folder per value
xenium_alignment_dir = OUT + "/xenium_alignment/{align_method}/max_width_{max_width}/{sample_id}/"
insitupy_dir = OUT + "/xenium_alignment/insitupy/max_width_{max_width}/{sample_id}/"

# local refinement after the affine (run, grid, work_max_px, max_shift_um); when run, the cell matching uses it
REFINE = config["xenium_alignment"]["insitupy"]["refine"]


def xenium_dir(wildcards):
    return SAMPLES[wildcards.sample_id]["xenium"]


# Xenium segmentations of each sample, as named in the Xenium processed/segmentation folder:
# {sample_id: {xenium_seg: (output dir, marker file)}}
XENIUM_SEGS = {s: xenium_segmentations(SAMPLES[s]["xenium"]) for s in XENIUM_SAMPLE_IDS}


rule align_to_xenium_insitupy:
    input:
        script=f"{workflow.basedir}/scripts/_xenium_alignment/insitupy_alignment.py",
        # from prepare_dapi
        dapi=OUT + "/image_preparation/{sample_id}/dapi.ome.tif",
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
        runtime=lambda wildcards: 45 if wildcards.max_width == 'codex' else 30,
        mem=lambda wildcards: {'codex': '160G', '16000': '100G'}.get(wildcards.max_width, '40G')
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

# local refinement: residual shift per tile after the affine (SIFT + RANSAC on both DAPIs), as control points
rule refine_alignment_insitupy:
    input:
        script=f"{workflow.basedir}/scripts/_xenium_alignment/local_refinement.py",
        # from align_to_xenium_insitupy
        registered=insitupy_dir + "dapi_registered.ome.tif",
        xenium=lambda wc: os.path.join(xenium_dir(wc), "experiment.xenium"),
    output:
        control_points=insitupy_dir + "local_refinement/control_points.csv",
        figure=insitupy_dir + "figures/local_refinement.png",
    params:
        xenium_dir=xenium_dir,
        grid=REFINE["grid"],
        work_max_px=REFINE["work_max_px"],
        max_shift_um=REFINE["max_shift_um"],
        max_neighbour_dev_um=REFINE["max_neighbour_dev_um"],
    resources:
        runtime=30,
        mem='16G',
    log:
        insitupy_dir + "logs/refine_alignment_insitupy.log"
    container:
        config["containers"]["pixi"]
    shell:
        """
        {PX} python -u {input.script} \
        --registered {input.registered} \
        --xenium_dir {params.xenium_dir} \
        --grid {params.grid} \
        --work_max_px {params.work_max_px} \
        --max_shift_um {params.max_shift_um} \
        --max_neighbour_dev_um {params.max_neighbour_dev_um} \
        --sample_id {wildcards.sample_id} \
        --out_control_points {output.control_points} \
        --out_figure {output.figure} \
        &> {log}
        """

# match aligned proteomics cells to Xenium cells; overlay and count plots per segmentation method
rule plot_alignment:
    input:
        script=f"{workflow.basedir}/scripts/_xenium_alignment/plot_alignment.py",
        # from quantify_cells_{seg_method}
        feature_table=OUT + "/segmentation/{seg_method}/{sample_id}/feature_table.csv",
        affine=xenium_alignment_dir + "affine_matrix.csv",
        # from refine_alignment_insitupy, when refine: run
        control_points=[xenium_alignment_dir + "local_refinement/control_points.csv"] if REFINE["run"] else [],
        xenium_cells=lambda wc: os.path.join(xenium_dir(wc), "cells.parquet"),
        xenium_experiment=lambda wc: os.path.join(xenium_dir(wc), "experiment.xenium"),
        # other Xenium segmentations of the sample (processed folder), cell numbers only
        xenium_segs=lambda wc: [marker for _, marker in XENIUM_SEGS[wc.sample_id].values()],
    output:
        matches=xenium_alignment_dir + "{seg_method}/cell_matches.csv",
        overlay=xenium_alignment_dir + "figures/{seg_method}/overlay.png",
        matched_xenium=xenium_alignment_dir + "figures/{seg_method}/matched_spatial_xenium.png",
        matched_proteomics=xenium_alignment_dir + "figures/{seg_method}/matched_spatial_proteomics.png",
        counts=xenium_alignment_dir + "figures/{seg_method}/counts.png",
        counts_table=xenium_alignment_dir + "figures/{seg_method}/counts.csv", # written by savefig_plotnine
        fraction_matched=xenium_alignment_dir + "figures/{seg_method}/fraction_matched.png",
    params:
        xenium_dir=xenium_dir,
        # name=output dir of each other Xenium segmentation
        xenium_segs=lambda wc: [f"{seg}={d}" for seg, (d, _) in XENIUM_SEGS[wc.sample_id].items()],
        max_dist=config["xenium_alignment"]["max_match_dist_um"],
        control_points=lambda wc, input: f"--control_points {input.control_points}" if REFINE["run"] else "",
    resources:
        runtime=30,
        mem='8G',
    log:
        xenium_alignment_dir + "logs/{seg_method}_plot_alignment.log"
    container:
        config["containers"]["pixi"]
    shell:
        """
        {PX} python {input.script} \
        --feature_table {input.feature_table} \
        --affine {input.affine} \
        {params.control_points} \
        --xenium_dir {params.xenium_dir} \
        --xenium_segs {params.xenium_segs:q} \
        --max_dist {params.max_dist} \
        --sample_id {wildcards.sample_id} \
        --out_matches {output.matches} \
        --out_overlay {output.overlay} \
        --out_matched_xenium {output.matched_xenium} \
        --out_matched_proteomics {output.matched_proteomics} \
        --out_counts {output.counts} \
        --out_fraction_matched {output.fraction_matched} \
        &> {log}
        """

# fraction of Xenium cells matched to a proteomics cell over all samples, segmentation methods and
# max_width values, for the default Xenium segmentation and Proseg (when the run has one)
xenium_alignment_summary_dir = OUT + "/xenium_alignment/{align_method}/summary/"
# default Xenium segmentation (matched in plot_alignment) and Proseg
SUMMARY_XENIUM_SEGS = [config["xenium_alignment"]["xenium_segmentation"], "proseg"]

rule plot_alignment_summary:
    input:
        script=f"{workflow.basedir}/scripts/_xenium_alignment/plot_alignment_summary.py",
        counts=expand(rules.plot_alignment.output.counts_table, max_width=ALIGN_MAX_WIDTHS, seg_method=SEG_METHODS,
                      sample_id=XENIUM_SAMPLE_IDS, allow_missing=True),
    output:
        table=xenium_alignment_summary_dir + "fraction_matched.csv",
        # one per Xenium segmentation / proteomics segmentation, file names set in the script
        barplots=expand(xenium_alignment_summary_dir + "figures/fraction_matched_per_sample_{xenium_seg}_{seg_method}.png",
                        xenium_seg=SUMMARY_XENIUM_SEGS, seg_method=SEG_METHODS, allow_missing=True),
        # one per Xenium segmentation, the proteomics segmentations side by side at summary_max_width
        barplots_seg_methods=expand(xenium_alignment_summary_dir + "figures/fraction_matched_per_seg_method_{xenium_seg}.png",
                                    xenium_seg=SUMMARY_XENIUM_SEGS, allow_missing=True),
        boxplot=xenium_alignment_summary_dir + "figures/fraction_matched_per_max_width.png",
    params:
        root=OUT + "/xenium_alignment",
        seg_root=OUT + "/segmentation", # feature tables, for the Proseg matching
        experiments=EXPERIMENTS_PATH, # Xenium run of each sample, for the Proseg folder
        max_dist=config["xenium_alignment"]["max_match_dist_um"],
        xenium_seg=config["xenium_alignment"]["xenium_segmentation"],
        summary_max_width=config["xenium_alignment"]["insitupy"]["summary_max_width"],
        refine="--refine" if REFINE["run"] else "", # Proseg matching with the local refinement
        figures_dir=xenium_alignment_summary_dir + "figures",
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
        --xenium_seg {params.xenium_seg} \
        --summary_max_width {params.summary_max_width} \
        {params.refine} \
        --out_table {output.table} \
        --out_figures_dir {params.figures_dir} \
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
        expand(rules.refine_alignment_insitupy.output, max_width=ALIGN_MAX_WIDTHS,
               sample_id=XENIUM_SAMPLE_IDS) if REFINE["run"] else [],
