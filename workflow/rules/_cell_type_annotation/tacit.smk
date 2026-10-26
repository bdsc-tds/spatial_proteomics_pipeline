# TACIT cell type annotation per segmentation method, then plots over all samples: the mesmer rules
# are the base, other methods derive from them with `use rule` (all inputs restated). Samples imaged
# on Xenium also get TACIT on 10x's own protein intensities (xenium_{xenium_segmentation} folder)

TACIT = config["cell_type_annotation"]["tacit"]


def tacit_signature(sample_id):
    """Signature CSV of the sample's panel: signatures/{panel}.csv next to the dataset config."""
    return os.path.join(os.path.dirname(CONFIG_PATH), "signatures", SAMPLES[sample_id]["panel"] + ".csv")


# samples whose panel has a signature
TACIT_SAMPLE_IDS = [s for s in SAMPLE_IDS if os.path.exists(tacit_signature(s))]
if len(TACIT_SAMPLE_IDS) < len(SAMPLE_IDS):
    print(f"WARNING: TACIT skips {len(SAMPLE_IDS) - len(TACIT_SAMPLE_IDS)} samples whose panel has no signature in "
          f"{os.path.dirname(CONFIG_PATH)}/signatures/", file=sys.stderr)

# temporary to run on one sample
TACIT_SAMPLE_IDS=TACIT_SAMPLE_IDS

# samples imaged on Xenium (technology xenium) whose panel has a signature and whose run has the
# XENIUM_CT_SEG segmentation (xenium_ct_seg_dir: Snakefile)
TACIT_XENIUM_SAMPLE_IDS = [s for s in TACIT_SAMPLE_IDS
                           if SAMPLES[s]["panel_meta"]["technology"] == "xenium" and xenium_ct_seg_dir(s)]

# cell type of every cell of one sample (cell_id, cell_type, TACIT's UMAP1/UMAP2)
rule tacit_mesmer:
    input:
        script=f"{workflow.basedir}/scripts/_cell_type_annotation/run_tacit.R",
        # from quantify_cells_mesmer
        feature_table=OUT + "/segmentation/mesmer/{sample_id}/feature_table.csv",
        signature=lambda wc: tacit_signature(wc.sample_id),
    output:
        cell_types=OUT + "/cell_type_annotation/tacit/mesmer/{sample_id}/cell_types.csv",
        ctr=OUT + "/cell_type_annotation/tacit/mesmer/{sample_id}/ctr.csv",
        thresholds=OUT + "/cell_type_annotation/tacit/mesmer/{sample_id}/thresholds.csv",
        knn=OUT + "/cell_type_annotation/tacit/mesmer/{sample_id}/knn_fractions.csv",
    params:
        cofactor=TACIT["cofactor"],
        r=TACIT["r"],
        p=TACIT["p"],
    threads: 2
    resources:
        runtime=240,
        mem='64G',
    log:
        OUT + "/cell_type_annotation/tacit/mesmer/{sample_id}/logs/tacit.log"
    container:
        config["containers"]["pixi"]
    shell:
        """
        {PX} -e tacit Rscript {input.script} \
        --feature_table {input.feature_table} \
        --signature {input.signature} \
        --cofactor {params.cofactor} \
        --r {params.r} \
        --p {params.p} \
        --threads {threads} \
        --out_csv {output.cell_types} \
        --out_ctr {output.ctr} \
        --out_thresholds {output.thresholds} \
        --out_knn {output.knn} \
        &> {log}
        """

use rule tacit_mesmer as tacit_cellpose with:
    input:
        script=f"{workflow.basedir}/scripts/_cell_type_annotation/run_tacit.R",
        # from quantify_cells_cellpose
        feature_table=OUT + "/segmentation/cellpose/{sample_id}/feature_table.csv",
        signature=lambda wc: tacit_signature(wc.sample_id),
    output:
        cell_types=OUT + "/cell_type_annotation/tacit/cellpose/{sample_id}/cell_types.csv",
        ctr=OUT + "/cell_type_annotation/tacit/cellpose/{sample_id}/ctr.csv",
        thresholds=OUT + "/cell_type_annotation/tacit/cellpose/{sample_id}/thresholds.csv",
        knn=OUT + "/cell_type_annotation/tacit/cellpose/{sample_id}/knn_fractions.csv",
    log:
        OUT + "/cell_type_annotation/tacit/cellpose/{sample_id}/logs/tacit.log"

use rule tacit_mesmer as tacit_virtues with:
    input:
        script=f"{workflow.basedir}/scripts/_cell_type_annotation/run_tacit.R",
        # from quantify_cells_virtues
        feature_table=OUT + "/segmentation/virtues/{sample_id}/feature_table.csv",
        signature=lambda wc: tacit_signature(wc.sample_id),
    output:
        cell_types=OUT + "/cell_type_annotation/tacit/virtues/{sample_id}/cell_types.csv",
        ctr=OUT + "/cell_type_annotation/tacit/virtues/{sample_id}/ctr.csv",
        thresholds=OUT + "/cell_type_annotation/tacit/virtues/{sample_id}/thresholds.csv",
        knn=OUT + "/cell_type_annotation/tacit/virtues/{sample_id}/knn_fractions.csv",
    log:
        OUT + "/cell_type_annotation/tacit/virtues/{sample_id}/logs/tacit.log"

# 10x protein intensities of the Xenium cells (read_xenium_sample) in the quantify_cells feature table format
rule xenium_protein_table:
    input:
        script=f"{workflow.basedir}/scripts/_cell_type_annotation/xenium_protein_table.py",
    output:
        feature_table=OUT + "/cell_type_annotation/tacit/xenium_" + XENIUM_CT_SEG + "/{sample_id}/feature_table.csv",
    params:
        # XENIUM_CT_SEG segmentation of the Xenium run (normalised_results/outs)
        xenium_dir=lambda wc: xenium_ct_seg_dir(wc.sample_id),
    threads: 1
    resources:
        runtime=60,
        mem='32G',
    log:
        OUT + "/cell_type_annotation/tacit/xenium_" + XENIUM_CT_SEG + "/{sample_id}/logs/xenium_protein_table.log"
    container:
        config["containers"]["pixi"]
    shell:
        """
        {PX} python -u {input.script} \
        --xenium_dir {params.xenium_dir} \
        --out_csv {output.feature_table} \
        &> {log}
        """

use rule tacit_mesmer as tacit_xenium with:
    input:
        script=f"{workflow.basedir}/scripts/_cell_type_annotation/run_tacit.R",
        # from xenium_protein_table
        feature_table=OUT + "/cell_type_annotation/tacit/xenium_" + XENIUM_CT_SEG + "/{sample_id}/feature_table.csv",
        signature=lambda wc: tacit_signature(wc.sample_id),
    output:
        cell_types=OUT + "/cell_type_annotation/tacit/xenium_" + XENIUM_CT_SEG + "/{sample_id}/cell_types.csv",
        ctr=OUT + "/cell_type_annotation/tacit/xenium_" + XENIUM_CT_SEG + "/{sample_id}/ctr.csv",
        thresholds=OUT + "/cell_type_annotation/tacit/xenium_" + XENIUM_CT_SEG + "/{sample_id}/thresholds.csv",
        knn=OUT + "/cell_type_annotation/tacit/xenium_" + XENIUM_CT_SEG + "/{sample_id}/knn_fractions.csv",
    params:
        cofactor=TACIT["xenium_cofactor"],
        r=TACIT["r"],
        p=TACIT["p"],
    log:
        OUT + "/cell_type_annotation/tacit/xenium_" + XENIUM_CT_SEG + "/{sample_id}/logs/tacit.log"

# cell type and mixed / unique summaries over all samples; per sample, the same as bar plots and a mean intensity
# heatmap. Only the cell type fraction plots are declared, the other figures are written next to them
rule plot_tacit_mesmer:
    input:
        script=f"{workflow.basedir}/scripts/_cell_type_annotation/plot_tacit.py",
        # from quantify_cells_mesmer
        feature_tables=expand(OUT + "/segmentation/mesmer/{sample_id}/feature_table.csv", sample_id=TACIT_SAMPLE_IDS),
        cell_types=expand(OUT + "/cell_type_annotation/tacit/mesmer/{sample_id}/cell_types.csv",
                          sample_id=TACIT_SAMPLE_IDS),
        knn_fractions=expand(OUT + "/cell_type_annotation/tacit/mesmer/{sample_id}/knn_fractions.csv",
                             sample_id=TACIT_SAMPLE_IDS),
        signatures=[tacit_signature(s) for s in TACIT_SAMPLE_IDS],
    output:
        fractions=OUT + "/cell_type_annotation/tacit/mesmer/summary/figures/cell_type_fractions.png",
        sample_fractions=expand(OUT + "/cell_type_annotation/tacit/mesmer/{sample_id}/figures/cell_type_fractions.png",
                                sample_id=TACIT_SAMPLE_IDS),
    params:
        sample_ids=TACIT_SAMPLE_IDS,
        cofactor=TACIT["cofactor"],
    threads: 1
    resources:
        runtime=120,
        mem='64G',
    log:
        OUT + "/cell_type_annotation/tacit/mesmer/summary/logs/plot_tacit.log"
    container:
        config["containers"]["pixi"]
    shell:
        """
        {PX} python -u {input.script} \
        --sample_ids {params.sample_ids:q} \
        --feature_tables {input.feature_tables} \
        --cell_types {input.cell_types} \
        --knn_fractions {input.knn_fractions} \
        --signatures {input.signatures} \
        --cofactor {params.cofactor} \
        --threads {threads} \
        --out_fractions {output.fractions} \
        --out_sample_fractions {output.sample_fractions} \
        &> {log}
        """

use rule plot_tacit_mesmer as plot_tacit_cellpose with:
    input:
        script=f"{workflow.basedir}/scripts/_cell_type_annotation/plot_tacit.py",
        # from quantify_cells_cellpose
        feature_tables=expand(OUT + "/segmentation/cellpose/{sample_id}/feature_table.csv", sample_id=TACIT_SAMPLE_IDS),
        cell_types=expand(OUT + "/cell_type_annotation/tacit/cellpose/{sample_id}/cell_types.csv",
                          sample_id=TACIT_SAMPLE_IDS),
        knn_fractions=expand(OUT + "/cell_type_annotation/tacit/cellpose/{sample_id}/knn_fractions.csv",
                             sample_id=TACIT_SAMPLE_IDS),
        signatures=[tacit_signature(s) for s in TACIT_SAMPLE_IDS],
    output:
        fractions=OUT + "/cell_type_annotation/tacit/cellpose/summary/figures/cell_type_fractions.png",
        sample_fractions=expand(OUT + "/cell_type_annotation/tacit/cellpose/{sample_id}/figures/cell_type_fractions.png",
                                sample_id=TACIT_SAMPLE_IDS),
    log:
        OUT + "/cell_type_annotation/tacit/cellpose/summary/logs/plot_tacit.log"

use rule plot_tacit_mesmer as plot_tacit_virtues with:
    input:
        script=f"{workflow.basedir}/scripts/_cell_type_annotation/plot_tacit.py",
        # from quantify_cells_virtues
        feature_tables=expand(OUT + "/segmentation/virtues/{sample_id}/feature_table.csv", sample_id=TACIT_SAMPLE_IDS),
        cell_types=expand(OUT + "/cell_type_annotation/tacit/virtues/{sample_id}/cell_types.csv",
                          sample_id=TACIT_SAMPLE_IDS),
        knn_fractions=expand(OUT + "/cell_type_annotation/tacit/virtues/{sample_id}/knn_fractions.csv",
                             sample_id=TACIT_SAMPLE_IDS),
        signatures=[tacit_signature(s) for s in TACIT_SAMPLE_IDS],
    output:
        fractions=OUT + "/cell_type_annotation/tacit/virtues/summary/figures/cell_type_fractions.png",
        sample_fractions=expand(OUT + "/cell_type_annotation/tacit/virtues/{sample_id}/figures/cell_type_fractions.png",
                                sample_id=TACIT_SAMPLE_IDS),
    log:
        OUT + "/cell_type_annotation/tacit/virtues/summary/logs/plot_tacit.log"

use rule plot_tacit_mesmer as plot_tacit_xenium with:
    input:
        script=f"{workflow.basedir}/scripts/_cell_type_annotation/plot_tacit.py",
        # from xenium_protein_table
        feature_tables=expand(OUT + "/cell_type_annotation/tacit/xenium_" + XENIUM_CT_SEG + "/{sample_id}/feature_table.csv",
                              sample_id=TACIT_XENIUM_SAMPLE_IDS),
        cell_types=expand(OUT + "/cell_type_annotation/tacit/xenium_" + XENIUM_CT_SEG + "/{sample_id}/cell_types.csv",
                          sample_id=TACIT_XENIUM_SAMPLE_IDS),
        knn_fractions=expand(OUT + "/cell_type_annotation/tacit/xenium_" + XENIUM_CT_SEG + "/{sample_id}/knn_fractions.csv",
                             sample_id=TACIT_XENIUM_SAMPLE_IDS),
        signatures=[tacit_signature(s) for s in TACIT_XENIUM_SAMPLE_IDS],
    output:
        fractions=OUT + "/cell_type_annotation/tacit/xenium_" + XENIUM_CT_SEG + "/summary/figures/cell_type_fractions.png",
        sample_fractions=expand(OUT + "/cell_type_annotation/tacit/xenium_" + XENIUM_CT_SEG + "/{sample_id}/figures/cell_type_fractions.png",
                                sample_id=TACIT_XENIUM_SAMPLE_IDS),
    params:
        sample_ids=TACIT_XENIUM_SAMPLE_IDS,
        cofactor=TACIT["xenium_cofactor"],
    log:
        OUT + "/cell_type_annotation/tacit/xenium_" + XENIUM_CT_SEG + "/summary/logs/plot_tacit.log"

# TACIT labels and plots of every method in segmentation: methods, and of the Xenium cells (samples imaged on Xenium)
rule tacit_all:
    input:
        expand(OUT + "/cell_type_annotation/tacit/{seg_method}/{sample_id}/cell_types.csv",
               seg_method=SEG_METHODS, sample_id=TACIT_SAMPLE_IDS),
        expand(OUT + "/cell_type_annotation/tacit/{seg_method}/summary/figures/cell_type_fractions.png",
               seg_method=SEG_METHODS),
        expand(OUT + "/cell_type_annotation/tacit/xenium_" + XENIUM_CT_SEG + "/{sample_id}/cell_types.csv",
               sample_id=TACIT_XENIUM_SAMPLE_IDS),
        # no plot without Xenium samples (plot_tacit.py needs at least one)
        [OUT + "/cell_type_annotation/tacit/xenium_" + XENIUM_CT_SEG + "/summary/figures/cell_type_fractions.png"]
        if TACIT_XENIUM_SAMPLE_IDS else [],
