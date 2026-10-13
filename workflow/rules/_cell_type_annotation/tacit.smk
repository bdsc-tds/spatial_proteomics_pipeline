# TACIT cell type annotation per segmentation method, then plots over all samples: the mesmer rules
# are the base, other methods derive from them with `use rule` (all inputs restated)

TACIT = config["cell_type_annotation"]["tacit"]


def tacit_signature(sample_id):
    """Signature CSV of the sample's panel: signatures/{panel}.csv next to the dataset config."""
    return os.path.join(os.path.dirname(CONFIG_PATH), "signatures", SAMPLES[sample_id]["panel"] + ".csv")


# samples whose panel has a signature
TACIT_SAMPLE_IDS = [s for s in SAMPLE_IDS if os.path.exists(tacit_signature(s))]
if len(TACIT_SAMPLE_IDS) < len(SAMPLE_IDS):
    print(f"WARNING: TACIT skips {len(SAMPLE_IDS) - len(TACIT_SAMPLE_IDS)} samples whose panel has no signature in "
          f"{os.path.dirname(CONFIG_PATH)}/signatures/", file=sys.stderr)

# cell type of every cell of one sample (cell_id, cell_type, TACIT's UMAP1/UMAP2)
rule tacit_mesmer:
    input:
        script=f"{workflow.basedir}/scripts/_cell_type_annotation/run_tacit.R",
        # from quantify_cells_mesmer
        feature_table=OUT + "/segmentation/mesmer/{sample_id}/feature_table.csv",
        signature=lambda wc: tacit_signature(wc.sample_id),
    output:
        cell_types=OUT + "/cell_type_annotation/tacit/mesmer/{sample_id}/cell_types.csv",
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
    log:
        OUT + "/cell_type_annotation/tacit/cellpose/{sample_id}/logs/tacit.log"

# cell type fractions over all samples; per sample, mean intensity heatmap and preprocess_protein UMAPs
rule plot_tacit_mesmer:
    input:
        script=f"{workflow.basedir}/scripts/_cell_type_annotation/plot_tacit.py",
        # from quantify_cells_mesmer
        feature_tables=expand(OUT + "/segmentation/mesmer/{sample_id}/feature_table.csv", sample_id=TACIT_SAMPLE_IDS),
        cell_types=expand(OUT + "/cell_type_annotation/tacit/mesmer/{sample_id}/cell_types.csv",
                          sample_id=TACIT_SAMPLE_IDS),
        signatures=[tacit_signature(s) for s in TACIT_SAMPLE_IDS],
    output:
        fractions=OUT + "/cell_type_annotation/tacit/mesmer/summary/figures/cell_type_fractions.png",
        fractions_table=OUT + "/cell_type_annotation/tacit/mesmer/summary/figures/cell_type_fractions.csv", # written by savefig_plotnine
        heatmaps=expand(OUT + "/cell_type_annotation/tacit/mesmer/{sample_id}/figures/mean_intensity_heatmap.png",
                        sample_id=TACIT_SAMPLE_IDS),
        umaps=expand(OUT + "/cell_type_annotation/tacit/mesmer/{sample_id}/figures/umap.png", sample_id=TACIT_SAMPLE_IDS),
        umap_coords=expand(OUT + "/cell_type_annotation/tacit/mesmer/{sample_id}/figures/umap.csv",
                           sample_id=TACIT_SAMPLE_IDS), # written by savefig
    params:
        sample_ids=TACIT_SAMPLE_IDS,
        cofactor=TACIT["cofactor"],
        exclude=TACIT["umap_exclude"],
    threads: 1
    resources:
        runtime=720,
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
        --signatures {input.signatures} \
        --cofactor {params.cofactor} \
        --exclude {params.exclude:q} \
        --threads {threads} \
        --out_fractions {output.fractions} \
        --out_heatmaps {output.heatmaps} \
        --out_umaps {output.umaps} \
        &> {log}
        """

use rule plot_tacit_mesmer as plot_tacit_cellpose with:
    input:
        script=f"{workflow.basedir}/scripts/_cell_type_annotation/plot_tacit.py",
        # from quantify_cells_cellpose
        feature_tables=expand(OUT + "/segmentation/cellpose/{sample_id}/feature_table.csv", sample_id=TACIT_SAMPLE_IDS),
        cell_types=expand(OUT + "/cell_type_annotation/tacit/cellpose/{sample_id}/cell_types.csv",
                          sample_id=TACIT_SAMPLE_IDS),
        signatures=[tacit_signature(s) for s in TACIT_SAMPLE_IDS],
    output:
        fractions=OUT + "/cell_type_annotation/tacit/cellpose/summary/figures/cell_type_fractions.png",
        fractions_table=OUT + "/cell_type_annotation/tacit/cellpose/summary/figures/cell_type_fractions.csv", # written by savefig_plotnine
        heatmaps=expand(OUT + "/cell_type_annotation/tacit/cellpose/{sample_id}/figures/mean_intensity_heatmap.png",
                        sample_id=TACIT_SAMPLE_IDS),
        umaps=expand(OUT + "/cell_type_annotation/tacit/cellpose/{sample_id}/figures/umap.png", sample_id=TACIT_SAMPLE_IDS),
        umap_coords=expand(OUT + "/cell_type_annotation/tacit/cellpose/{sample_id}/figures/umap.csv",
                           sample_id=TACIT_SAMPLE_IDS), # written by savefig
    log:
        OUT + "/cell_type_annotation/tacit/cellpose/summary/logs/plot_tacit.log"

# TACIT labels and plots of every method in segmentation: methods
rule cell_type_annotation_all:
    input:
        expand(OUT + "/cell_type_annotation/tacit/{seg_method}/{sample_id}/cell_types.csv",
               seg_method=SEG_METHODS, sample_id=TACIT_SAMPLE_IDS),
        expand(OUT + "/cell_type_annotation/tacit/{seg_method}/summary/figures/cell_type_fractions.png",
               seg_method=SEG_METHODS),
