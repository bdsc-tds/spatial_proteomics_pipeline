### --- Cell type annotation: TACIT labels per segmentation method, plots over all samples ---
# SEG_METHODS: Snakefile; signatures: signatures/{panel}.csv next to the dataset config

# --- Rules ---
include: "_cell_type_annotation/tacit.smk"

# --- Map run switches to rules ---
run_rules_cell_type_annotation = {
    "cell_type_annotation": rules.cell_type_annotation_all.input,
}
