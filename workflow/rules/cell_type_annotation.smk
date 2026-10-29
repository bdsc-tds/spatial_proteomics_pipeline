### --- Cell type annotation: cell types per segmentation method (TACIT), TACIT plots over all samples ---
# one rule file per method, included only if listed in config cell_type_annotation: methods (CT_METHODS: Snakefile);
# signatures (TACIT): signatures/{panel}.csv next to the dataset config

# variables defined in other .smk files:
# CT_METHODS: defined in Snakefile

# --- Rules ---
for _ct_method in CT_METHODS:
    include: f"_cell_type_annotation/{_ct_method}.smk"

# --- Map run switches to rules ---
run_rules_cell_type_annotation = {
    "cell_type_annotation": [t for m in CT_METHODS for t in getattr(rules, f"{m}_all").input],
}
