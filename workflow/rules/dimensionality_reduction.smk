### --- Dimensionality reduction: per-cell embeddings of foundation models, for every segmentation method's mask ---
# one rule file per method, included only if listed in config dimensionality_reduction: methods (DR_METHODS: Snakefile);
# samples imaged on Xenium also embed 10x's own cells, exported as a label mask first (xenium_cells.smk)

# variables defined in other .smk files:
# DR_METHODS: defined in Snakefile

# --- Rules ---
include: "_dimensionality_reduction/xenium_cells.smk"
for _dr_method in DR_METHODS:
    include: f"_dimensionality_reduction/{_dr_method}.smk"

# --- Map run switches to rules ---
run_rules_dimensionality_reduction = {
    "dimensionality_reduction": [t for m in DR_METHODS for t in getattr(rules, f"{m}_embeddings_all").input],
}
