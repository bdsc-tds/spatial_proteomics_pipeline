### --- Segmentation: masks per method, per-cell quantification, GeoJSON ---
# one rule file per method, included only if listed in config segmentation: methods (SEG_METHODS: Snakefile)

# --- Rules ---
for _seg_method in SEG_METHODS:
    include: f"_segmentation/{_seg_method}.smk"
include: "_segmentation/quantify_cells.smk"

# --- Map run switches to rules ---
run_rules_segmentation = {
    "segmentation": rules.segmentation_all.input,
}
