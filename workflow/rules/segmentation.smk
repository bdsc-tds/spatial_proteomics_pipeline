### --- Segmentation: cell masks per method, per-cell quantification, GeoJSON export ---
# Methods come from config segmentation: methods

# --- Rules ---
include: "_segmentation/mesmer.smk"

# --- Map run switches to rules ---
run_rules_segmentation = {
    "segmentation": rules.segmentation_all.input,
}
