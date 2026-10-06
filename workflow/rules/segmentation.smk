### --- Segmentation: cell masks per method, per-cell quantification, GeoJSON export ---
# Methods come from config segmentation: methods

SEG_METHODS = config["segmentation"]["methods"]

wildcard_constraints:
    seg_method="|".join(re.escape(m) for m in SEG_METHODS),

# --- Rules ---
include: "_segmentation/mesmer.smk"

# --- Map run switches to rules ---
run_rules_segmentation = {
    "segmentation": rules.segmentation_all.input,
}
