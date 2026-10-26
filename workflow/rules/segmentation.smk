### --- Segmentation: masks per method, per-cell quantification, GeoJSON ---
# one rule file per method, included only if listed in config segmentation: methods (SEG_METHODS: Snakefile);
# VirTues' also when only cell_type_annotation: methods lists it (its class map labels every mask)

# --- Rules ---
for _seg_method in SEG_METHODS + [m for m in CT_METHODS if m == "virtues" and m not in SEG_METHODS]:
    include: f"_segmentation/{_seg_method}.smk"
include: "_segmentation/quantify_cells.smk"

# --- Map run switches to rules ---
run_rules_segmentation = {
    "segmentation": rules.segmentation_all.input,
}
