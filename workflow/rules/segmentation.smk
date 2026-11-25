### --- Segmentation: masks per method, per-cell quantification, GeoJSON ---
# one rule file per entry of config segmentation: methods (SEG_FILES: Snakefile); the xenium file imports the
# masks of every xenium_{seg} method (SEG_METHODS: Snakefile)

# variables defined in other .smk files:
# SEG_FILES: defined in Snakefile

# --- Rules ---
for _seg_file in SEG_FILES:
    include: f"_segmentation/{_seg_file}.smk"
include: "_segmentation/quantify_cells.smk"

# --- Map run switches to rules ---
run_rules_segmentation = {
    "segmentation": rules.segmentation_all.input,
}
