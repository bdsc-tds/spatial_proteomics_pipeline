### --- Alignment of the proteomics DAPI to the Xenium DAPI of the same slide, and QC plots ---
# Methods come from config xenium_alignment: methods

# --- Rules ---
include: "_xenium_alignment/insitupy.smk"

# --- Map run switches to rules ---
run_rules_xenium_alignment = {
    "xenium_alignment": rules.xenium_alignment_all.input,
}
