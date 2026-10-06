### --- Alignment of the proteomics DAPI to the Xenium DAPI of the same slide, and QC plots ---
# Methods come from config xenium_alignment: methods

ALIGN_METHODS = config["xenium_alignment"]["methods"]
# max_width values to align with (one output folder each), to compare registrations
ALIGN_MAX_WIDTHS = config["xenium_alignment"]["insitupy"]["max_width"]

wildcard_constraints:
    align_method="|".join(re.escape(m) for m in ALIGN_METHODS),
    max_width=r"\d+|codex",

# --- Rules ---
include: "_xenium_alignment/insitupy.smk"

# --- Map run switches to rules ---
run_rules_xenium_alignment = {
    "xenium_alignment": rules.xenium_alignment_all.input,
}
