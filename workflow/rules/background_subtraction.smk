### --- Background subtraction (COMET): only for panels with _background_subtraction: true ---
# When switched on, every downstream step reads the background-subtracted image (source_image)

# --- Rules ---
include: "_background_subtraction/backsub.smk"

# --- Map run switches to rules ---
run_rules_background_subtraction = {
    "background_subtraction": rules.background_subtraction_all.input,
}
