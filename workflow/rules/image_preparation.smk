### --- Image preparation: channel table, DAPI image and segmentation input per sample ---

# --- Rules ---
include: "_image_preparation/prepare_image.smk"

# --- Map run switches to rules ---
run_rules_image_preparation = {
    "image_preparation": rules.image_preparation_all.input,
}
