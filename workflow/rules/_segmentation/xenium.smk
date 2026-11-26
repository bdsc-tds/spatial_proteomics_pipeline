# Xenium segmentations as proteomics segmentation methods (xenium_{seg}): 10x cells warped onto the proteomics
# image grid OUTSIDE the pipeline, imported here (checked against the image, then copied). The masks are then
# quantified, labelled and embedded like every other method (quantify_cells.smk and the downstream steps).
# Methods and samples: SEG_METHOD_INFO (Snakefile), from config segmentation: xenium

# variables defined in other .smk files:
# OUT:               defined in Snakefile
# source_image:      defined in Snakefile
# XENIUM_SEG_CONFIG: defined in Snakefile
# XENIUM_SEG_NAMES:  defined in Snakefile

if XENIUM_SEG_NAMES and not XENIUM_SEG_CONFIG.get("mask_file"):
    raise WorkflowError("segmentation: xenium: mask_file is needed (path template of the warped masks, with {name} "
                        "and {segmentation}, {sample_id} or {sample})")


def xenium_mask_file(name):
    """Input function: file of the external warp output of one sample and Xenium segmentation; mask_file template
    with {name} (file name), {segmentation}, {sample_id} and {sample} (last level of the sample id)."""
    def _file(wc):
        return XENIUM_SEG_CONFIG["mask_file"].format(name=name, segmentation=wc.xenium_seg, sample_id=wc.sample_id,
                                                     sample=wc.sample_id.split("/")[-1])
    return _file


xenium_seg_dir = OUT + "/segmentation/xenium_{xenium_seg}/{sample_id}/"

# warped cell / nuclei / cytoplasm masks and their label -> 10x cell_id table, checked against the image
rule import_xenium_mask:
    input:
        script=f"{workflow.basedir}/scripts/_segmentation/import_xenium_mask.py",
        # external warp of the 10x masks onto the proteomics image grid
        cells=xenium_mask_file("cell_boundaries_mask.tif"),
        nuclei=xenium_mask_file("nuclei_boundaries_mask.tif"),
        cytoplasm=xenium_mask_file("cytoplasm_boundaries_mask.tif"),
        cell_ids=xenium_mask_file("cell_ids.csv"),
        # raw or background-subtracted image (source_image, Snakefile): its size only
        image=source_image,
    output:
        cells=xenium_seg_dir + "cell_boundaries_mask.tif",
        nuclei=xenium_seg_dir + "nuclei_boundaries_mask.tif",
        cytoplasm=xenium_seg_dir + "cytoplasm_boundaries_mask.tif",
        cell_ids=xenium_seg_dir + "cell_ids.csv",
    wildcard_constraints:
        xenium_seg="|".join(re.escape(s) for s in XENIUM_SEG_NAMES) or "(?!)",  # no segmentation: matches nothing
    threads: 1
    resources:
        # one mask row block at a time
        runtime=60,
        mem='16G',
    log:
        xenium_seg_dir + "logs/import_xenium_mask.log"
    container:
        config["containers"]["pixi"]
    shell:
        """
        {PX} python -u {input.script} \
        --image {input.image} \
        --cells {input.cells} \
        --nuclei {input.nuclei} \
        --cytoplasm {input.cytoplasm} \
        --cell_ids {input.cell_ids} \
        --out_cells {output.cells} \
        --out_nuclei {output.nuclei} \
        --out_cytoplasm {output.cytoplasm} \
        --out_cell_ids {output.cell_ids} \
        &> {log}
        """
