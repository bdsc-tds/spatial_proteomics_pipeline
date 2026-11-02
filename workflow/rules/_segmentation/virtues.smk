# VirTues (github.com/bunnelab/virtues) whole-cell masks, from every channel mapped
# to a protein in virtues/{panel}.csv. Quantified in quantify_cells.smk

# variables defined in other .smk files:
# os:                    imported in Snakefile
# re:                    imported in Snakefile
# OUT:                   defined in Snakefile
# CONFIG_PATH:           defined in Snakefile
# SAMPLES:               defined in Snakefile
# SAMPLE_IDS:            defined in Snakefile
# source_image:          defined in Snakefile
# panel_meta:            defined in Snakefile
# VIRTUES_GPU_PARTITION: defined in Snakefile
# GPU_EXTRA:             defined in Snakefile

VIRTUES = config["segmentation"]["virtues"]
virtues_dir = OUT + "/segmentation/virtues/{sample_id}/"
virtues_panel_dir = OUT + "/segmentation/virtues/panels/{panel}/"


def virtues_markers(panel):
    """Marker table of a panel (channel_name, uniprot_id): virtues/{panel}.csv next to the dataset config."""
    return os.path.join(os.path.dirname(CONFIG_PATH), "virtues", panel + ".csv")


# every panel needs a marker table: VirTues identifies channels by their protein
VIRTUES_PANELS = sorted({SAMPLES[s]["panel"] for s in SAMPLE_IDS})
_missing_markers = [virtues_markers(p) for p in VIRTUES_PANELS if not os.path.exists(virtues_markers(p))]
if _missing_markers:
    raise WorkflowError(f"VirTues needs a marker table per panel, missing: {_missing_markers}")

wildcard_constraints:
    panel="|".join(re.escape(p) for p in VIRTUES_PANELS),

# ESM-2 embedding of every protein of the panel (sequences from UniProt)
rule virtues_marker_embeddings:
    input:
        script=f"{workflow.basedir}/scripts/_segmentation/virtues_marker_embeddings.py",
        markers=lambda wc: virtues_markers(wc.panel),
    output:
        embeddings=virtues_panel_dir + "marker_embeddings.pt",
    params:
        model_dir=VIRTUES["model_dir"],
    threads: 4
    resources:
        runtime=60,
        mem='32G',
    log:
        virtues_panel_dir + "logs/virtues_marker_embeddings.log"
    container:
        config["containers"]["pixi"]
    shell:
        """
        {PX} -e gpu python -u {input.script} \
        --markers {input.markers} \
        --model_dir {params.model_dir} \
        --threads {threads} \
        --out_pt {output.embeddings} \
        &> {log}
        """

# mapped channels at VirTues' resolution, clipped and log1p-transformed; per-channel sums for the panel's mean / std
rule virtues_prepare_image:
    input:
        script=f"{workflow.basedir}/scripts/_segmentation/virtues_prepare_image.py",
        # raw or background-subtracted image (source_image, Snakefile)
        image=source_image,
        markers=lambda wc: virtues_markers(SAMPLES[wc.sample_id]["panel"]),
    output:
        image=virtues_dir + "image_1um.npy",
        stats=virtues_dir + "channel_stats.csv",
    params:
        pixel_size=lambda wc: panel_meta(wc)["pixel_size_um"],
        resolution_um=VIRTUES["resolution_um"],
        upper_quantile=VIRTUES["upper_quantile"],
    threads: 1
    resources:
        runtime=120,
        mem='64G',
    log:
        virtues_dir + "logs/virtues_prepare_image.log"
    container:
        config["containers"]["pixi"]
    shell:
        """
        {PX} python -u {input.script} \
        --image {input.image} \
        --markers {input.markers} \
        --pixel_size_fallback {params.pixel_size} \
        --resolution_um {params.resolution_um} \
        --upper_quantile {params.upper_quantile} \
        --threads {threads} \
        --out_image {output.image} \
        --out_stats {output.stats} \
        &> {log}
        """

# mean and std of every channel over the tissue of all samples of the panel (spora's global statistics)
rule virtues_standardization_stats:
    input:
        script=f"{workflow.basedir}/scripts/_segmentation/virtues_standardization_stats.py",
        # from virtues_prepare_image
        stats=lambda wc: expand(OUT + "/segmentation/virtues/{sample_id}/channel_stats.csv",
                                sample_id=[s for s in SAMPLE_IDS if SAMPLES[s]["panel"] == wc.panel]),
    output:
        stats=virtues_panel_dir + "standardization.csv",
    threads: 1
    resources:
        runtime=10,
        mem='4G',
    log:
        virtues_panel_dir + "logs/virtues_standardization_stats.log"
    container:
        config["containers"]["pixi"]
    shell:
        """
        {PX} python -u {input.script} \
        --channel_stats {input.stats} \
        --out_csv {output.stats} \
        &> {log}
        """

# VirTues + segmentation head: whole-cell mask on the image's pixel grid
rule segmentation_virtues:
    input:
        script=f"{workflow.basedir}/scripts/_segmentation/virtues_segmentation.py",
        # from virtues_prepare_image
        image_1um=virtues_dir + "image_1um.npy",
        channel_stats=virtues_dir + "channel_stats.csv",
        # from virtues_standardization_stats and virtues_marker_embeddings (the sample's panel)
        standardization=lambda wc: OUT + f"/segmentation/virtues/panels/{SAMPLES[wc.sample_id]['panel']}/standardization.csv",
        embeddings=lambda wc: OUT + f"/segmentation/virtues/panels/{SAMPLES[wc.sample_id]['panel']}/marker_embeddings.pt",
        # raw or background-subtracted image (source_image, Snakefile): its size only
        image=source_image,
    output:
        cells=virtues_dir + "cell_boundaries_mask.tif",
    params:
        model_dir=VIRTUES["model_dir"],
        model=VIRTUES["model"],
        blur_sigma=VIRTUES["blur_sigma"],
    threads: 1
    resources:
        runtime=240,
        mem='250G',
        partition=VIRTUES_GPU_PARTITION,
        slurm_extra=GPU_EXTRA,
    log:
        virtues_dir + "logs/segmentation_virtues.log"
    container:
        config["containers"]["pixi"]
    shell:
        """
        {PX} -e gpu python -u {input.script} \
        --image_1um {input.image_1um} \
        --channel_stats {input.channel_stats} \
        --standardization {input.standardization} \
        --embeddings {input.embeddings} \
        --image {input.image} \
        --model_dir {params.model_dir} \
        --model {params.model} \
        --blur_sigma {params.blur_sigma} \
        --out_cells {output.cells} \
        &> {log}
        """
