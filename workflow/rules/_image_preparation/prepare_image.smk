# Channel table, DAPI image and segmentation input of each sample (from qptiff / OME-TIFF).
# Two rules so that changing the segmentation markers does not re-run the DAPI (and the Xenium alignment)

image_preparation_dir = OUT + "/image_preparation/{sample_id}/"

# channel table + first nuclear channel alone, for the alignment to Xenium
rule prepare_dapi:
    input:
        script=f"{workflow.basedir}/scripts/_image_preparation/prepare_dapi.py",
        image=source_image,
    output:
        channels=image_preparation_dir + "channels.csv",
        dapi=image_preparation_dir + "dapi.ome.tif",
    params:
        # only the first nuclear channel, so other marker changes leave this rule untouched
        dapi_channel=lambda wc: shell_list(panel_meta(wc)["nuclear_channels"][:1]),
        pixel_size=lambda wc: panel_meta(wc)["pixel_size_um"],
    # threads decode the compressed image tiles
    threads: 1
    resources:
        # one full-resolution channel plane (up to ~24k x 24k uint16) in memory
        runtime=240,
        mem='32G',
    log:
        image_preparation_dir + "logs/prepare_dapi.log"
    container:
        config["containers"]["pixi"]
    shell:
        """
        {PX} python {input.script} \
        --image {input.image} \
        --dapi_channel {params.dapi_channel} \
        --pixel_size_fallback {params.pixel_size} \
        --out_channels {output.channels} \
        --out_dapi {output.dapi} \
        --threads {threads} \
        &> {log}
        """

# nuclear and membrane max projections, input of the segmentation methods
rule prepare_segmentation_input:
    input:
        script=f"{workflow.basedir}/scripts/_image_preparation/prepare_segmentation_input.py",
        image=source_image,
    output:
        # deleted once all segmentation jobs have read it
        segmentation_input=temp(image_preparation_dir + "segmentation_input.ome.tif"),
    params:
        nuclear_channels=lambda wc: shell_list(panel_meta(wc)["nuclear_channels"]),
        membrane_channels=lambda wc: shell_list(panel_meta(wc)["membrane_channels"]),
        pixel_size=lambda wc: panel_meta(wc)["pixel_size_um"],
    # threads decode the compressed image tiles
    threads: 1
    resources:
        # one full-resolution channel plane (up to ~24k x 24k uint16) in memory at a time
        runtime=240,
        mem='32G',
    log:
        image_preparation_dir + "logs/prepare_segmentation_input.log"
    container:
        config["containers"]["pixi"]
    shell:
        """
        {PX} python {input.script} \
        --image {input.image} \
        --nuclear_channels {params.nuclear_channels} \
        --membrane_channels {params.membrane_channels} \
        --pixel_size_fallback {params.pixel_size} \
        --out_segmentation_input {output.segmentation_input} \
        --threads {threads} \
        &> {log}
        """

rule image_preparation_all:
    input:
        # temp segmentation input left out on purpose
        expand(rules.prepare_dapi.output, sample_id=SAMPLE_IDS),
