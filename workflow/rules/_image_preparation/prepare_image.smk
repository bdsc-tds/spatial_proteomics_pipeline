# rule to read one sample's multiplexed image (qptiff / OME-TIFF, raw or background-subtracted) and
# write the technology-agnostic inputs of the downstream steps

image_preparation_dir = OUT + "/image_preparation/{sample_id}/"

rule prepare_image:
    input:
        script=f"{workflow.basedir}/scripts/_image_preparation/prepare_image.py",
        image=source_image,
    output:
        channels=image_preparation_dir + "channels.csv",
        dapi=image_preparation_dir + "dapi.ome.tif",
        # only read by segmentation, deleted once every segmentation job has used it
        segmentation_input=temp(image_preparation_dir + "segmentation_input.ome.tif"),
    params:
        nuclear_channel=lambda wc: shlex.quote(panel_meta(wc)["nuclear_channel"]),
        membrane_channels=lambda wc: shell_list(panel_meta(wc)["membrane_channels"]),
        pixel_size=lambda wc: panel_meta(wc)["pixel_size_um"],
    # threads decode the compressed image tiles
    threads: 4
    resources:
        # one full-resolution channel plane (up to ~24k x 24k uint16) in memory at a time
        runtime=240,
        mem='32G',
    log:
        image_preparation_dir + "logs/prepare_image.log"
    container:
        config["containers"]["pixi"]
    shell:
        """
        {PX} python {input.script} \
        --image {input.image} \
        --nuclear_channel {params.nuclear_channel} \
        --membrane_channels {params.membrane_channels} \
        --pixel_size_fallback {params.pixel_size} \
        --out_channels {output.channels} \
        --out_dapi {output.dapi} \
        --out_segmentation_input {output.segmentation_input} \
        --threads {threads} \
        &> {log}
        """

rule image_preparation_all:
    input:
        # the temporary segmentation input is left out so that it is not kept as a final target
        expand([rules.prepare_image.output.channels, rules.prepare_image.output.dapi], sample_id=SAMPLE_IDS),
