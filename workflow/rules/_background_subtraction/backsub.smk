# Subtract the autofluorescence background channels from the raw image with backsub (COMET)

# variables defined in other .smk files:
# OUT:                        defined in Snakefile
# CONFIG_PATH:                defined in Snakefile
# SAMPLES:                    defined in Snakefile
# BACKSUB_IMAGE:              defined in Snakefile
# SAMPLE_IDS:                 defined in Snakefile
# use_background_subtraction: defined in Snakefile

background_subtraction_dir = OUT + "/background_subtraction/{sample_id}/"

rule background_subtraction:
    input:
        script=f"{workflow.basedir}/scripts/_background_subtraction/run_backsub.py",
        image=lambda wc: SAMPLES[wc.sample_id]["image"],
        # raw -> clean channel names of the panel: channel_names/{panel}.csv next to the dataset config
        channel_names=lambda wc: os.path.join(os.path.dirname(CONFIG_PATH), "channel_names",
                                              SAMPLES[wc.sample_id]["panel"] + ".csv"),
    output:
        image=BACKSUB_IMAGE,
        markers=background_subtraction_dir + "markers.csv",
        markerout=background_subtraction_dir + "markerout.csv",
    threads: 1
    resources:
        runtime=90,
        mem='24G',
    log:
        background_subtraction_dir + "logs/backsub.log"
    container:
        config["containers"]["pixi"]
    shell:
        """
        {PX} python {input.script} \
        --img_path {input.image} \
        --output_dir $(dirname {output.image}) \
        --channel_names {input.channel_names:q} \
        &> {log}
        """

rule background_subtraction_all:
    input:
        expand(
            rules.background_subtraction.output.image,
            sample_id=[s for s in SAMPLE_IDS if use_background_subtraction(s)],
        ),
