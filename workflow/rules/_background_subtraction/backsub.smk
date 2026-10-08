# Subtract the autofluorescence background channels from the raw image with backsub (COMET)

background_subtraction_dir = OUT + "/background_subtraction/{sample_id}/"

rule background_subtraction:
    input:
        script=f"{workflow.basedir}/scripts/_background_subtraction/run_backsub.py",
        image=lambda wc: SAMPLES[wc.sample_id]["image"],
    output:
        image=BACKSUB_IMAGE,
        markers=background_subtraction_dir + "markers.csv",
        markerout=background_subtraction_dir + "markerout.csv",
    threads: 4
    resources:
        runtime=480,
        mem='64G',
    log:
        background_subtraction_dir + "logs/backsub.log"
    container:
        config["containers"]["pixi"]
    shell:
        """
        {PX} python {input.script} \
        --img_path {input.image} \
        --output_dir $(dirname {output.image}) \
        &> {log}
        """

rule background_subtraction_all:
    input:
        expand(
            rules.background_subtraction.output.image,
            sample_id=[s for s in SAMPLE_IDS if use_background_subtraction(s)],
        ),
