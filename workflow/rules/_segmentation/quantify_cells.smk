# Per-cell quantification and GeoJSON export of every segmentation method: quantification is one rule on
# {seg_method}; GeoJSON is one rule per method (generated below), as the number of mask compartments differs

# variables defined in other .smk files:
# source_image:    defined in Snakefile
# OUT:             defined in Snakefile
# panel_meta:      defined in Snakefile
# SEG_METHODS:     defined in Snakefile
# SEG_METHOD_INFO: defined in Snakefile
# seg_targets:     defined in Snakefile
# xenium_qc_mask:  defined in Snakefile

# area, centroid and mean intensity of every channel per cell (every pixel of the cell, the same for every method)
rule quantify_cells:
    input:
        script=f"{workflow.basedir}/scripts/_segmentation/quantify_cells.py",
        # raw or background-subtracted image (source_image, Snakefile)
        image=source_image,
        # from the method's segmentation rule (segmentation_{seg_method}, import_xenium_mask)
        mask=OUT + "/segmentation/{seg_method}/{sample_id}/cell_boundaries_mask.tif",
        # Xenium only: saturation QC masks of the Xenium run (xenium_qc_mask)
        qc_mask=xenium_qc_mask,
    output:
        feature_table=OUT + "/segmentation/{seg_method}/{sample_id}/feature_table.csv",
    params:
        pixel_size=lambda wc: panel_meta(wc)["pixel_size_um"],
        # Xenium only: QC-flagged pixels are left out of the mean
        qc_mask=lambda wc: f"--qc_mask {xenium_qc_mask(wc)}" if xenium_qc_mask(wc) else "",
    threads: 1
    resources:
        # one channel plane and one mask row block in memory at a time
        runtime=480,
        mem='32G',
    log:
        OUT + "/segmentation/{seg_method}/{sample_id}/logs/quantify_cells.log"
    container:
        config["containers"]["pixi"]
    shell:
        """
        {PX} python {input.script} \
        --image {input.image} \
        --mask {input.mask} \
        --pixel_size_fallback {params.pixel_size} \
        --out_csv {output.feature_table} \
        {params.qc_mask} \
        --threads {threads} \
        &> {log}
        """

# shells of the GeoJSON rules below (indented as in the former per-method rules, so the code of existing outputs
# is unchanged): three compartments in one job, or the cell alone
GEOJSON_SHELL_3 = """
        ( {PX} python {input.script} {input.cells} {output.cells} && \
        {PX} python {input.script} {input.nuclei} {output.nuclei} && \
        {PX} python {input.script} {input.cytoplasm} {output.cytoplasm} ) \
        &> {log}
        """
GEOJSON_SHELL_CELL = """
        {PX} python {input.script} {input.cells} {output.cells} &> {log}
        """

# mask compartments as GeoJSON polygons (e.g. for QuPath): cell, nuclei and cytoplasm in one job, or the cell
# alone for methods predicting whole cells only (VirTues)
for _m in SEG_METHODS:
    _seg_dir = OUT + f"/segmentation/{_m}/{{sample_id}}/"
    if SEG_METHOD_INFO[_m]["compartments"] == ["cell", "nuclei", "cytoplasm"]:
        rule:
            name: f"mask_to_geojson_{_m}"
            input:
                script=f"{workflow.basedir}/scripts/_segmentation/mask_to_geojson.py",
                cells=_seg_dir + "cell_boundaries_mask.tif",
                nuclei=_seg_dir + "nuclei_boundaries_mask.tif",
                cytoplasm=_seg_dir + "cytoplasm_boundaries_mask.tif",
            output:
                cells=_seg_dir + "cell_boundaries.geojson",
                nuclei=_seg_dir + "nuclei_boundaries.geojson",
                cytoplasm=_seg_dir + "cytoplasm_boundaries.geojson",
            resources:
                runtime=120,
                mem='60G',
            log:
                _seg_dir + "logs/mask_to_geojson.log"
            container:
                config["containers"]["pixi"]
            shell:
                GEOJSON_SHELL_3
    else:
        rule:
            name: f"mask_to_geojson_{_m}"
            input:
                script=f"{workflow.basedir}/scripts/_segmentation/mask_to_geojson.py",
                cells=_seg_dir + "cell_boundaries_mask.tif",
            output:
                cells=_seg_dir + "cell_boundaries.geojson",
            resources:
                runtime=120,
                mem='60G',
            log:
                _seg_dir + "logs/mask_to_geojson.log"
            container:
                config["containers"]["pixi"]
            shell:
                GEOJSON_SHELL_CELL

# outputs of every segmentation method, on its own samples (nuclei / cytoplasm only where the method has them)
rule segmentation_all:
    input:
        seg_targets(OUT + "/segmentation/{seg_method}/{sample_id}/feature_table.csv"),
        seg_targets(OUT + "/segmentation/{seg_method}/{sample_id}/cell_boundaries.geojson"),
        [t for m in SEG_METHODS for c in SEG_METHOD_INFO[m]["compartments"] if c != "cell"
         for t in seg_targets(OUT + "/segmentation/{seg_method}/{sample_id}/{compartment}_boundaries.geojson",
                              methods=[m], compartment=c)],
