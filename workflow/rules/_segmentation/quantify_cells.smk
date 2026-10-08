# Per-cell quantification and GeoJSON export, one rule per segmentation method: the mesmer rules are
# the base, other methods derive from them with `use rule` (all inputs restated)

# area, centroid and mean intensity of every channel per cell
rule quantify_cells_mesmer:
    input:
        script=f"{workflow.basedir}/scripts/_segmentation/quantify_cells.py",
        # raw or background-subtracted image (source_image, Snakefile)
        image=source_image,
        mask=OUT + "/segmentation/mesmer/{sample_id}/cell_boundaries_mask.npy",
    output:
        feature_table=OUT + "/segmentation/mesmer/{sample_id}/feature_table.csv",
    params:
        pixel_size=lambda wc: panel_meta(wc)["pixel_size_um"],
    threads: 1
    resources:
        # one channel plane and one mask row block in memory at a time
        runtime=480,
        mem='32G',
    log:
        OUT + "/segmentation/mesmer/{sample_id}/logs/quantify_cells.log"
    container:
        config["containers"]["pixi"]
    shell:
        """
        {PX} python {input.script} \
        --image {input.image} \
        --mask {input.mask} \
        --pixel_size_fallback {params.pixel_size} \
        --out_csv {output.feature_table} \
        --threads {threads} \
        &> {log}
        """

use rule quantify_cells_mesmer as quantify_cells_cellpose with:
    input:
        script=f"{workflow.basedir}/scripts/_segmentation/quantify_cells.py",
        # raw or background-subtracted image (source_image, Snakefile)
        image=source_image,
        mask=OUT + "/segmentation/cellpose/{sample_id}/cell_boundaries_mask.npy",
    output:
        feature_table=OUT + "/segmentation/cellpose/{sample_id}/feature_table.csv",
    log:
        OUT + "/segmentation/cellpose/{sample_id}/logs/quantify_cells.log"

# cell, nuclei and cytoplasm masks as GeoJSON polygons (e.g. for QuPath)
rule mask_to_geojson_mesmer:
    input:
        script=f"{workflow.basedir}/scripts/_segmentation/mask_to_geojson.py",
        cells=OUT + "/segmentation/mesmer/{sample_id}/cell_boundaries_mask.npy",
        nuclei=OUT + "/segmentation/mesmer/{sample_id}/nuclei_boundaries_mask.npy",
        cytoplasm=OUT + "/segmentation/mesmer/{sample_id}/cytoplasm_boundaries_mask.npy",
    output:
        cells=OUT + "/segmentation/mesmer/{sample_id}/cell_boundaries.geojson",
        nuclei=OUT + "/segmentation/mesmer/{sample_id}/nuclei_boundaries.geojson",
        cytoplasm=OUT + "/segmentation/mesmer/{sample_id}/cytoplasm_boundaries.geojson",
    resources:
        runtime=480,
        mem='64G',
    log:
        OUT + "/segmentation/mesmer/{sample_id}/logs/mask_to_geojson.log"
    container:
        config["containers"]["pixi"]
    shell:
        """
        ( {PX} python {input.script} {input.cells} {output.cells} && \
        {PX} python {input.script} {input.nuclei} {output.nuclei} && \
        {PX} python {input.script} {input.cytoplasm} {output.cytoplasm} ) \
        &> {log}
        """

use rule mask_to_geojson_mesmer as mask_to_geojson_cellpose with:
    input:
        script=f"{workflow.basedir}/scripts/_segmentation/mask_to_geojson.py",
        cells=OUT + "/segmentation/cellpose/{sample_id}/cell_boundaries_mask.npy",
        nuclei=OUT + "/segmentation/cellpose/{sample_id}/nuclei_boundaries_mask.npy",
        cytoplasm=OUT + "/segmentation/cellpose/{sample_id}/cytoplasm_boundaries_mask.npy",
    output:
        cells=OUT + "/segmentation/cellpose/{sample_id}/cell_boundaries.geojson",
        nuclei=OUT + "/segmentation/cellpose/{sample_id}/nuclei_boundaries.geojson",
        cytoplasm=OUT + "/segmentation/cellpose/{sample_id}/cytoplasm_boundaries.geojson",
    log:
        OUT + "/segmentation/cellpose/{sample_id}/logs/mask_to_geojson.log"

# outputs of every method in segmentation: methods
rule segmentation_all:
    input:
        expand(OUT + "/segmentation/{seg_method}/{sample_id}/feature_table.csv",
               seg_method=SEG_METHODS, sample_id=SAMPLE_IDS),
        expand(OUT + "/segmentation/{seg_method}/{sample_id}/{compartment}_boundaries.geojson",
               seg_method=SEG_METHODS, sample_id=SAMPLE_IDS, compartment=["cell", "nuclei", "cytoplasm"]),
