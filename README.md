# Spatial proteomics pipeline

A [Snakemake](https://snakemake.readthedocs.io/) workflow for multiplexed spatial proteomics images
(CODEX / PhenoCycler, COMET / Lunaphore, …). For each sample it prepares the image, segments cells,
quantifies every marker per cell and, when the same slide was also run on 10x Xenium, aligns the
proteomics DAPI to the Xenium DAPI and matches the cells of both modalities.

The pipeline is technology-agnostic: everything technology-specific (image format, channel names,
pixel size, which channels to segment on) is read from the image or set per panel in
`experiments.yml`.

## Data layout

Raw images, one per sample folder:

```
{raw}/{tissue}/{panel}/{donor}/{sample}/{sample}.{qptiff|ome.tiff|...}
```

e.g. `codex/raw/NSCLC/human_IO60_metabolism/donor_01/donor_01/donor_01.qptiff`.
Sample folders without an image yet (e.g. only a crop geojson) are ignored. The image may also sit
in a subfolder of the sample folder, e.g. Xenium protein data (`config/configs_xenium_mm/`), where
the sample folder is the Xenium output folder and the image is the multi-file OME-TIFF opened
through `morphology_focus/ch0000_dapi.ome.tif` (see `_image_ext` below).

The sample ID used throughout is `tissue/panel/donor/sample`. All results go to a `processed/`
folder next to `raw/` (`output_path` in the config), one folder per analysis:

```
processed/
  image_preparation/{sample_id}/                    channels.csv, dapi.ome.tif, segmentation_input.ome.tif, logs/
  background_subtraction/{sample_id}/               image.ome.tiff, markers.csv, markerout.csv, logs/
  segmentation/{seg_method}/{sample_id}/            *_boundaries_mask.tif, feature_table.csv, *_boundaries.geojson, logs/  (mesmer, cellpose, virtues)
  segmentation/virtues/{sample_id}/                 + image_1um.npy, channel_stats.csv (cell mask only, no nuclei / cytoplasm)
  segmentation/virtues/panels/{panel}/              marker_embeddings.pt, standardization.csv
  xenium_alignment/insitupy/max_width_{max_width}/{sample_id}/  affine_matrix.csv, {seg_method}_cell_matches.csv, figures/, logs/
  cell_type_annotation/{ct_method}/{seg_method}/{sample_id}/  cell_types.csv, logs/  (tacit; also xenium_{seg} for 10x's cells)
  dimensionality_reduction/{dr_method}/{seg_method}/{sample_id}/  cell_embeddings.h5, logs/  (virtues, kronos2; also xenium_{seg} for 10x's cells)
  dimensionality_reduction/xenium_{seg}/{sample_id}/  cell_boundaries_mask.tif, cell_ids.csv (10x's cells as a label mask)
  snkmt.db                                          run monitoring database
```

## Setup

The pipeline is designed to run on HPC clusters and has been tested with SLURM. Every rule runs
inside one Apptainer image that holds the pixi environments.

1. Install the `snakemake` and `default` environments locally:
   ```bash
   pixi install -e snakemake -e default
   ```
2. Build the container from `pixi.toml` / `pixi.lock` (a SLURM job, up to 4 h); the image path must match
   `containers: pixi:` in the config:
   ```bash
   sbatch reproducibility/build_container.sh
   ```
   Changes to `pixi.toml` only reach the rules after the image is rebuilt.

## Configuration

One folder per dataset under `config/`, e.g. `config/configs_codex/`, holding a `config.yml` and an
`experiments.yml`.

### config.yml

`output_path`, the container, the `run:` switches and method settings (segmentation / alignment
methods, Mesmer / Cellpose model folders and Cellpose models and diameters, matching distance, GPU
partition). Only the segmentation methods listed in `segmentation: methods` need a settings block.

`config/config.yml` is a symlink to the active dataset's config (see [config/README.md](config/README.md)).

### experiments.yml

The samples, as tissue → panel → donor → `[samples]`. Keys starting with `_` are metadata:

```yaml
_base_path: /path/to/codex/raw                 # raw images
_xenium_base_path: /path/to/xenium/raw         # optional: Xenium runs of the same slides
NSCLC:
  human_IO60_metabolism:
    _technology: codex
    _image_ext: qptiff                          # image file extension in each sample folder
    _pixel_size_um: 0.5073566                   # fallback when the image carries no pixel size
    _nuclear_channels: [DAPI]                   # nuclear signal for segmentation; the first one is used for alignment
    _membrane_channels: [Pan-Cytokeratin, E-cadherin, CD45, Vimentin, CD31, SMA]   # whole-cell signal for segmentation
    _background_subtraction: false              # COMET only
    donor_01: [donor_01]          # donor: [samples]
```

`_image_ext` is how the image path inside the sample folder ends: an extension (`qptiff`,
`ome.tiff`), optionally behind a subfolder (`images/ome.tiff`), or a full file name
(`morphology_focus/ch0000_dapi.ome.tif` for Xenium). Exactly one file must match.

#### Creating / updating experiments.yml

The file can be written by hand, or generated from the raw data folder with `update-experiments`
(before the first run and whenever new samples arrive):

```bash
pixi run update-experiments --config config/configs_codex/config.yml --dry-run   # show changes
pixi run update-experiments --config config/configs_codex/config.yml             # write them
```

(first time, when experiments.yml does not exist yet, also pass `--base_path <raw>` and optionally
`--xenium_base_path <xenium raw>`). The script:

- adds new tissues / panels / donors / samples and never overwrites what is already there, so
  hand-edited panel metadata and comments are kept;
- gives a new panel a block of `_` keys set to `null` (or a copy of another panel's with
  `--panel_template <panel>`); the pipeline refuses to run until they are filled in;
- previews which Xenium run each sample will be aligned to;
- reports samples whose image disappeared (removed only with `--prune`).

Images and Xenium runs are located when the pipeline starts, with Xenium runs expected as
`{_xenium_base_path}/{tissue}/{panel}/{donor}/{sample}/`. A sample is matched to the run whose donor
contains its donor and whose donor or sample contains its sample (case-insensitive, `-` ≡ `_`, whole
`_`-separated tokens); samples with no match or several matches get no alignment.

## Running

```bash
./run.sh -n                                              # dry run: list the jobs
./run.sh                                                 # run with the active config
./run.sh <rule or output file>                           # only that target
```

`run.sh` loads `singularityce`, starts Snakemake from the `snakemake` pixi environment with the SLURM
profile (`profiles/slurm/`, one job per rule, cluster-generic executor) and `--sdm apptainer`, and
stores the monitoring database in `{output_path}/snkmt.db`. 

Some rules run on the `gpu` partition:  (set `gpu: available: false` in the config to use CPUs). 

## Steps

| Switch (`run:`) | Rules | What it does |
|---|---|---|
| `image_preparation` | `prepare_dapi`, `prepare_segmentation_input` | Writes the channel table, a DAPI image for alignment and a two-channel (nuclear, membrane) segmentation input from the raw image. |
| `background_subtraction` | `background_subtraction` | COMET only: subtracts the autofluorescence channels with [backsub](https://pypi.org/project/backsub/), and all later steps read the subtracted image. |
| `segmentation` | `segmentation_{method}`, `quantify_cells_{method}`, `mask_to_geojson_{method}` | Segments cells with each method in `segmentation: methods` (Mesmer, Cellpose 3, VirTues) and writes masks, per-cell mean intensities of every channel and GeoJSON outlines. |
| `xenium_alignment` | `align_to_xenium_insitupy`, `plot_alignment` | Registers the proteomics DAPI onto the Xenium DAPI with [insitupy](https://github.com/SpatialPathology/InSituPy) and matches cells to Xenium cells (samples with a Xenium run only). |
| `cell_type_annotation` | `tacit_{seg_method}`, `tacit_xenium` | Assigns cell types with TACIT (`signatures/{panel}.csv`) to every segmentation method's cells, and to 10x's own cells for Xenium samples. |
| `dimensionality_reduction` | `{dr_method}_embeddings_{seg_method}`, `xenium_cell_mask`, `{dr_method}_embeddings_xenium` | Computes per-cell embeddings with VirTues and/or KRONOS2 for every segmentation method's cells, and for 10x's own cells on Xenium samples. |

**VirTues** ([bunnelab/virtues](https://github.com/bunnelab/virtues)) needs an A100 GPU and a
`virtues/{panel}.csv` next to the dataset config mapping each channel to its protein (`channel_name,
uniprot_id`, blank = not used; DAPI as Histone H3 P68431). The weights and ESM-2 are downloaded to
`segmentation: virtues: model_dir` on first use.

**KRONOS2** ([MahmoodLab/KRONOS2](https://huggingface.co/MahmoodLab/KRONOS2)) needs a `kronos2/{panel}.csv`
next to the dataset config mapping each channel to a KRONOS2 marker (`channel_name, kronos2_marker`, blank =
not used; names from `marker_metadata.csv` in the model folder). The model is gated: request access on
Hugging Face, then download it once to `dimensionality_reduction: kronos2: model_dir`:
```bash
export CONDA_OVERRIDE_CUDA=12.9
pixi run -e gpu hf auth login    # read token
pixi run -e gpu hf download MahmoodLab/KRONOS2 --exclude "demo_image/*" --local-dir <model_dir>
```

A switched-off step is not loaded: later steps read its existing outputs as plain input files, so it
is never re-run (Snakemake stops with a missing input error if they do not exist yet). Background
subtraction is the exception: its switch also chooses the image every step reads.

## Extending

- **New panel / technology:** put the images in `{raw}/{tissue}/{panel}/{donor}/{sample}/`, run
  `update-experiments` and fill in the new panel's `_` keys. qptiff and OME-TIFF are read by
  `workflow/scripts/utils/readwrite.py`; another format only needs support there.
- **New segmentation method:** add `workflow/rules/_segmentation/<method>.smk` with a
  `segmentation_<method>` rule writing `cell_`, `nuclei_` and `cytoplasm_boundaries_mask.tif` (tiled,
  compressed, all labelled with the cell id, see `scripts/_segmentation/label_masks.py`) to
  `segmentation/<method>/{sample_id}/`, add `quantify_cells_<method>` and `mask_to_geojson_<method>`
  to `_segmentation/quantify_cells.smk` (`use rule quantify_cells_mesmer as quantify_cells_<method>
  with:`, changing the paths), and list the method in `segmentation: methods`.
- **New dataset:** copy `config/configs_codex/` to `config/configs_<name>/` and adapt `output_path`,
  `_base_path` and the switches.

## Repository layout

```
run.sh                                   launch the pipeline on SLURM
pixi.toml, pixi.lock                     environments (default, spacec, tacit, gpu, snakemake)
reproducibility/                         container definition and build script
profiles/slurm/config.v8+.yaml           Snakemake SLURM profile
config/configs_<dataset>/                config.yml + experiments.yml per dataset
workflow/
  Snakefile                              samples, helpers, run switches, rule all
  rules/<topic>.smk                      includes + run switch -> rules map
  rules/_<topic>/<method>.smk            the rules, ending with the step's _all rule
  rules/_segmentation/quantify_cells.smk quantify_cells / mask_to_geojson of every segmentation method
  scripts/_<topic>/                      one argparse script per rule
  scripts/utils/                         experiments.yml handling, image reading, alignment helpers
```

## Known limitations

- Xenium runs are matched to samples by donor and sample name (see above); data whose names do not
  contain each other may need an explicit mapping.
- The COMET dataset config (`config/configs_comet/`) is not written yet.

## License

This project is licensed under the terms of the [MIT License](LICENCE).

## Citation

If you use this pipeline in your research, please cite it using the metadata in [CITATION.cff](CITATION.cff).
