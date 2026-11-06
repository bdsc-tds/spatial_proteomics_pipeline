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

The nuclear and membrane signals given to the segmentation methods are each the per-pixel maximum of their channels
(raw values), after removing every channel's pixel value offset read from the OME metadata (+100 on
Xenium protein channels, 0 for DAPI, stains and formats without one). Values below the offset are kept negative
(float32 segmentation input), as in the per-cell mean intensities and the KRONOS2 embeddings; only VirTues clips
them at 0.

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

Each sample's image and Xenium run are not stored: they are found when the pipeline starts. Xenium
runs are expected as `{_xenium_base_path}/{tissue}/{panel}/{donor}/{sample}/` (the
xenium_analysis_pipeline layout, `sample` being the Xenium output folder). Names are
compared case-insensitively with `-` and `_` treated alike, and one name "contains" another when it
holds it as whole `_`-separated tokens. A proteomics sample is matched to the run whose donor contains
the proteomics donor and whose donor or sample contains the proteomics sample, e.g.
`donor_01/donor_01` →
`donor_01/output-XETG00000__slide_01__donor_01__…`. Tissue and panel
names are ignored. Samples without a match, or with several, get no alignment. The pipeline also warns about raw samples missing
from experiments.yml.

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
| `image_preparation` | `prepare_dapi`, `prepare_segmentation_input` | Reads the raw image (qptiff or OME-TIFF) one channel at a time; `prepare_dapi` writes the channel table and a DAPI image (alignment input), `prepare_segmentation_input` a two-channel segmentation input (nuclear and membrane max projections, offsets removed). Split so that changing segmentation markers does not re-run the alignment |
| `background_subtraction` | `background_subtraction` | COMET only: subtracts the autofluorescence channels with [backsub](https://pypi.org/project/backsub/). When on, all later steps read the subtracted image |
| `segmentation` | `segmentation_{method}`, `quantify_cells_{method}`, `mask_to_geojson_{method}` | Nucleus / whole-cell / cytoplasm masks for each method in `segmentation: methods`: Mesmer (deepcell, through SPACEc) and/or Cellpose 3 (through SPACEc: `nuclei` model on the nuclear signal, `cyto3` on membrane + nuclear, diameters set in µm) and/or VirTues (see below; whole cells only); per-cell area, centroid and mean intensity of **every** channel; masks as GeoJSON (e.g. for QuPath) |
| `xenium_alignment` | `align_to_xenium_insitupy`, `plot_alignment` | Affine registration of the proteomics DAPI onto the Xenium DAPI with [insitupy](https://github.com/SpatialPathology/InSituPy); cells matched to Xenium cells (mutual nearest neighbours), overlay and count plots. Only for samples with a Xenium run |
| `cell_type_annotation` | `tacit_{seg_method}`, `tacit_xenium` | Cell types of every segmentation method's cells for each method in `cell_type_annotation: methods`: TACIT (marker signatures, `signatures/{panel}.csv`). Samples imaged on Xenium also get cell types for 10x's own cells (`cell_type_annotation: xenium_segmentation`) |
| `dimensionality_reduction` | `{dr_method}_embeddings_{seg_method}`, `xenium_cell_mask`, `{dr_method}_embeddings_xenium` | Per-cell embeddings of every segmentation method's cells for each method in `dimensionality_reduction: methods`: VirTues (see below; needs `virtues` in `segmentation: methods`) and/or KRONOS2 (see below). Samples imaged on Xenium also embed 10x's own cells, exported as a label mask first |

**VirTues** ([bunnelab/virtues](https://github.com/bunnelab/virtues), `gpu` environment, A100 GPU): a
foundation model that reads any panel, identifying each channel by the ESM-2 embedding of its protein
(`virtues/{panel}.csv` next to the dataset config: `channel_name, uniprot_id`, blank = not used; DAPI
as Histone H3 P68431). The mapped channels are resampled to 1 µm/px, clipped at their 99th percentile
in the tissue, log1p-transformed, z-scored with the panel's mean / std (spora's standardization) and
blurred (3 x 3 Gaussian). The segmentation head then predicts whole cells; the VirTues cell mask is
upsampled to the image's pixel grid like the other methods'. Cell embeddings follow VirTues' cell
tokens: the backbone embeds 128 px crops (stride 42, image zero-padded by 120 px) into one 512-d token
per 8 x 8 px patch, and each cell gets the mean of the tokens it overlaps, weighted by its number of
1 µm pixels in each patch. Masks are brought to 1 µm by centre sampling (each 1 µm pixel takes the
label of the full-resolution pixel at its centre); cells covering no pixel centre (slivers under ~1 µm²)
get no embedding. The sums are accumulated crop by crop, as VirTues' own `compute_cell_tokens` runs
out of memory on whole slides.

**KRONOS2** ([MahmoodLab/KRONOS2](https://huggingface.co/MahmoodLab/KRONOS2), gated, CC BY-NC-ND 4.0;
`gpu` environment): a marker-aware ViT that identifies each channel by its own marker name
(`kronos2/{panel}.csv` next to the dataset config: `channel_name, kronos2_marker`, blank = not used; the
names are those of `marker_metadata.csv` in the model folder, an unknown name stops the run). Cell
embeddings follow the cell recipe of [CORAL](https://github.com/mahmoodlab/CORAL), at the image's native
resolution: a `box_size` px box (64) around each cell's rounded centroid, pixel value offset removed in float
(values below it stay negative), pixels outside the cell set to 0, divided by the dtype maximum, z-scored with KRONOS2's pretraining statistics; the embedding is the 768-d
CLS token. The rules load the model from a local copy (`dimensionality_reduction: kronos2: model_dir`) and
need no token. Download it once, after requesting access on Hugging Face:
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
