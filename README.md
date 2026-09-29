# Spatial proteomics pipeline

A [Snakemake](https://snakemake.readthedocs.io/) workflow for multiplexed spatial proteomics images
(CODEX / PhenoCycler, COMET / Lunaphore, …). For each sample it prepares the image, segments cells,
quantifies every marker per cell and, when the same slide was also run on 10x Xenium, aligns the
proteomics DAPI to the Xenium DAPI and matches the cells of both modalities.

The pipeline is technology-agnostic: everything technology-specific (image format, channel names,
pixel size, which channels to segment on) is read from the image or set per panel in
`experiments.yml`.

## Steps

| Switch (`run:`) | Rules | What it does |
|---|---|---|
| `image_preparation` | `prepare_image` | Reads the raw image (qptiff or OME-TIFF) one channel at a time; writes the channel table, a DAPI image and a small segmentation input (nuclear + membrane channels) |
| `background_subtraction` | `background_subtraction` | COMET only: subtracts the autofluorescence channels with [backsub](https://pypi.org/project/backsub/). When on, all later steps read the subtracted image |
| `segmentation` | `segmentation_mesmer`, `quantify_cells`, `mask_to_geojson` | Mesmer (deepcell, through SPACEc) nucleus / whole-cell / cytoplasm masks; per-cell area, centroid and mean intensity of **every** channel; masks as GeoJSON (e.g. for QuPath) |
| `xenium_alignment` | `align_to_xenium_insitupy`, `plot_alignment` | Affine registration of the proteomics DAPI onto the Xenium DAPI with [insitupy](https://github.com/SpatialPathology/InSituPy); cells matched to Xenium cells (mutual nearest neighbours), overlay and count plots. Only for samples with a Xenium run |

Switches select what is built; upstream steps still run when a switched-on step needs them.

## Data layout

Raw images, one per sample folder:

```
{raw}/{tissue}/{panel}/{donor}/{sample}/{sample}.{qptiff|ome.tiff|...}
```

e.g. `codex/raw/NSCLC/human_IO60_metabolism/donor_01/donor_01/donor_01.qptiff`.
Sample folders without an image yet (e.g. only a crop geojson) are ignored.

The sample ID used throughout is `tissue/panel/donor/sample`. All results go to a `processed/`
folder next to `raw/` (`output_path` in the config), one folder per analysis:

```
processed/
  image_preparation/{sample_id}/                    channels.csv, dapi.ome.tif, segmentation_input.ome.tif, logs/
  background_subtraction/{sample_id}/               image.ome.tiff, markers.csv, markerout.csv, logs/
  segmentation/mesmer/{sample_id}/                  *_boundaries_mask.npy, feature_table.csv, *_boundaries.geojson, logs/
  xenium_alignment/insitupy/{sample_id}/            affine_matrix.csv, mesmer_cell_matches.csv, figures/, logs/
  snkmt.db                                          run monitoring database
```

## Setup

The pipeline runs on the Curnagl SLURM cluster; every rule runs inside one Apptainer image that holds
the pixi environments:

| pixi environment | Used by |
|---|---|
| `default` | image preparation, background subtraction, quantification, GeoJSON, insitupy alignment, plots, `update-experiments` |
| `spacec` | Mesmer segmentation (deepcell / TensorFlow 2.8 need python 3.10 and numpy < 2) |
| `snakemake` | the Snakemake driver itself (runs outside the container) |

1. Install the environments locally (the `snakemake` one launches the pipeline, `default` runs
   `update-experiments`):
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

One folder per dataset under `config/`, e.g. `config/configs_codex/`:

- **`config.yml`**: `output_path`, the container, the `run:` switches and method settings
  (segmentation / alignment methods, Mesmer model folder, matching distance, GPU partition).
- **`experiments.yml`**: the samples, as tissue → panel → donor → `[samples]`. Keys starting with `_`
  are metadata:

  ```yaml
  _base_path: /path/to/codex/raw                 # raw images
  _xenium_base_path: /path/to/xenium/raw         # optional: Xenium runs of the same slides
  NSCLC:
    human_IO60_metabolism:
      _technology: codex
      _image_ext: qptiff                          # image file extension in each sample folder
      _pixel_size_um: 0.5073566                   # fallback when the image carries no pixel size
      _nuclear_channel: DAPI
      _membrane_channels: [Pan-Cytokeratin, E-cadherin, CD45, Vimentin, CD31, SMA]   # Mesmer whole-cell signal
      _background_subtraction: false              # COMET only
      donor_01: [donor_01]          # donor: [samples]
  ```

`config/config.yml` is a symlink to the active dataset's config (see [config/README.md](config/README.md)).

### Creating / updating experiments.yml

Do not write the sample list by hand: generate it from the raw data folder, before the first run and
whenever new samples arrive:

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

Each sample's image and Xenium run are not stored: they are found when the pipeline starts. The
Xenium run is the folder below `_xenium_base_path` (holding `experiment.xenium`) whose name, or the
sample field of its `output-<instrument>__<slide>__<sample>__…` name, equals the sample name exactly.
Samples without a match simply get no alignment. The pipeline also warns about raw samples missing
from experiments.yml.

## Running

```bash
./run.sh -n                                              # dry run: list the jobs
./run.sh                                                 # run with the active config
./run.sh --configfile config/configs_comet/config.yml    # another dataset
./run.sh <rule or output file>                           # only that target
```

`run.sh` loads `singularityce`, starts Snakemake from the `snakemake` pixi environment with the SLURM
profile (`profiles/slurm/`, one job per rule, cluster-generic executor) and `--sdm apptainer`, and
stores the monitoring database in `{output_path}/snkmt.db`. The Snakemake driver keeps running while
jobs are queued, so start it in a persistent session (e.g. `tmux`/`screen`). Each rule's messages are
in the `logs/` folder next to its outputs.

Mesmer runs on the `gpu` partition (set `gpu: available: false` in the config to use CPUs). Its
log states whether TensorFlow sees the GPU.

## Extending

- **New panel / technology:** put the images in `{raw}/{tissue}/{panel}/{donor}/{sample}/`, run
  `update-experiments` and fill in the new panel's `_` keys. qptiff and OME-TIFF are read by
  `workflow/scripts/utils/readwrite.py`; another format only needs support there.
- **New dataset:** copy `config/configs_codex/` to `config/configs_<name>/` and adapt `output_path`,
  `_base_path` and the switches.
- **New segmentation method:** add `workflow/rules/_segmentation/<method>.smk` with a rule writing
  `{output_path}/segmentation/<method>/{sample_id}/{cell,nuclei,cytoplasm}_boundaries_mask.npy`, include
  it in `workflow/rules/segmentation.smk` and add the method to `segmentation: methods`;
  quantification, GeoJSON export and the alignment plots then work for it unchanged.

## Repository layout

```
run.sh                                   launch the pipeline on SLURM
pixi.toml, pixi.lock                     environments (default, spacec, snakemake)
reproducibility/                         container definition and build script
profiles/slurm/config.v8+.yaml           Snakemake SLURM profile
config/configs_<dataset>/                config.yml + experiments.yml per dataset
workflow/
  Snakefile                              samples, helpers, run switches, rule all
  rules/<topic>.smk                      includes + run switch -> rules map
  rules/_<topic>/<method>.smk            the rules, ending with the step's _all rule
  scripts/_<topic>/                      one argparse script per rule
  scripts/utils/                         experiments.yml handling, image reading, alignment helpers
```

## Known limitations

- Xenium runs are matched to samples by exact sample name; data with other naming may need an
  explicit mapping.
- The COMET dataset config (`config/configs_comet/`) is not written yet.

## License

This project is licensed under the terms of the [MIT License](LICENCE).

## Citation

If you use this pipeline in your research, please cite it using the metadata in [CITATION.cff](CITATION.cff).
