# Lunaphore–Xenium Image Processing Pipeline

A [Snakemake](https://snakemake.readthedocs.io/) workflow for processing COMET Lunaphore multiplexed protein imaging data (post-Xenium samples), aligning it to matched 10x Xenium spatial transcriptomics data, and performing cell segmentation with downstream conversion to GeoJSON for use in tools like QuPath.

## Overview

The pipeline takes raw COMET Lunaphore spatial proteomics images (`.ome.tiff`), and for each sample:

1. **Background subtraction** — removes background signal from the raw image and extracts marker information.
2. **Image cropping** — crops the image to a region of interest defined by user-supplied coordinates.
3. **Image alignment** — computes an affine transformation matrix to register the cropped image with the corresponding Xenium output.
4. **Cell segmentation** — segments cells, nuclei, and cytoplasm from the background-subtracted image.
5. **Plotting** — generates aligned-image plots and cell-count barplots comparing COMET and Xenium data.
6. **Mask → GeoJSON conversion** — converts the cell/nuclei/cytoplasm segmentation masks (`.npy`) into GeoJSON files.

## Directory structure expected

Raw images must follow this layout under `main_dir`:

```
{main_dir}/lunaphore_protein/raw/{data_type}/uncropped/{hid}/{rid}/{sample_name}.ome.tiff
```

Where:
- `data_type` — sample category
- `hid` — Id for the sample according to Xenium 
- `rid` — Id for the sample according to COMET

## Requirements

- [Snakemake](https://snakemake.readthedocs.io/) = 9.20.0 (with `conda`/`mamba` support)
- Conda/mamba for environment management
- Python 3, plus the following per-rule conda environments (provided under `envs/`):
  - `envs/backsub.yaml`
  - `envs/insitupy.yaml`
  - `envs/spacec.yaml`
  - `envs/squid_lazy_expl.yaml`

### Activate the snakemake environment
The environments above are used as individual rules, but you still need snakemake itself installed to launch the workflow in the first place. A `sp-pipeline.yaml` is provided in the `envs` folder exactly for this purpose:

```bash
conda env create -f envs/sp-pipeline.yaml
conda activate sp-pipeline
```

This install the right version of snakemake to launch the program. Once activated, you can run the pipeline as described below.

## Configuration

All settings are provided via `config.yaml`. Required keys:

| Key | Description |
|---|---|
| `main_dir` | Root directory containing the `lunaphore_protein/` raw data folder |
| `align_methods` | List of alignment method names to run/expand over |
| `segmentation_methods` | List of segmentation method names to run/expand over |
| `paths` | Dictionary of path templates (relative to `main_dir`) used by `get_path()`, keyed by the names below |

### Path keys used by `paths`

| Path key | Produced/used by |
|---|---|
| `img` | Input raw image glob pattern |
| `backsub_img`, `backsub_markers`, `backsub_markerout`, `backsub_log` | `background_substraction` |
| `cropped_coordinates`, `cropped_img`, `cropped_geojson`, `cropped_log` | `crop_image` |
| `xenium_outs` | Xenium output directory (input to alignment) |
| `xenium_adata` | Xenium AnnData object (input to plotting) |
| `align_matrix`, `align_log` | `align_image` |
| `plot_aligned_samples`, `plot_aligned_bar`, `plot_aligned_log` | `plot_alignment_and_data` |
| `seg_csv`, `seg_npy_cells`, `seg_npy_nuclei`, `seg_npy_cytoplasm`, `seg_txt`, `seg_log` | `segmentation_mesmer` |
| `conversion_geojson_cells_file`, `conversion_geojson_nuclei_file`, `conversion_geojson_cytoplasm_file`, `conversion_log` | `npy_to_geojson` |

Each path template is expected to contain the wildcards used to `.format()` it (e.g. `{data_type}`, `{hid}`, `{rid}`, `{align}`) as appropriate to the rule that consumes it.


## Running the pipeline

```bash
snakemake --profile profile/
```

Useful flags:
- `-n` / `--dry-run` — preview the jobs that would run
- `--until $rule_name` - is going to run only the required rules to get the outputs of $rule_name
- For cluster/HPC execution, combine with a Snakemake [executor profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) (e.g. Slurm) to make use of the dynamic memory resources described below.

## Rules

| Rule | Description | Conda env |
|---|---|---|
| `background_substraction` | Subtract background from raw protein image; output cleaned image + marker CSVs | `envs/backsub.yaml` |
| `crop_image` | Crop image using supplied coordinates CSV | `envs/insitupy.yaml` |
| `align_image` | Align cropped image to Xenium output via affine transform | `envs/insitupy.yaml` |
| `plot_alignment_and_data` | Plot aligned images and cell-count barplots | `envs/squid_lazy_expl.yaml` |
| `segmentation_mesmer` | Segment cells/nuclei/cytoplasm using Mesmer | `envs/spacec.yaml` |
| `npy_to_geojson` | Convert `.npy` segmentation masks to GeoJSON | `envs/squid_lazy_expl.yaml` |

## Retry / failure handling

- All compute-heavy rules use `get_mem_mb_exp` as their memory resource function, which scales memory on retry (via Snakemake's `--retries` specified in the profile file), so transient out-of-memory failures are automatically retried with more RAM.
- The `align_image` rule additionally tracks failures per sample: if alignment fails on the 4th attempt, an `error.txt` marker file is written to the sample's output directory and the sample is skipped on subsequent runs (a message is printed instead of re-attempting) until that marker is removed.

## Outputs

For each processed sample (`data_type`, `hid`, `rid`), and depending on the configured `align_methods` / `segmentation_methods`, the pipeline produces:

- Background-subtracted image and marker tables
- Cropped image and ROI GeoJSON
- Affine alignment matrix (per alignment method)
- Aligned-image plots and cell-count barplots (per alignment × segmentation method)
- Segmentation outputs: cell CSV, and `.npy` masks for cells/nuclei/cytoplasm, plus a channel list
- GeoJSON files for cells, nuclei, and cytoplasm masks (for viewing in QuPath or similar tools)

## Repository layout

```
.
├── Snakefile
├── config.yaml
├── envs/
│   ├── backsub.yaml
│   ├── insitupy.yaml
│   ├── spacec.yaml
│   └── squid_lazy_expl.yaml
└── scripts/
    ├── background_substraction/
    │   └── background_substraction.py
    ├── alignment/
    │   ├── image_cropping.py
    │   └── create_alignment.py
    ├── segmentation/
    │   └── mesmer_segmentation.py
    └── utils/
        ├── plot_aligned_samples.py
        └── convert_mask.py
        └── utils.py
```

## License

This project is licensed under the terms of the [MIT License](LICENCE).

## Citation

If you use this pipeline in your research, please cite it using the metadata in [CITATION.cff](CITATION.cff)