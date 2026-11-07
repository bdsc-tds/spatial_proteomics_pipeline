# Configuration

One folder per dataset, holding a `config.yml`, its `experiments.yml` and the per-panel tables:

```
<dataset folder>/
  config.yml                  output_path, container, run: switches, method settings
  experiments.yml             samples (tissue -> panel -> donor -> [samples]) and per-panel metadata
  signatures/{panel}.csv      TACIT cell type signatures (cell_type, then one column per marker)
  virtues/{panel}.csv         channel -> UniProt ID for VirTues (channel_name, uniprot_id)
  kronos2/{panel}.csv         channel -> KRONOS2 marker (channel_name, kronos2_marker)
  channel_names/{panel}.csv   COMET only: raw channel name -> clean name after background subtraction
```

Every path in the folder (`experiments:` in `config.yml` and the panel tables) is read relative to
`config.yml`, so the folder can live anywhere, e.g. in a private repository next to your analysis code.

`example/` is a template with the tables of five panels: CODEX `human_IO60_metabolism` and `indepth_45plex`,
COMET `spyre_io_core_14plx`, and Xenium protein `hImmune_v1_Colon_100g` and `hImmune_v1_Multi_100g`. The
signatures of the three CODEX/COMET panels only list the markers and still need cell type rows.

## Adding a dataset

1. Copy the template: `cp -r config/example <dataset folder>`.
2. In `config.yml`, set `output_path`, the container and model paths, and the `run:` switches.
3. Remove the copied `experiments.yml` and generate a new one from the raw data:
   ```bash
   pixi run update-experiments --config <dataset folder>/config.yml \
       --base_path <raw> [--xenium_base_path <xenium raw>] --dry-run
   ```
   then run it again without `--dry-run` and fill in the `_` metadata of each new panel.
4. Add the panel tables of your panels (keep or adapt the ones of `example/`).

## Choosing the dataset

Pass the config of the dataset to run:

```bash
./run.sh --configfile <dataset folder>/config.yml
```

Without `--configfile`, the profile reads `config/config.yml`. Everything in `config/` except this README
and `example/` is ignored by git, so dataset folders (or symlinks to them) can be kept there, with
`config/config.yml` as a symlink to the active one:

```bash
ln -s /path/to/my_dataset config/my_dataset
ln -sfn my_dataset/config.yml config/config.yml
```

See the main [README](../README.md) for the meaning of every key.
