# Configuration

One folder per dataset, each with two files:

```
config/
  configs_codex/
    config.yml         output_path, container, run: switches, method settings
    experiments.yml    samples (tissue -> panel -> donor -> [samples]) and per-panel metadata
  config.yml -> configs_codex/config.yml     the active dataset (symlink)
```

`experiments:` in `config.yml` points to its `experiments.yml`, relative to the `config.yml` file.

## Choosing the dataset

The SLURM profile reads `config/config.yml`. Point the symlink to the dataset to run:

```bash
ln -sfn configs_codex/config.yml config/config.yml
```

or leave the symlink and pass another config for one run:

```bash
./run.sh --configfile config/configs_comet/config.yml
```

## Adding a dataset

1. Copy an existing folder: `cp -r config/configs_codex config/configs_<name>`.
2. In `config.yml`, set `output_path` (the `processed/` folder next to the raw data) and the `run:`
   switches.
3. Remove the copied `experiments.yml` and generate a new one from the raw data:
   ```bash
   pixi run update-experiments --config config/configs_<name>/config.yml \
       --base_path <raw> [--xenium_base_path <xenium raw>] --dry-run
   ```
   then run it again without `--dry-run` and fill in the `_` metadata of each new panel.

See the main [README](../README.md) for the meaning of every key.
