#!/bin/bash -l
# Run the pipeline on SLURM with the active dataset config (config/config.yml), e.g.
#   ./run.sh -n                                             # dry run
#   ./run.sh --configfile config/configs_comet/config.yml   # another dataset
# Extra arguments are passed on to snakemake.

# the config in use: --configfile if given, else the profile default
CONFIG=config/config.yml
for ((i = 1; i <= $#; i++)); do
    if [[ "${!i}" == "--configfile" ]]; then j=$((i + 1)); CONFIG="${!j}"; fi
done

# the snkmt monitoring database lives next to the results of that dataset
OUTPUT_PATH=$(sed -n 's/^output_path:[[:space:]]*//p' "$CONFIG" | sed 's/[[:space:]]*#.*//; s/["'\'']//g')
if [[ -z "$OUTPUT_PATH" ]]; then
    echo "ERROR: no output_path in $CONFIG"
    exit 1
fi
mkdir -p "$OUTPUT_PATH" log_files/smk_submission

# do not remove --sdm apptainer, otherwise the pipeline will not use the pixi container
module load singularityce
pixi run -e snakemake \
    snakemake --profile profiles/slurm/ --sdm apptainer \
    --logger-snkmt-db "$OUTPUT_PATH/snkmt.db" "$@"
