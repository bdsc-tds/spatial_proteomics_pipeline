#!/bin/bash -l
# Run the pipeline on SLURM with the active dataset config (config/config.yml), e.g.
#   ./run.sh -n                                              # dry run
#   ./run.sh --configfile /path/to/my_dataset/config.yml     # another dataset
#   ./run.sh --profile /path/to/my_profile                   # own profile instead of profiles/slurm/
# Extra arguments are passed on to snakemake.

# the config in use: --configfile if given, else the profile default
CONFIG=config/config.yml
# the profile: --profile if given, else $SNAKEMAKE_PROFILE (read by snakemake itself),
# else profiles/local/ (gitignored, e.g. a symlink to your own profile), else profiles/slurm/
PROFILE_ARGS=(--profile profiles/slurm/)
[[ -d profiles/local ]] && PROFILE_ARGS=(--profile profiles/local/)
# printenv: only an exported variable reaches snakemake (a plain shell variable from ~/.bashrc does not)
[[ -n "$(printenv SNAKEMAKE_PROFILE)" ]] && PROFILE_ARGS=()
for ((i = 1; i <= $#; i++)); do
    if [[ "${!i}" == "--configfile" ]]; then j=$((i + 1)); CONFIG="${!j}"; fi
    if [[ "${!i}" == "--profile" ]]; then PROFILE_ARGS=(); fi
done

# the snkmt monitoring database lives next to the results of that dataset
OUTPUT_PATH=$(sed -n 's/^output_path:[[:space:]]*//p' "$CONFIG" | sed 's/[[:space:]]*#.*//; s/["'\'']//g')
if [[ -z "$OUTPUT_PATH" ]]; then
    echo "ERROR: no output_path in $CONFIG"
    exit 1
fi
mkdir -p "$OUTPUT_PATH" log_files/smk_submission

# apptainer (or singularity) must be on PATH; on clusters providing it as an environment module, load it
if ! command -v apptainer &>/dev/null && ! command -v singularity &>/dev/null; then
    module load apptainer 2>/dev/null || module load singularityce
fi

# do not remove --sdm apptainer, otherwise the pipeline will not use the pixi container
pixi run -e snakemake \
    snakemake "${PROFILE_ARGS[@]}" --sdm apptainer \
    --logger-snkmt-db "$OUTPUT_PATH/snkmt.db" "$@"
