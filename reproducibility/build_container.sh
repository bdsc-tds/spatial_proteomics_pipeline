#!/bin/bash

# Build the pipeline's pixi container (default / spacec / tacit / gpu / snakemake environments):
#   sbatch reproducibility/build_container.sh
# The image is written to OUTPUT_DIR (default: resources/ in the repository; override with
# OUTPUT_DIR=/path sbatch ...), which must match `containers: pixi:` in the dataset config. The build
# runs in BUILD_DIR (default: $TMPDIR, else /tmp), which needs ~30 GB and must allow --fakeroot builds.

## Choose paths/names of output (stdout) and error (stderr) logs (relative to the submission folder)
#SBATCH --output build_spatial_proteomics_container_%j.out
#SBATCH --error build_spatial_proteomics_container_%j.err

#SBATCH --job-name=build_spatial_proteomics_container
#SBATCH --partition cpu

#SBATCH --time 02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=32GB

# apptainer (or singularity) must be on PATH; on clusters providing it as an environment module, load it
if ! command -v apptainer &>/dev/null && ! command -v singularity &>/dev/null; then
    module load apptainer 2>/dev/null || module load singularityce
fi
# the build command: apptainer, or singularity on older installations
BUILDER=$(command -v apptainer || command -v singularity)

# Exit immediately if a command exits with a non-zero status
set -e

# ================= CONFIGURATION =================
# 1. Auto-detect Root: first candidate that actually contains pixi.toml wins.
# Under sbatch the script is executed from a spool copy (/var/spool/slurm/...), so its own location
# is useless there and SLURM_SUBMIT_DIR is what points back at the repo; when run directly, pwd or
# the parent of this script's directory does. Explicit "$1" always overrides.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)"
PROJECT_ROOT=""
for candidate in "$1" "$SLURM_SUBMIT_DIR" "$(pwd)" "$(dirname "$SCRIPT_DIR")"; do
    if [[ -n "$candidate" && -f "$candidate/pixi.toml" ]]; then
        PROJECT_ROOT="$candidate"
        break
    fi
done

if [[ -z "$PROJECT_ROOT" ]]; then
    echo "ERROR: could not locate a directory containing pixi.toml."
    echo "       Tried: '$1', '\$SLURM_SUBMIT_DIR'='$SLURM_SUBMIT_DIR', pwd='$(pwd)', '$(dirname "$SCRIPT_DIR")'"
    echo "       Pass the project root explicitly: sbatch $0 /path/to/spatial_proteomics_pipeline"
    exit 1
fi
echo ">>> Using project root: $PROJECT_ROOT"

# 2. Define expected file locations
DEF_REL_PATH="reproducibility/pixi.def"
OUTPUT_DIR="${OUTPUT_DIR:-$PROJECT_ROOT/resources}"
IMAGE_NAME="pixi.sif"

FULL_DEF_PATH="${PROJECT_ROOT}/${DEF_REL_PATH}"

# 3. Sanity Checks
if [[ ! -f "$PROJECT_ROOT/pixi.toml" ]]; then
    echo "ERROR: pixi.toml not found in $PROJECT_ROOT"
    exit 1
fi

if [[ ! -f "$FULL_DEF_PATH" ]]; then
    echo "ERROR: Definition file not found at $FULL_DEF_PATH"
    exit 1
fi

if [[ ! -f "$PROJECT_ROOT/external/virtues/pyproject.toml" ]]; then
    echo "ERROR: VirTues submodule missing; run: git submodule update --init"
    exit 1
fi

# 4. Setup Scratch
SCRATCH_WORK_DIR="${BUILD_DIR:-${TMPDIR:-/tmp}}/build_pixi_$(date +%s)"
# =================================================

# Save starting location (though script execution creates a subshell, this is good practice)
ORIGINAL_DIR=$(pwd)

echo ">>> Setting up build environment in: $SCRATCH_WORK_DIR"
mkdir -p "$SCRATCH_WORK_DIR"

# -------------------------------------------------
# STAGE 1: Copy files to Scratch
# -------------------------------------------------
echo ">>> Staging files..."
cp "${PROJECT_ROOT}/pixi.toml" "$SCRATCH_WORK_DIR/"
cp "${PROJECT_ROOT}/pixi.lock" "$SCRATCH_WORK_DIR/"
# VirTues source exactly as in the pinned submodule commit (no .git, no local edits)
mkdir -p "$SCRATCH_WORK_DIR/external/virtues"
git -C "${PROJECT_ROOT}/external/virtues" archive HEAD | tar -x -C "$SCRATCH_WORK_DIR/external/virtues"
cp "$FULL_DEF_PATH" "$SCRATCH_WORK_DIR/singularity.def"


# -------------------------------------------------
# STAGE 2: Build
# -------------------------------------------------
echo ">>> Changing directory to scratch..."
cd "$SCRATCH_WORK_DIR"

echo ">>> Starting container build with $BUILDER (this may take time)..."
# --fakeroot needs a local or scratch file system (it often fails on network/project storage)
"$BUILDER" build --fakeroot --force "$IMAGE_NAME" singularity.def

# -------------------------------------------------
# STAGE 3: Finalize
# -------------------------------------------------
echo ">>> Moving image back to project: $OUTPUT_DIR/$IMAGE_NAME"
mkdir -p "$OUTPUT_DIR"
mv "$IMAGE_NAME" "$OUTPUT_DIR/"

echo ">>> Cleaning up scratch..."
cd "$ORIGINAL_DIR"
rm -rf "$SCRATCH_WORK_DIR"

echo ">>> Done. Your container is ready at: $OUTPUT_DIR/$IMAGE_NAME"
