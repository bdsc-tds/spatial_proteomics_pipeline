#!/bin/bash --login

# script to run image registration with insitupy

## Choose paths/names of output (stdout) and error (stderr) logs
#SBATCH --output cropped_logs/out.out
#SBATCH --error cropped_logs/err.err

##SBATCH --partition gpu
##SBATCH --gres gpu:1
##SBATCH --gres-flags enforce-binding
#SBATCH --partition cpu

#SBATCH --time 04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=256GB

conda init
conda activate insitupy
python image_cropping.py\
    '/path/to/data'\
    '/path/to/data'\
    '/path/to/data'

python create_alignment.py\
    /path/to/data\
    '/path/to/data'\
    '/path/to/data'\
    --min_good_matches_per_area 1
# '/path/to/data'\
# '/path/to/data'\