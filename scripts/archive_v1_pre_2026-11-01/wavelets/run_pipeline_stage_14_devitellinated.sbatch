#!/usr/bin/env bash
#SBATCH --job-name=st14dev_wavelets
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/sciclone/home/cabare01/CALCIUM/logs/%x-%j.out
#SBATCH --error=/sciclone/home/cabare01/CALCIUM/logs/%x-%j.err
#SBATCH --mail-type=FAIL,BEGIN,END # when to email you
#SBATCH --mail-user=cabare01@wm.edu # who to email

# sbatch --export=ALL,CALCIUM_RAW_TRACE_ROOT=/sciclone/home/cabare01/data10/merged_data/stage_14_devitellinated,CALCIUM_OUTPUT_ROOT=/sciclone/home/cabare01/data10/HPC_spectral/stage_14_devitellinated \
#   /sciclone/home/cabare01/CALCIUM/scripts/wavelets/run_pipeline_stage_14_devitellinated.sbatch

set -euo pipefail

: "${CALCIUM_RAW_TRACE_ROOT:=/sciclone/home/cabare01/data10/merged_data/stage_14_devitellinated}"
: "${CALCIUM_OUTPUT_ROOT:=/sciclone/home/cabare01/data10/HPC_spectral/stage_14_devitellinated}"
export CALCIUM_RAW_TRACE_ROOT CALCIUM_OUTPUT_ROOT

module load miniforge3
source /sciclone/apps/miniforge3-24.9.2-0/etc/profile.d/conda.sh
conda activate wavelet-analysis

SCRIPT=/sciclone/home/cabare01/CALCIUM/scripts/wavelets/pipeline.py
python "${SCRIPT}"
