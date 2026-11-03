#!/usr/bin/env bash
#SBATCH --job-name=st13dev_wavelets
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/sciclone/home/cabare01/CALCIUM/logs/%x-%j.out
#SBATCH --error=/sciclone/home/cabare01/CALCIUM/logs/%x-%j.err
#SBATCH --mail-type=FAIL,BEGIN,END # when to email you
#SBATCH --mail-user=cabare01@wm.edu # who to email

# sbatch --export=ALL,CALCIUM_OUTPUT_ROOT=/sciclone/home/cabare01/data10/HPC_spectral_v3/stage_13_devitellinated \
#   /sciclone/home/cabare01/CALCIUM/scripts/wavelets/run_pipeline_stage_13_devitellinated.sbatch

set -euo pipefail

: "${CALCIUM_RAW_TRACE_ROOT:=/sciclone/home/cabare01/data10/merged_data}"   # searched recursively
: "${CALCIUM_GROUP:=stage_13_devitellinated}"   # embryos chosen from data10/Embryo_Ref_Table.csv
: "${CALCIUM_OUTPUT_ROOT:=/sciclone/home/cabare01/data10/HPC_spectral_v3/stage_13_devitellinated}"
: "${CALCIUM_MAKE_TRACE_PLOTS:=false}"   # set true to write one dF/F PNG per cell (~32 GB per run in v1)
export CALCIUM_GROUP CALCIUM_RAW_TRACE_ROOT CALCIUM_OUTPUT_ROOT CALCIUM_MAKE_TRACE_PLOTS

module load miniforge3
source /sciclone/apps/miniforge3-24.9.2-0/etc/profile.d/conda.sh
conda activate wavelet-analysis

SCRIPT=/sciclone/home/cabare01/CALCIUM/scripts/wavelets/pipeline.py
python "${SCRIPT}"
