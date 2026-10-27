#!/usr/bin/env bash
#SBATCH --job-name=calcium_iii
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/sciclone/home/cabare01/CALCIUM/logs/%x-%j.out
#SBATCH --error=/sciclone/home/cabare01/CALCIUM/logs/%x-%j.err

# Submit with, for example:
# sbatch --export=ALL,CALCIUM_WORKFLOW=spectral,CALCIUM_RAW_TRACE_ROOT=/path/to/traces,CALCIUM_OUTPUT_ROOT=/path/to/new/run \
#   /sciclone/home/cabare01/CALCIUM/scripts/calcium_iii_workflow.sbatch
# The MATLAB entry point refuses an output path that already exists.

set -ueo pipefail

: "${CALCIUM_WORKFLOW:?Set CALCIUM_WORKFLOW before submission}"
: "${CALCIUM_OUTPUT_ROOT:?Set CALCIUM_OUTPUT_ROOT to a new path before submission}"

module load matlab/R2024a
SCRIPT_DIR=/sciclone/home/cabare01/CALCIUM/scripts
matlab -batch "addpath('${SCRIPT_DIR}'); calcium_iii_hpc_workflows"
