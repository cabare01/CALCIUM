#!/usr/bin/env bash
#SBATCH --job-name=calcium_iii_pywt
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/sciclone/home/cabare01/CALCIUM/logs/%x-%j.out
#SBATCH --error=/sciclone/home/cabare01/CALCIUM/logs/%x-%j.err

set -euo pipefail

# Submit after activating/loading a Python environment with numpy, pandas,
# scipy, matplotlib, and PyWavelets. Set CALCIUM_WORKFLOW and all required
# input variables, plus a new CALCIUM_OUTPUT_ROOT, before using sbatch.
# Example:
# sbatch --export=ALL,CALCIUM_WORKFLOW=spectral,CALCIUM_RAW_TRACE_ROOT=/path/to/traces,CALCIUM_OUTPUT_ROOT=/path/to/new/run \
#   /sciclone/home/cabare01/CALCIUM/scripts/calcium_iii_pywt_workflow.sbatch

: "${CALCIUM_WORKFLOW:?Set CALCIUM_WORKFLOW before submission}"
: "${CALCIUM_OUTPUT_ROOT:?Set CALCIUM_OUTPUT_ROOT to a new path before submission}"

PYTHON_BIN="${CALCIUM_PYTHON:-python3}"
SCRIPT=/sciclone/home/cabare01/CALCIUM/scripts/calcium_iii_pywavelets_workflows.py
"${PYTHON_BIN}" "${SCRIPT}"
