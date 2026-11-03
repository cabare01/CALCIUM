#!/usr/bin/env bash
#SBATCH --job-name=stage14_preprocess
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --output=/sciclone/home/cabare01/CALCIUM/logs/%x-%j.out
#SBATCH --error=/sciclone/home/cabare01/CALCIUM/logs/%x-%j.err
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=cabare01@wm.edu

# Submit from Bora with:
# sbatch /sciclone/home/cabare01/CALCIUM/scripts/raw_prep/preprocess_stage_14.sbatch

set -euo pipefail

SCRIPT=/sciclone/home/cabare01/CALCIUM/scripts/raw_prep/preprocess_stage_14.py
python3 "${SCRIPT}"
