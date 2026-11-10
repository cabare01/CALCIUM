#!/usr/bin/env bash
#SBATCH --job-name=pdt_wavelet_panels
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --output=/sciclone/home/cabare01/CALCIUM/logs/%x-%j.out
#SBATCH --error=/sciclone/home/cabare01/CALCIUM/logs/%x-%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=cabare01@wm.edu

# SfN poster panels 2 and 3: all 204 PDT setups vs wavelet features, red noise
# and membrane removal. Reads ~3 GB of PDT per-cell tables + the 2.9 GB
# magnitude file (chunked). The log records the join diagnostics.
# Submit from Bora:
#   bash -lc 'module load slurm/bora; sbatch /sciclone/home/cabare01/CALCIUM/scripts/step6_figures/pdt_wavelet_panels.sbatch'
# Override the output folder with --export=ALL,OUT=/path (must be a new folder).

set -euo pipefail
module load miniforge3
source /sciclone/apps/miniforge3-24.9.2-0/etc/profile.d/conda.sh
conda activate wavelet-analysis

OUT=${OUT:-/sciclone/home/cabare01/scr10/CALCIUM/data/analysis_v3/step6_figures/pdt_wavelet_panels}
if [ -e "${OUT}" ] && [ -n "$(ls -A "${OUT}")" ]; then
  echo "Output folder ${OUT} already has files; refusing to overwrite. Set OUT to a new folder." >&2
  exit 1
fi

python -u /sciclone/home/cabare01/CALCIUM/scripts/step6_figures/pdt_wavelet_panels.py --out "${OUT}"
# poster-style (v2) panels 2, 3A, 3B from the tables above -> ${OUT}/poster_v2
python -u /sciclone/home/cabare01/CALCIUM/scripts/step6_figures/pdt_wavelet_panels_styled.py --in "${OUT}"
