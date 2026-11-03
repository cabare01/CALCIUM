#!/usr/bin/env bash
#SBATCH --job-name=spatial_combine
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=/sciclone/home/cabare01/CALCIUM/logs/%x-%j.out
#SBATCH --error=/sciclone/home/cabare01/CALCIUM/logs/%x-%j.err
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=cabare01@wm.edu

# Run after all four spatial_mapping_stage_*.sbatch jobs have finished, e.g.
#   sbatch --dependency=afterok:<id1>:<id2>:<id3>:<id4> combine_group_maps.sbatch
# Optional: EXCLUDE="200114SP_E1_st_13001,180405SP_E4_st_13" for a sensitivity run
# (set OUTPUT to a different folder as well).

set -euo pipefail

module load miniforge3
source /sciclone/apps/miniforge3-24.9.2-0/etc/profile.d/conda.sh
conda activate wavelet-analysis

MAPPING_ROOT=${MAPPING_ROOT:-/sciclone/home/cabare01/scr10/CALCIUM/data/mapping_v3}
OUTPUT=${OUTPUT:-${MAPPING_ROOT}/combined}
EXCLUDE=${EXCLUDE:-}

python /sciclone/home/cabare01/CALCIUM/scripts/mapping/combine_group_maps.py \
  --mapping-root "${MAPPING_ROOT}" \
  --stage14-sox2 /sciclone/data10/cabare01/boundary_info/SuperEmbryoBoundaries/1_new_Stage_14_Sox2_Boundary.csv \
  --stage15-sox2 /sciclone/data10/cabare01/boundary_info/SuperEmbryoBoundaries/1_new_Stage_15_Sox2_Boundary.csv \
  --output-dir "${OUTPUT}" \
  --min-embryos "${MIN_EMBRYOS:-2}" \
  --exclude-embryos "${EXCLUDE}"
