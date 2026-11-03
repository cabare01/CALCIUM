#!/usr/bin/env bash
#SBATCH --job-name=spatial_S13
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --output=/sciclone/home/cabare01/CALCIUM/logs/%x-%j.out
#SBATCH --error=/sciclone/home/cabare01/CALCIUM/logs/%x-%j.err
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=cabare01@wm.edu

# Submit from Bora with:
# sbatch /sciclone/home/cabare01/CALCIUM/scripts/mapping/spatial_mapping_stage_13.sbatch

set -euo pipefail

module load miniforge3
source /sciclone/apps/miniforge3-24.9.2-0/etc/profile.d/conda.sh
conda activate wavelet-analysis

SCRIPT=/sciclone/home/cabare01/CALCIUM/scripts/mapping/spatial_spectral_mapping.py
python "${SCRIPT}" \
  --input /sciclone/home/cabare01/CALCIUM/data/master_tables/group_tables_v3/stage_13_CW.csv \
  --group S13 \
  --expected-embryos 46 \
  --boundary-root /sciclone/home/cabare01/data10/boundary_info/Boundaries \
  --stage14-sox2 /sciclone/data10/cabare01/boundary_info/SuperEmbryoBoundaries/1_new_Stage_14_Sox2_Boundary.csv \
  --stage15-sox2 /sciclone/data10/cabare01/boundary_info/SuperEmbryoBoundaries/1_new_Stage_15_Sox2_Boundary.csv \
  --output-dir /sciclone/home/cabare01/scr10/CALCIUM/data/mapping_v3
