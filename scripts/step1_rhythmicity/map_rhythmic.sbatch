#!/usr/bin/env bash
#SBATCH --job-name=map_rhythmic
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=/sciclone/home/cabare01/CALCIUM/logs/%x-%j.out
#SBATCH --error=/sciclone/home/cabare01/CALCIUM/logs/%x-%j.err
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=cabare01@wm.edu

# Step 1c: rerun the unchanged v3 mapping + combine step on rhythmic cells only.
# Inputs come from build_rhythmic_tables.py; output goes to a NEW folder, so the
# all-cell maps in mapping_v3 are untouched. Submit from Bora after step 1b.

set -euo pipefail
module load miniforge3
source /sciclone/apps/miniforge3-24.9.2-0/etc/profile.d/conda.sh
conda activate wavelet-analysis

TABLES=/sciclone/home/cabare01/scr10/CALCIUM/data/analysis_v3/step1_rhythmicity/tables_rhythmic
OUT=/sciclone/home/cabare01/scr10/CALCIUM/data/mapping_v3_rhythmic
MAP=/sciclone/home/cabare01/CALCIUM/scripts/mapping
S14=/sciclone/data10/cabare01/boundary_info/SuperEmbryoBoundaries/1_new_Stage_14_Sox2_Boundary.csv
S15=/sciclone/data10/cabare01/boundary_info/SuperEmbryoBoundaries/1_new_Stage_15_Sox2_Boundary.csv

declare -A FOLDER=([S13]=stage_13 [S13_devitellinated]=stage_13_devitellinated
                   [S14]=stage_14 [S14_devitellinated]=stage_14_devitellinated)
for g in S13 S13_devitellinated S14 S14_devitellinated; do
  python "${MAP}/spatial_spectral_mapping.py" \
    --input "${TABLES}/${FOLDER[$g]}_CW.csv" --group "$g" \
    --boundary-root /sciclone/home/cabare01/data10/boundary_info/Boundaries \
    --stage14-sox2 "${S14}" --stage15-sox2 "${S15}" --output-dir "${OUT}"
done
python "${MAP}/combine_group_maps.py" --mapping-root "${OUT}" \
  --stage14-sox2 "${S14}" --stage15-sox2 "${S15}" \
  --output-dir "${OUT}/combined" --min-embryos 2
