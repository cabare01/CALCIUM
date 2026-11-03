#!/usr/bin/env bash
#SBATCH --job-name=analysis_step
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=/sciclone/home/cabare01/CALCIUM/logs/%x-%j.out
#SBATCH --error=/sciclone/home/cabare01/CALCIUM/logs/%x-%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=cabare01@wm.edu

# Generic wrapper for steps 1b-6 (each is a single Python script, minutes to an hour).
# Submit from Bora, naming the job after the step, e.g.:
#   S=/sciclone/home/cabare01/CALCIUM/scripts
#   sbatch --job-name=step2_flow   --export=ALL,CMD="python $S/step2_cell_flow/cell_flow.py" $S/analysis_common/run_step.sbatch
#   sbatch --job-name=step4_effects --export=ALL,CMD="python $S/step4_embryo_stats/group_effects.py --lmm" $S/analysis_common/run_step.sbatch

set -euo pipefail
: "${CMD:?set CMD to the python command to run}"
module load miniforge3
source /sciclone/apps/miniforge3-24.9.2-0/etc/profile.d/conda.sh
conda activate wavelet-analysis
echo "Running: ${CMD}"
eval "${CMD}"
