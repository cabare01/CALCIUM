#!/usr/bin/env bash
#SBATCH --job-name=rhythm
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/sciclone/home/cabare01/CALCIUM/logs/%x-%j.out
#SBATCH --error=/sciclone/home/cabare01/CALCIUM/logs/%x-%j.err
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=cabare01@wm.edu

# Step 1 red-noise test, one job per group. Submit from Bora (module load slurm/bora):
#   for g in S13 S13_devitellinated S14 S14_devitellinated; do
#     sbatch --job-name=rhythm_$g --export=ALL,GROUP=$g run_rhythmicity.sbatch
#   done
# Pilot (one embryo, 50 cells) before the full run:
#   sbatch --time=01:00:00 --export=ALL,GROUP=S14,EXTRA="--embryos 180209SP_E4_st_14 --max-cells 50" run_rhythmicity.sbatch
# Finished embryos are skipped, so a timed-out job can simply be resubmitted.

set -euo pipefail
: "${GROUP:?set GROUP=S13|S13_devitellinated|S14|S14_devitellinated}"
: "${N_SURROGATES:=1000}"
: "${EXTRA:=}"

module load miniforge3
source /sciclone/apps/miniforge3-24.9.2-0/etc/profile.d/conda.sh
conda activate wavelet-analysis
# One BLAS thread per worker process; parallelism comes from the process pool.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

python /sciclone/home/cabare01/CALCIUM/scripts/step1_rhythmicity/rhythmicity_test.py \
  --group "${GROUP}" --n-surrogates "${N_SURROGATES}" \
  --processes "${SLURM_CPUS_PER_TASK}" ${EXTRA}
