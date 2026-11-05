#!/bin/bash
# SLURM wrapper used for the paper's runs on the authors' cluster (Apptainer image from mastro.def/pottr.def); adapt partition and paths to your site.
#SBATCH --job-name=mastro-power-cell
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=allgroups
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
# =============================================================================
# Power experiment with one WY threshold per (N, f) cell, calibrated on an
# implanted cohort of that cell (implant_experiment_multitree.py --calibration
# cell), on one DEI compute node inside ~/images/mastro.sif. Grid of the paper:
# 60 patients, M = 10, k = 4, sigma = 2, N in {2,..,30}, f in {0.1,..,1.0},
# 200 trials per cell, theta in {0.5, 1.0}, 500 calibration resamples per cell.
#
# Submit from the code directory (<repo>/MAESTRI):
#   cd ~/MAESTRI/MAESTRI
#   sbatch scripts/dei/run_power_cell.slurm
#
# Resumable: the per-cell thresholds are cached in
# <OUT_BASE>/implant_theta*/_cal_m<WY_M>/ and the per-trial minimum p-values in
# <OUT_BASE>/implant_theta*/trials.csv. A second submission with another WY_M
# and the same OUT_BASE recalibrates only and reuses every trial.
# =============================================================================
set -euo pipefail

SIF=${SIF:-$HOME/images/mastro.sif}
cd "$SLURM_SUBMIT_DIR"
[ -f scripts/run_exp_power.sh ] || { echo "submit from <repo>/MAESTRI" >&2; exit 1; }
[ -f "$SIF" ] || { echo "image $SIF not found" >&2; exit 1; }

export PAR=${PAR:-$SLURM_CPUS_PER_TASK}
export CALIBRATION=cell
export WY_M=${WY_M:-500}
export F_LIST=${F_LIST:-0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
export OUT_BASE=${OUT_BASE:-results/power_v3_cell}
echo "node=$(hostname) cpus=$SLURM_CPUS_PER_TASK PAR=$PAR WY_M=$WY_M OUT_BASE=$OUT_BASE start=$(date -Is)"
srun apptainer exec "$SIF" bash scripts/run_exp_power.sh
echo "end=$(date -Is)"
