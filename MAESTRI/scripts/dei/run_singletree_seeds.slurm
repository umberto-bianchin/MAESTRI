#!/bin/bash
# SLURM wrapper used for the paper's runs on the authors' cluster (Apptainer image from mastro.def/pottr.def); adapt partition and paths to your site.
#SBATCH --job-name=singletree-seeds
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --partition=lovelace
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=0-9
# =============================================================================
# scripts/run_exp_singletree_seeds.sh (MASTRO on one tree drawn per patient,
# 10 draws) on the DEI cluster, inside ~/images/mastro.sif. One array task per
# seed, both sigmas; DATASET is breastCancer (default) or tracerx.
#
# Submit from the code directory (<repo>/MAESTRI):
#   cd ~/MAESTRI/MAESTRI
#   DATASET=tracerx sbatch scripts/dei/run_singletree_seeds.slurm
#
# Every task ends with the aggregation over the seeds already on disk; the one
# that matters is the run after the last task, so once all are done:
#   python3 analyze_singletree_seeds.py --root results/singletree_seeds_tracerx \
#     --sigma_list "2 5" --alpha 0.05
# Resumable: a (sigma, seed) with its outputs is skipped.
# =============================================================================
set -euo pipefail

SIF=${SIF:-$HOME/images/mastro.sif}
cd "$SLURM_SUBMIT_DIR"
[ -f scripts/run_exp_singletree_seeds.sh ] || { echo "submit from <repo>/MAESTRI" >&2; exit 1; }
[ -f "$SIF" ] || { echo "image $SIF not found" >&2; exit 1; }

export SEEDS=${SLURM_ARRAY_TASK_ID}
export PAR=${PAR:-$SLURM_CPUS_PER_TASK}
echo "node=$(hostname -s) seed=$SEEDS cpus=$SLURM_CPUS_PER_TASK DATASET=${DATASET:-breastCancer}"
srun apptainer exec "$SIF" bash scripts/run_exp_singletree_seeds.sh
