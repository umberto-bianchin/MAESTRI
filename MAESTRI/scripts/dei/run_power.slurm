#!/bin/bash
# SLURM wrapper used for the paper's runs on the authors' cluster (Apptainer image from mastro.def/pottr.def); adapt partition and paths to your site.
#SBATCH --job-name=mastro-power
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=allgroups
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
# =============================================================================
# scripts/run_exp_power.sh with the script defaults, i.e. the paper's per-cell
# configuration, on one DEI compute node inside ~/images/mastro.sif (built from
# scripts/dei/mastro.def). PAR follows --cpus-per-task; 30 cpus fit the 32-core
# runners (2 cores per node are reserved). run_power_cell.slurm sets the same
# values explicitly.
#
# Submit from the code directory (<repo>/MAESTRI); the knobs of run_exp_power.sh
# pass through the environment:
#   cd ~/MAESTRI/MAESTRI
#   sbatch scripts/dei/run_power.slurm
#
# Resumable: the per-cell thresholds and the per-trial minimum p-values are
# cached under <OUT_BASE>/implant_theta*/.
# =============================================================================
set -euo pipefail

SIF=${SIF:-$HOME/images/mastro.sif}
cd "$SLURM_SUBMIT_DIR"
[ -f scripts/run_exp_power.sh ] || { echo "submit from <repo>/MAESTRI" >&2; exit 1; }
[ -f "$SIF" ] || { echo "image $SIF not found" >&2; exit 1; }

export PAR=${PAR:-$SLURM_CPUS_PER_TASK}
echo "node=$(hostname) cpus=$SLURM_CPUS_PER_TASK PAR=$PAR MTREES=${MTREES:-10} OUT_BASE=${OUT_BASE:-results/power_v3_cell}"
srun apptainer exec "$SIF" bash scripts/run_exp_power.sh
