#!/bin/bash
# SLURM wrapper used for the paper's runs on the authors' cluster (Apptainer image from mastro.def/pottr.def); adapt partition and paths to your site.
#SBATCH --job-name=mastro-discovery
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=allgroups
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=94
#SBATCH --mem=500G
#SBATCH --time=2-00:00:00
# =============================================================================
# scripts/run_exp_discovery_<COHORT>.sh on one DEI compute node, inside
# ~/images/mastro.sif (built from scripts/dei/mastro.def). PAR follows
# --cpus-per-task; COHORT is breastcancer (default) or tracerx.
#
# Submit from the code directory (<repo>/MAESTRI); the knobs of the discovery
# script pass through the environment:
#   cd ~/MAESTRI/MAESTRI
#   sbatch scripts/dei/run_discovery.slurm
#   COHORT=tracerx sbatch --cpus-per-task=30 --mem=64G scripts/dei/run_discovery.slurm
#
# Mining and the observed p-values resume; the Westfall-Young step does not, it
# writes its outputs only at the end.
# =============================================================================
set -euo pipefail

SIF=${SIF:-$HOME/images/mastro.sif}
COHORT=${COHORT:-breastcancer}
cd "$SLURM_SUBMIT_DIR"
[ -f "scripts/run_exp_discovery_${COHORT}.sh" ] || { echo "submit from <repo>/MAESTRI" >&2; exit 1; }
[ -f "$SIF" ] || { echo "image $SIF not found" >&2; exit 1; }

export PAR=${PAR:-$SLURM_CPUS_PER_TASK}
# WY mining scratch in RAM on the node, not in the 50 GB home quota: at sigma=2
# a single null dataset can write tens of GB of LCM output. tmpfs counts
# against --mem.
export WORKDIR=${WORKDIR:-/dev/shm/$USER-$SLURM_JOB_ID}
mkdir -p "$WORKDIR"
trap 'rm -rf "$WORKDIR"' EXIT
echo "node=$(hostname) cpus=$SLURM_CPUS_PER_TASK PAR=$PAR COHORT=$COHORT EXACT_WHEN_CHEAPER=${EXACT_WHEN_CHEAPER:-1} MC_SAMPLES=${MC_SAMPLES:-default} OUT_ROOT=${OUT_ROOT:-default} WORKDIR=$WORKDIR"
srun apptainer exec "$SIF" bash "scripts/run_exp_discovery_${COHORT}.sh"
