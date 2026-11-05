#!/bin/bash
# SLURM wrapper used for the paper's runs on the authors' cluster (Apptainer image from mastro.def/pottr.def); adapt partition and paths to your site.
#SBATCH --job-name=timing
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --partition=lovelace
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --time=1-00:00:00
#SBATCH --array=0-4
# =============================================================================
# POTTR against MAESTRI wall-clock timing on whole lovelace nodes. One array
# task per cohort, a whole node each (--exclusive, --mem=0 = all the node's RAM),
# so that no other job shares the cores while timing. POTTR runs with its
# default settings (no solution pool).
#
#   task 0: cap1   task 1: cap2   task 2: cap10   task 3: drop5   (breast)
#   task 4: TRACERx, the matched cohort of run_pottr_tracerx.sh (469 trees,
#           built there with PREP_ONLY=1 when missing), POTTR also at 1 thread
#
# Each task writes its own OUTDIR, so that the timing.csv rebuilds of
# concurrent tasks never race; timing_table.py --indir results/timing_dei
# merges them afterwards.
#
# Needs, once: a named academic Gurobi license per node in ~/gurobi/<node>.lic,
# POTTR cloned next to the repo (~/POTTR) and ../data/breastCancer.npy. The x86
# lcm53/lcm binary must be built inside the image (it is built on first use).
#
# Submit from <repo>/MAESTRI:
#   sbatch scripts/dei/run_timing.slurm                 # all five cohorts
#   sbatch --array=1 scripts/dei/run_timing.slurm       # cap2 only
#   POTTR_THREADS="20 30" sbatch scripts/dei/run_timing.slurm
#
# Resumable: a cell with its JSON record is skipped.
# =============================================================================
set -euo pipefail

SIF=${SIF:-$HOME/images/pottr.sif}
cd "$SLURM_SUBMIT_DIR"
[ -f scripts/run_timing_pottr_vs_mastro.sh ] || { echo "submit from <repo>/MAESTRI" >&2; exit 1; }
[ -f "$SIF" ] || { echo "image $SIF not found" >&2; exit 1; }

export GRB_LICENSE_FILE=$HOME/gurobi/$(hostname -s).lic
[ -f "$GRB_LICENSE_FILE" ] || { echo "Gurobi license $GRB_LICENSE_FILE not found" >&2; exit 1; }

COHORTS=(1 2 10 d5 tracerx)
CAP=${COHORTS[$SLURM_ARRAY_TASK_ID]}
export POTTR_THREADS=${POTTR_THREADS:-"$SLURM_CPUS_PER_TASK"}
if [ "$CAP" = tracerx ]; then
  LABEL=tracerx
  export COHORT_DIR=${TX_COHORT:-results/pottr_cmp/tracerx_cohort} COHORT_NAME=tracerx
  [ -f "$COHORT_DIR/inputs/graphs_all.txt" ] || \
    CDIR="$COHORT_DIR" PREP_ONLY=1 POTTR_ENV="" srun apptainer exec "$SIF" bash scripts/run_pottr_tracerx.sh
  POTTR_THREADS="1 $POTTR_THREADS"
else
  LABEL=$([ "${CAP#d}" != "$CAP" ] && echo "drop${CAP#d}" || echo "cap${CAP}")
  export CAP_LIST="$CAP"
fi

export OUTDIR=${OUTDIR_BASE:-results/timing_dei}/${LABEL}
export K_LIST=${K_LIST:-"2 5 10 20 50"}
export CORES=$SLURM_CPUS_PER_TASK
export TIMEOUT=${TIMEOUT:-14400}
export MAX_CONSEC_TIMEOUTS=0
export POOL=0
export POTTR_ENV=""
echo "node=$(hostname -s) cohort=$LABEL cpus=$SLURM_CPUS_PER_TASK pottr_threads=$POTTR_THREADS k=[$K_LIST] timeout=${TIMEOUT}s out=$OUTDIR"
srun apptainer exec "$SIF" bash scripts/run_timing_pottr_vs_mastro.sh
