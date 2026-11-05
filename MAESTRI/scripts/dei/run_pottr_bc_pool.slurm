#!/bin/bash
# SLURM wrapper used for the paper's runs on the authors' cluster (Apptainer image from mastro.def/pottr.def); adapt partition and paths to your site.
#SBATCH --job-name=pottr-bcpool
#SBATCH --output=%x.k%a.out
#SBATCH --error=%x.k%a.err
#SBATCH --partition=lovelace
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=450G
#SBATCH --time=5-00:00:00
#SBATCH --array=2-101%8
# =============================================================================
# Breast POTTR pool sweep (POOL=50000, SWEEP_ONLY) on the lovelace nodes,
# inside ~/images/pottr.sif (built from scripts/dei/pottr.def). One array task
# per k, a whole node each (30 usable cores = 30 Gurobi threads); %8 keeps one
# task per lovelace node and the scheduler hands the next k to the first node
# that frees up, low k (the slow ones) first. 101 is the last k at which a pair
# is shared.
#
# The output dir follows the thread count: once the pool is full, its content
# depends on it.
#
# Needs, once:
#   - the matched cohort in results/pottr_cmp/dags and results/pottr_cmp/bc_inputs
#     (breastcancer_to_pottr.py --max_trees 2 --seed 0, as stage (A) of
#     run_pottr_breastcancer.sh builds it)
#   - a named academic Gurobi license per node in ~/gurobi/<node>.lic
#   - POTTR cloned next to the repo (~/POTTR)
#
# Gurobi's pool enumeration keeps only 4-7 threads busy (k=51-58 at 30
# threads: CPU efficiency 14-24%, below the DEI 60% auto-cancel policy), so
# pack 5 k per node at 6 threads each. The paper's run is the first line below;
# k = 2 did not finish at 6 threads, and pottr_compare_report.py reads it from
# the 30-thread run of the second line.
#
# Submit from <repo>/MAESTRI:
#   sbatch --cpus-per-task=6 --mem=90G --array=2-101%40 scripts/dei/run_pottr_bc_pool.slurm
#   sbatch --array=2 scripts/dei/run_pottr_bc_pool.slurm     # k = 2, 30 threads
# Stages B2-D then run once, outside SLURM: SKIP_PREP=1 bash
# scripts/run_pottr_breastcancer.sh. Stage B skips every k that has output; at 6
# threads k = 2 has none, so stage B starts it unless that output is in place.
#
# Resumable: a k with converted_graphs.txt is skipped.
# =============================================================================
set -euo pipefail

SIF=${SIF:-$HOME/images/pottr.sif}
POTTR_OUT=${POTTR_OUT:-results/pottr_cmp/pottr_bc_pool_t${SLURM_CPUS_PER_TASK}}
cd "$SLURM_SUBMIT_DIR"
[ -f scripts/run_pottr_breastcancer.sh ] || { echo "submit from <repo>/MAESTRI" >&2; exit 1; }
[ -f "$SIF" ] || { echo "image $SIF not found" >&2; exit 1; }
[ -f results/pottr_cmp/dags/P0-0_bc.txt ] || { echo "cohort results/pottr_cmp/dags missing" >&2; exit 1; }

# A named academic license is bound to one node: there is no fallback, the WLS
# license allows only 2 concurrent sessions
export GRB_LICENSE_FILE=$HOME/gurobi/$(hostname -s).lic
[ -f "$GRB_LICENSE_FILE" ] || { echo "Gurobi license $GRB_LICENSE_FILE not found" >&2; exit 1; }

K=${SLURM_ARRAY_TASK_ID}
echo "node=$(hostname -s) k=$K cores=$SLURM_CPUS_PER_TASK out=$POTTR_OUT license=$GRB_LICENSE_FILE"
SWEEP_ONLY=1 SKIP_PREP=1 POOL=50000 CORES=$SLURM_CPUS_PER_TASK POTTR_ENV="" \
  POTTR_OUT="$POTTR_OUT" K_LIST="$K" \
  srun apptainer exec "$SIF" bash scripts/run_pottr_breastcancer.sh
