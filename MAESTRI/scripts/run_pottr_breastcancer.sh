#!/bin/bash
# POTTR on breastCancer, then its trajectories re-scored under the MAESTRI test
# on the same trees. POTTR builds a conflict graph over every pair of trees
# across patients: on all 37809 trees that is about 7e8 pairs and the stage does
# not finish, so the driver first builds a deduplicated cohort capped at
# MAX_TREES distinct trees per patient. Stages: (A) matched cohort, (B) POTTR
# k-sweep, (B2) POTTR's own p-values, (C) MAESTRI p-values, (D) figure.
# Defaults are the paper configuration.
set -euo pipefail
cd "$(dirname "$0")/.."                     # -> <repo>/MAESTRI, the code directory
CODE_DIR="$(pwd)"

NPY=${NPY:-../data/breastCancer.npy}
POTTR_REPO=${POTTR_REPO:-../../POTTR}        # POTTR clone next to this repo (holds code/)
POTTR_ENV=${POTTR_ENV:-pottr_env}            # conda env with gurobi (environment.yaml)
DAGS_DIR=${DAGS_DIR:-results/pottr_cmp/dags}            # POTTR input, one file per tree
ENS_DIR=${ENS_DIR:-results/pottr_cmp/bc_inputs}         # matched multi-tree inputs

# Paper run: whole cohort, 2 distinct trees per patient, k = 2..101 (101 is the
# last k at which a pair is shared), Gurobi solution pool of 50000 on 6 threads.
MAX_TREES=${MAX_TREES:-2}                    # distinct trees/patient cap (drives O(T^2))
N_PATIENTS=${N_PATIENTS:-0}                  # 0 = all eligible patients
MULTITREE_ONLY=${MULTITREE_ONLY:-0}          # 1 = only patients with >=2 distinct trees
K_MIN=${K_MIN:-2}                            # POTTR recurrence levels to sweep;
K_MAX=${K_MAX:-101}                          #   K_LIST defaults to K_MIN..K_MAX
K_LIST=${K_LIST:-"$(seq "$K_MIN" "$K_MAX" | tr '\n' ' ')"}
CORES=${CORES:-6}                            # Gurobi threads; the pool content depends on it
# >0: Gurobi solution pool of this size (the POTTR paper uses 50000), run
# through pottr_pool_logged.py, which also writes k<k>/pool_log.json
POOL=${POOL:-50000}
POTTR_OUT=${POTTR_OUT:-results/pottr_cmp/pottr_bc_pool_t${CORES}}   # <- k<k>/ subdirs written here
SIG_OUT=${SIG_OUT:-${POTTR_OUT}_significance.csv}
FORCE_POTTR_SIG=${FORCE_POTTR_SIG:-1}        # 1 = re-apply run_POTTR's significance
POTTR_SIG_MAX_NODES=${POTTR_SIG_MAX_NODES:-6} #   gate per trajectory, not per k
# 1 = DAGS_DIR and ENS_DIR already hold a matched cohort (run_pottr_tracerx.sh)
SKIP_PREP=${SKIP_PREP:-0}
# 1 = keep POTTR's '-?-' pairs as same-node relations (TRACERx encoding)
KEEP_CLUSTER=${KEEP_CLUSTER:-0}
MC_SAMPLES=${MC_SAMPLES:-2000}               # patients above 8 trees only; the capped
EXACT_WHEN_CHEAPER=${EXACT_WHEN_CHEAPER:-0}  #   breast cohort has none
THETA=${THETA:-1.0}
NULL=${NULL:-perm}
SEED=${SEED:-0}

MT_FLAG=""
[ "$MULTITREE_ONLY" = "1" ] && MT_FLAG="--multitree_only"

# ---- (A) build matched, tractable cohort ------------------------------------
if [ "$SKIP_PREP" = "1" ]; then
  echo "=== (A) using prebuilt cohort: dags=$DAGS_DIR inputs=$ENS_DIR ==="
else
  echo "=== (A) prepare matched POTTR + multi-tree cohort ==="
  python3 breastcancer_to_pottr.py \
    --npy "$NPY" --out "$DAGS_DIR" --multitree_out "$ENS_DIR" \
    --max_trees "$MAX_TREES" --n_patients "$N_PATIENTS" $MT_FLAG --seed "$SEED"
fi

# absolute paths for POTTR (run from POTTR/code with a different CWD)
DAGS_ABS="$(cd "$(dirname "$DAGS_DIR")" && pwd)/$(basename "$DAGS_DIR")"
POTTR_OUT_ABS="$(mkdir -p "$POTTR_OUT" && cd "$POTTR_OUT" && pwd)"

# Optional conda env for POTTR; POTTR_ENV="" uses the current python3. A
# missing env only prints a warning.
if [ -n "$POTTR_ENV" ] && command -v conda >/dev/null 2>&1; then
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"
  if conda env list | awk '{print $1}' | grep -qx "$POTTR_ENV"; then
    conda activate "$POTTR_ENV"
  else
    echo "[warn] conda env '$POTTR_ENV' not found; using current python3"
  fi
fi
# Fail fast with a clear message if POTTR's deps are not importable.
( cd "${POTTR_REPO}/code" && python3 -c "import gurobipy, networkx, pandas" ) \
  || { echo "[err] POTTR deps missing (gurobipy/networkx/pandas) in python3."; \
       echo "      Set POTTR_ENV=<conda env> or install them for system python3."; \
       exit 1; }

# ---- (B) POTTR k-sweep ------------------------------------------------------
echo "=== (B) POTTR k-sweep: k in [$K_LIST] ==="
for K in $K_LIST; do
  KOUT="${POTTR_OUT_ABS}/k${K}"
  if [ -f "${KOUT}/converted_graphs.txt" ]; then
    echo "[SKIP] k=$K already has converted_graphs.txt"
    continue
  fi
  echo "--- POTTR k=$K -> ${KOUT} ---"
  if [ "$POOL" -gt 0 ]; then
    ( cd "${POTTR_REPO}/code" && \
      python3 "${CODE_DIR}/pottr_pool_logged.py" --pool_log "${KOUT}/pool_log.json" \
        -o "$KOUT" -d "$DAGS_ABS" -k "$K" -c "$CORES" -parallel -pool "$POOL" )
  else
    ( cd "${POTTR_REPO}/code" && \
      python3 run_POTTR.py -o "$KOUT" -d "$DAGS_ABS" -k "$K" -c "$CORES" -parallel )
  fi
done

# 1 = stop after the sweep: several instances can share POTTR_OUT with disjoint
# K_LISTs, then one run with the full K_LIST skips every k and does B2-D once
if [ "${SWEEP_ONLY:-0}" = "1" ]; then
  echo "=== sweep only: done k in [$K_LIST] ==="
  exit 0
fi

KMIN=$(echo $K_LIST | tr ' ' '\n' | sort -n | head -1)
KMAX=$(echo $K_LIST | tr ' ' '\n' | sort -n | tail -1)

# ---- (B2) POTTR's own p-value for every trajectory it can score --------------
# run_POTTR.py runs its test only when the largest trajectory of that k has
# fewer than 13 nodes (compute_num_automorph enumerates n! permutations), so
# whole k directories, small trajectories included, come out without a p-value.
# Re-apply the gate per trajectory instead; POTTR's own files are left alone.
SIG_GLOBAL=""
if [ "$FORCE_POTTR_SIG" = "1" ]; then
  echo "=== (B2) POTTR significance test, per trajectory (<= $POTTR_SIG_MAX_NODES nodes) ==="
  SIG_GLOBAL="${POTTR_OUT}/significance_forced.txt"
  python3 pottr_force_significance.py \
    --pottr_dir "$POTTR_OUT" --k_range "${KMIN},${KMAX}" \
    --dags "$DAGS_ABS" --pottr_repo "$POTTR_REPO" \
    --max_nodes "$POTTR_SIG_MAX_NODES" --cores "$CORES" \
    --out "$SIG_GLOBAL"
fi

# ---- (C) re-score POTTR trajectories under the MAESTRI test -----------------
echo "=== (C) MAESTRI significance of POTTR trajectories ==="
python3 pottr_significance.py \
  --pottr_dir "$POTTR_OUT" --k_range "${KMIN},${KMAX}" \
  ${SIG_GLOBAL:+--pottr_sig_global "$SIG_GLOBAL"} \
  --graphs_all "${ENS_DIR}/graphs_all.txt" \
  -w "${ENS_DIR}/weights_uniform.txt" --owner "${ENS_DIR}/owner.txt" \
  --theta "$THETA" --null "$NULL" --n_jobs "$CORES" --seed "$SEED" \
  --mc_samples "$MC_SAMPLES" \
  $( [ "$EXACT_WHEN_CHEAPER" = 1 ] && echo --exact_when_cheaper ) \
  $( [ "$KEEP_CLUSTER" = 1 ] && echo --keep_cluster ) \
  --out "$SIG_OUT"

# ---- (D) figure ------------------------------------------------------------
PLOT_OUT=${PLOT_OUT:-results/pottr_cmp/pottr_comparison.pdf}
echo "=== (D) figure ==="
python3 plot_pottr_comparison.py --csv "$SIG_OUT" --out "$PLOT_OUT"

echo "=== done: POTTR dags=$DAGS_ABS  pottr_out=$POTTR_OUT_ABS  sig=$SIG_OUT  plot=$PLOT_OUT ==="
