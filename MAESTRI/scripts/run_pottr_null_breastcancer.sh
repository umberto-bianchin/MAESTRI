#!/bin/bash
# POTTR on null cohorts: the matched breastCancer cohort of
# run_pottr_breastcancer.sh, relabeled per patient by one random bijection of
# its alterations shared by all its trees (permute_pottr_cohort.py). Per seed:
# POTTR k-sweep, POTTR's own p-values, MAESTRI p-values of the same
# trajectories; then a per-seed count below alpha (summarize_pottr_null.py).
# Defaults are the paper configuration.
set -euo pipefail
cd "$(dirname "$0")/.."                     # -> <repo>/MAESTRI, the code directory

NPY=${NPY:-../data/breastCancer.npy}
POTTR_REPO=${POTTR_REPO:-../../POTTR}        # POTTR clone next to this repo (holds code/)
POTTR_ENV=${POTTR_ENV:-pottr_env}            # conda env with gurobi (environment.yaml)
OUT_BASE=${OUT_BASE:-results/pottr_null}
MAX_TREES=${MAX_TREES:-2}                    # same cohort as the real-data run
SEEDS=${SEEDS:-"1 2 3 4 5 6 7 8 9 10"}
K_LIST=${K_LIST:-"2 3 5 8 12 17 22 29"}
CORES=${CORES:-20}
POTTR_SIG_MAX_NODES=${POTTR_SIG_MAX_NODES:-6}
THETA=${THETA:-1.0}
NULL=${NULL:-perm}
SUMMARY=${SUMMARY:-1}                        # 0 = skip the summary (e.g. one job per seed)

# ---- (A) the real cohort, built once ----------------------------------------
# Deterministic (seed 0). The lock (flock, from util-linux) keeps concurrent
# jobs from rebuilding it on top of each other; a job that finds it already
# built skips the step.
echo "=== (A) prepare the matched real cohort ==="
mkdir -p "$OUT_BASE"
(
  flock 9
  if [ -f "${OUT_BASE}/real/inputs/graphs_all.txt" ] \
     && [ -n "$(ls "${OUT_BASE}/real/dags" 2>/dev/null)" ]; then
    echo "[SKIP] real cohort already built"
  else
    python3 breastcancer_to_pottr.py \
      --npy "$NPY" --out "${OUT_BASE}/real/dags" --multitree_out "${OUT_BASE}/real/inputs" \
      --max_trees "$MAX_TREES" --seed 0
  fi
) 9>"${OUT_BASE}/.real.lock"

if [ -n "$POTTR_ENV" ] && command -v conda >/dev/null 2>&1; then
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"
  if conda env list | awk '{print $1}' | grep -qx "$POTTR_ENV"; then
    conda activate "$POTTR_ENV"
  else
    echo "[warn] conda env '$POTTR_ENV' not found; using current python3"
  fi
fi
( cd "${POTTR_REPO}/code" && python3 -c "import gurobipy, networkx, pandas" ) \
  || { echo "[err] POTTR deps missing (gurobipy/networkx/pandas) in python3."; exit 1; }

KMIN=$(echo $K_LIST | tr ' ' '\n' | sort -n | head -1)
KMAX=$(echo $K_LIST | tr ' ' '\n' | sort -n | tail -1)

for SEED in $SEEDS; do
  SD="${OUT_BASE}/seed${SEED}"
  echo "=== seed ${SEED} ==="
  python3 permute_pottr_cohort.py \
    --dags_in "${OUT_BASE}/real/dags" --ens_in "${OUT_BASE}/real/inputs" \
    --dags_out "${SD}/dags" --ens_out "${SD}/inputs" --seed "$SEED"

  DAGS_ABS="$(cd "${SD}/dags" && pwd)"
  POTTR_OUT_ABS="$(mkdir -p "${SD}/pottr" && cd "${SD}/pottr" && pwd)"

  # ---- (B) POTTR k-sweep ----------------------------------------------------
  for K in $K_LIST; do
    KOUT="${POTTR_OUT_ABS}/k${K}"
    if [ -f "${KOUT}/converted_graphs.txt" ]; then
      echo "[SKIP] seed=${SEED} k=$K already has converted_graphs.txt"
      continue
    fi
    echo "--- POTTR seed=${SEED} k=$K ---"
    ( cd "${POTTR_REPO}/code" && \
      python3 run_POTTR.py -o "$KOUT" -d "$DAGS_ABS" -k "$K" -c "$CORES" -parallel )
  done

  # ---- (B2) POTTR's own p-value, per trajectory -----------------------------
  python3 pottr_force_significance.py \
    --pottr_dir "${SD}/pottr" --k_range "${KMIN},${KMAX}" \
    --dags "$DAGS_ABS" --pottr_repo "$POTTR_REPO" \
    --max_nodes "$POTTR_SIG_MAX_NODES" --cores "$CORES" \
    --out "${SD}/pottr/significance_forced.txt"

  # ---- (C) MAESTRI on the same trajectories ---------------------------------
  python3 pottr_significance.py \
    --pottr_dir "${SD}/pottr" --k_range "${KMIN},${KMAX}" \
    --pottr_sig_global "${SD}/pottr/significance_forced.txt" \
    --graphs_all "${SD}/inputs/graphs_all.txt" \
    -w "${SD}/inputs/weights_uniform.txt" --owner "${SD}/inputs/owner.txt" \
    --theta "$THETA" --null "$NULL" --n_jobs "$CORES" --seed "$SEED" \
    --out "${SD}/significance.csv"
done

# ---- (D) false discoveries per seed -----------------------------------------
[ "$SUMMARY" = 1 ] || exit 0
echo "=== (D) summary ==="
python3 summarize_pottr_null.py --csv "${OUT_BASE}"/seed*/significance.csv \
  --out "${OUT_BASE}/summary.csv"
