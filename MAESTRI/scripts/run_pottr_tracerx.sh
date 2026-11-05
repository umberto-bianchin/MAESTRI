#!/bin/bash
# POTTR against MAESTRI on TRACERx, the counterpart of run_pottr_breastcancer.sh.
# No cap: both methods see all 469 trees, written from graphs_tracerx.txt and
# owner_tracerx.txt. '-?-' means same clone here, so POTTR's cluster pairs are
# kept (KEEP_CLUSTER=1). POTTR follows its paper's protocol: pool of 50000,
# k = 2..35 (beyond that only one-alteration trajectories). Needs Gurobi.
set -euo pipefail
cd "$(dirname "$0")/.."                     # -> <repo>/MAESTRI, the code directory

GRAPHS=${GRAPHS:-../data/TRACERx/graphs_tracerx.txt}
OWNER=${OWNER:-../data/TRACERx/owner_tracerx.txt}
CDIR=${CDIR:-results/pottr_cmp/tracerx_cohort}
SIGMA=${SIGMA:-2}

# ---- matched cohort ---------------------------------------------------------
mkdir -p "$CDIR/dags"
if [ -z "$(ls -A "$CDIR/dags")" ]; then
  awk -v d="$CDIR/dags" 'NR==FNR { own[FNR] = $0; next }
    { p = own[FNR]; f = d "/" p "-" (n[p]++) "_tracerx.txt"; print > f; close(f) }' \
    "$OWNER" "$GRAPHS"
fi
N_DAGS=$(ls "$CDIR/dags" | wc -l | tr -d ' ')
[ "$N_DAGS" = 469 ] || { echo "[err] $N_DAGS dags in $CDIR/dags, expected 469" >&2; exit 1; }

if [ ! -f "$CDIR/inputs/graphs_all.txt" ]; then
  python3 run_pipeline.py --graphs "$GRAPHS" --owner "$OWNER" --sigma "$SIGMA" --theta_list 1.0 \
    --seed 0 --outdir "$CDIR"
fi
# PREP_ONLY=1 stops here, for the experiments that only need the cohort (timing)
[ "${PREP_ONLY:-0}" = "1" ] && { echo "=== prep only: cohort in $CDIR ==="; exit 0; }

# ---- POTTR k-sweep, POTTR p-values, MAESTRI re-scoring -----------------------
SKIP_PREP=1 KEEP_CLUSTER=1 EXACT_WHEN_CHEAPER=1 CORES=${CORES:-20} \
POOL=${POOL:-50000} K_LIST=${K_LIST:-"$(seq 2 35 | tr '\n' ' ')"} \
MC_SAMPLES=${MC_SAMPLES:-20000} \
DAGS_DIR="$CDIR/dags" ENS_DIR="$CDIR/inputs" \
POTTR_OUT=${POTTR_OUT:-results/pottr_cmp/pottr_tracerx} \
SIG_OUT=${SIG_OUT:-results/pottr_cmp/pottr_tracerx_significance.csv} \
PLOT_OUT=${PLOT_OUT:-results/pottr_cmp/pottr_tracerx_comparison.pdf} \
  bash scripts/run_pottr_breastcancer.sh
