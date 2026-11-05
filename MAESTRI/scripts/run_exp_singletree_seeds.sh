#!/bin/bash
# Seed sweep of the single-tree baseline: for each seed, sample one tree per
# patient, mine and score the baseline family (run_pipeline.py
# --single_tree_only), correct it with Westfall-Young, then aggregate across
# seeds (analyze_singletree_seeds.py). DATASET=tracerx for TRACERx.
# Defaults are the paper configuration.
set -euo pipefail
# Every path below is relative to the code directory (<repo>/MAESTRI), which
# holds this scripts/ folder; the datasets are in <repo>/data (../data).
cd "$(dirname "$0")/.."

SEEDS=${SEEDS:-"0 1 2 3 4 5 6 7 8 9"}
SIGMA_LIST=${SIGMA_LIST:-"2 5"}
M=${M:-2000}                 # WY resamples
PAR=${PAR:-4}
NULL=${NULL:-perm}
MC_SAMPLES=${MC_SAMPLES:-3000}
MC_CUTOFF=${MC_CUTOFF:-8}
ALPHA=${ALPHA:-0.05}
# breastCancer (the npy) or tracerx (graphs + owner, as the TRACERx experiments)
DATASET=${DATASET:-breastCancer}
case "$DATASET" in
  breastCancer)
    NPY=../data/breastCancer.npy
    [ -f "$NPY" ] || { echo "ERROR: $NPY not found" >&2; exit 1; }
    INPUT=(--npy "$NPY")
    ROOT=${ROOT:-results/singletree_seeds} ;;
  tracerx)
    INPUT=(--graphs ../data/TRACERx/graphs_tracerx.txt --owner ../data/TRACERx/owner_tracerx.txt)
    ROOT=${ROOT:-results/singletree_seeds_tracerx} ;;
  *) echo "unknown DATASET $DATASET" >&2; exit 1 ;;
esac

for SIGMA in $SIGMA_LIST; do
  for SEED in $SEEDS; do
    OUT="$ROOT/sigma${SIGMA}/seed${SEED}"

    if [ -f "$OUT/significance/alg0_mastro_random_pvalues_exp.csv" ]; then
      echo "=== [SKIP] mining sigma=$SIGMA seed=$SEED (exists) ==="
    else
      echo "=== [1] single-tree mining + significance: sigma=$SIGMA seed=$SEED ==="
      python3 run_pipeline.py \
        "${INPUT[@]}" --sigma "$SIGMA" --seed "$SEED" \
        --single_tree_only \
        --significance --sig_null "$NULL" \
        --sig_mc_cutoff "$MC_CUTOFF" --sig_mc_samples "$MC_SAMPLES" --sig_n_jobs "$PAR" \
        --outdir "$OUT"
    fi

    PVAL_ALG0="$OUT/significance/alg0_mastro_random_pvalues_exp.csv"
    WYOUT="$ROOT/sigma${SIGMA}/seed${SEED}_wy"

    if [ -f "$WYOUT/wy_thresholds.txt" ]; then
      echo "=== [SKIP] WY sigma=$SIGMA seed=$SEED (exists) ==="
    else
      echo "=== [2] WY on single-tree family: sigma=$SIGMA seed=$SEED ==="
      python3 run_wy_correction_multitree.py \
        --graphs_all "$OUT/inputs/graphs_sampled.txt" \
        -w           "$OUT/inputs/weights_sampled.txt" \
        --owner      "$OUT/inputs/owner_sampled.txt" \
        --sigma      "$SIGMA" -M "$M" \
        --test exp --null "$NULL" --par "$PAR" \
        --mc_cutoff "$MC_CUTOFF" --mc_samples "$MC_SAMPLES" \
        --outdir "$WYOUT" \
        --pvalues_exp "$PVAL_ALG0"
    fi
  done
done

echo "=== [3] aggregate across seeds ==="
python3 analyze_singletree_seeds.py \
  --root "$ROOT" --sigma_list "$SIGMA_LIST" --alpha "$ALPHA"

echo "=== single-tree seed sweep complete -> $ROOT/summary.txt ==="
