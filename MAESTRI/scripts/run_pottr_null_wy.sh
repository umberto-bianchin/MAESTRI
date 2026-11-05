#!/bin/bash
# Westfall-Young thresholds of the null cohorts of
# run_pottr_null_breastcancer.sh (same per-patient null, resampled from each
# null cohort itself), and the count of MAESTRI p-values below them. Needs
# seed<S>/inputs and seed<S>/significance.csv from that driver.
# The 'perm' null of a null cohort is again a uniform per-patient relabeling of
# the real cohort, so every seed shares one null distribution: the thresholds
# are computed on the first seed of SEEDS and applied to all seeds (--wy_file).
set -euo pipefail
cd "$(dirname "$0")/.."                     # -> <repo>/MAESTRI, the code directory

OUT_BASE=${OUT_BASE:-results/pottr_null}
SEEDS=${SEEDS:-"1"}
SIGMA=${SIGMA:-2}
THETA=${THETA:-1.0}                          # the THETA of the null-cohort run
# SIGMA=2: every POTTR trajectory (k >= 2) is in the expected-support family
NULL=${NULL:-perm}
M=${M:-2000}                                 # WY resamples, as in the discovery runs
PAR=${PAR:-20}
MC_SAMPLES=${MC_SAMPLES:-20000}              # unused while every patient has <= MC_CUTOFF trees
MC_CUTOFF=${MC_CUTOFF:-8}
WY_NAME="wy_sigma${SIGMA}_theta${THETA}"

for SEED in $SEEDS; do
  SD="${OUT_BASE}/seed${SEED}"
  [ -f "${SD}/inputs/graphs_all.txt" ] && [ -f "${SD}/significance.csv" ] \
    || { echo "[err] ${SD}: run run_pottr_null_breastcancer.sh first"; exit 1; }

  echo "=== WY seed=${SEED} sigma=${SIGMA} theta=${THETA} M=${M} par=${PAR} ==="
  python3 run_wy_correction_multitree.py \
    --graphs_all "${SD}/inputs/graphs_all.txt" \
    -w           "${SD}/inputs/weights_uniform.txt" \
    --owner      "${SD}/inputs/owner.txt" \
    --sigma "$SIGMA" -M "$M" \
    --test both --null "$NULL" --theta "$THETA" --par "$PAR" \
    --min_mine_sigma 2 \
    --mc_cutoff "$MC_CUTOFF" --mc_samples "$MC_SAMPLES" \
    --outdir "${SD}/${WY_NAME}"
done

WY_SEED=${SEEDS%% *}                          # first seed: its thresholds serve all seeds
echo "=== summary with WY of seed ${WY_SEED} ==="
python3 summarize_pottr_null.py --csv "${OUT_BASE}"/seed*/significance.csv \
  --wy_file "${OUT_BASE}/seed${WY_SEED}/${WY_NAME}/wy_thresholds.txt" \
  --out "${OUT_BASE}/summary.csv"
