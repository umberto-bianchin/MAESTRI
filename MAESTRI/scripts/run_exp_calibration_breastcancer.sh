#!/bin/bash
# Calibration on breastCancer: empirical FWER of the Westfall-Young thresholds
# of both tests, from a pool of global-null datasets split repeatedly into a
# calibration and a validation half (empirical_fwer_multitree.py). The heaviest
# experiment; defaults are the paper configuration.
set -euo pipefail
# Every path below is relative to the code directory (<repo>/MAESTRI), which
# holds this scripts/ folder; the datasets are in <repo>/data (../data).
cd "$(dirname "$0")/.."

PAR=${PAR:-4}
SIGMA_LIST=${SIGMA_LIST:-"2"}
THETA_LIST=${THETA_LIST:-"0.5 1.0"}
NULL=${NULL:-perm}
N_DATASETS=${N_DATASETS:-500}
M_CAL=${M_CAL:-250}
N_TRIALS=${N_TRIALS:-500}
ALPHAS=${ALPHAS:-0.01,0.05,0.1}
MC_SAMPLES=${MC_SAMPLES:-20000}  # per-patient MC draws; dominant cost knob
MC_CUTOFF=${MC_CUTOFF:-8}        # M_i above which a patient uses MC
OUT_ROOT=${OUT_ROOT:-results/calibration_exact}
EXACT_WHEN_CHEAPER=${EXACT_WHEN_CHEAPER:-1}  # exact above the cutoff when placements <= MC_SAMPLES; 0 = always MC
EXACT_ARG=()
[ "$EXACT_WHEN_CHEAPER" = 1 ] && EXACT_ARG=(--exact_when_cheaper)
WORKDIR=${WORKDIR:-}             # scratch for null-dataset mining; empty = <outdir>/resamples
WORKDIR_ARG=()
[ -n "$WORKDIR" ] && WORKDIR_ARG=(--workdir "$WORKDIR")
NPY=../data/breastCancer.npy
THETA_CSV=$(echo "$THETA_LIST" | tr ' ' ',')

[ -f "$NPY" ] || { echo "ERROR: $NPY not found" >&2; exit 1; }

for SIGMA in $SIGMA_LIST; do
  INP=$OUT_ROOT/breastCancer_inputs_sigma${SIGMA}
  echo "=== [1] build inputs: sigma=$SIGMA ==="
  # Only the input files are used below: --min_mine_sigma skips the theta
  # candidate mining at floor(theta*sigma) (sigma=1 for theta=0.5), which blows up
  python3 run_pipeline.py --npy "$NPY" --sigma "$SIGMA" --seed 0 \
    --theta_list "$THETA_CSV" --min_mine_sigma "$SIGMA" --outdir "$INP"

  for THETA in $THETA_LIST; do
    echo "=== [2] empirical FWER: sigma=$SIGMA theta=$THETA null=$NULL ==="
    python3 empirical_fwer_multitree.py \
      --graphs_all "$INP/inputs/graphs_all.txt" \
      -w           "$INP/inputs/weights_uniform.txt" \
      --owner      "$INP/inputs/owner.txt" \
      --sigma "$SIGMA" --theta "$THETA" --null "$NULL" --test both \
      --n_datasets "$N_DATASETS" --m "$M_CAL" --n_trials "$N_TRIALS" \
      --alphas "$ALPHAS" \
      --mc_cutoff "$MC_CUTOFF" --mc_samples "$MC_SAMPLES" --par "$PAR" \
      ${EXACT_ARG[@]+"${EXACT_ARG[@]}"} \
      ${WORKDIR_ARG[@]+"${WORKDIR_ARG[@]}"} \
      --outdir $OUT_ROOT/breastCancer_sigma${SIGMA}_theta${THETA}
  done
done
echo "=== calibration (breastCancer) complete -> $OUT_ROOT ==="
