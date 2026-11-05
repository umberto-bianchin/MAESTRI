#!/bin/bash
# Discovery on TRACERx: same pipeline as run_exp_discovery_breastcancer.sh,
# with the cohort given as a transaction file plus an owner file
# (--graphs / --owner). Defaults are the paper configuration.
set -euo pipefail
# Every path below is relative to the code directory (<repo>/MAESTRI), which
# holds this scripts/ folder; the datasets are in <repo>/data (../data).
cd "$(dirname "$0")/.."

M=${M:-2000}
PAR=${PAR:-4}
SIGMA_LIST=${SIGMA_LIST:-"2 5"}
THETA_LIST=${THETA_LIST:-"0.5 1.0"}
NULL=${NULL:-perm}
MC_SAMPLES=${MC_SAMPLES:-20000}  # per-patient MC draws; dominant cost knob
MC_CUTOFF=${MC_CUTOFF:-8}        # M_i above which a patient uses MC
EXACT_WHEN_CHEAPER=${EXACT_WHEN_CHEAPER:-1}  # exact above the cutoff when placements <= MC_SAMPLES; 0 = always MC
# Floor on the mining threshold of the theta-consensus family: when
# floor(theta*sigma) drops below it, the family is mined at sigma instead.
# Mining at support 1 is intractable on TRACERx (one transaction carries 1006
# items), hence 2, as in the breastCancer script. Empty = no floor.
MIN_MINE_SIGMA=${MIN_MINE_SIGMA-2}
OUT_ROOT=${OUT_ROOT:-results/exact_tracerx}
PLOT_ROOT=${PLOT_ROOT:-$OUT_ROOT/fdr_plots}
GRAPHS=../data/TRACERx/graphs_tracerx.txt
OWNER=../data/TRACERx/owner_tracerx.txt
THETA_CSV=$(echo "$THETA_LIST" | tr ' ' ',')
WORKDIR=${WORKDIR:-}             # scratch for WY mining; empty = <wy outdir>/resamples
WORKDIR_ARG=()
[ -n "$WORKDIR" ] && WORKDIR_ARG=(--workdir "$WORKDIR")
EXACT_SIG=(); EXACT_WY=()
if [ "$EXACT_WHEN_CHEAPER" = 1 ]; then
  EXACT_SIG=(--sig_exact_when_cheaper); EXACT_WY=(--exact_when_cheaper)
fi
WEIGHTS=${WEIGHTS:-}             # per-tree weights file (run_pipeline.py --weights); empty = uniform
WEIGHTS_ARG=()
[ -n "$WEIGHTS" ] && WEIGHTS_ARG=(--weights "$WEIGHTS")
MIN_MINE_ARG=()
[ -n "$MIN_MINE_SIGMA" ] && MIN_MINE_ARG=(--min_mine_sigma "$MIN_MINE_SIGMA")

[ -f "$GRAPHS" ] || { echo "ERROR: $GRAPHS not found" >&2; exit 1; }
mkdir -p "$PLOT_ROOT" "$OUT_ROOT"

for SIGMA in $SIGMA_LIST; do
  OUT="$OUT_ROOT/tracerx_sigma${SIGMA}"

  echo "=== [1] observed mining + significance: sigma=$SIGMA ==="
  python3 run_pipeline.py \
    --graphs "$GRAPHS" --owner "$OWNER" --sigma "$SIGMA" --seed 0 \
    --theta_list "$THETA_CSV" \
    ${MIN_MINE_ARG[@]+"${MIN_MINE_ARG[@]}"} \
    --significance --sig_null "$NULL" \
    --sig_mc_cutoff "$MC_CUTOFF" --sig_mc_samples "$MC_SAMPLES" --sig_n_jobs "$PAR" \
    ${EXACT_SIG[@]+"${EXACT_SIG[@]}"} \
    ${WEIGHTS_ARG[@]+"${WEIGHTS_ARG[@]}"} \
    --outdir "$OUT"

  SIGDIR="$OUT/significance"
  PVAL_EXP="$SIGDIR/alg1_expected_uniform_pvalues_exp.csv"

  for THETA in $THETA_LIST; do
    PVAL_THETA="$SIGDIR/alg3_theta${THETA}_pvalues_theta.csv"
    WYOUT="$OUT_ROOT/tracerx_wy_sigma${SIGMA}_theta${THETA}"

    echo "=== [2] WY/FDR: sigma=$SIGMA theta=$THETA M=$M par=$PAR null=$NULL ==="
    python3 run_wy_correction_multitree.py \
      --graphs_all "$OUT/inputs/graphs_all.txt" \
      -w           "$OUT/inputs/weights.txt" \
      --owner      "$OUT/inputs/owner.txt" \
      --sigma      "$SIGMA" -M "$M" \
      --test both --null "$NULL" --theta "$THETA" --par "$PAR" \
      ${MIN_MINE_ARG[@]+"${MIN_MINE_ARG[@]}"} \
      --mc_cutoff "$MC_CUTOFF" --mc_samples "$MC_SAMPLES" \
      ${EXACT_WY[@]+"${EXACT_WY[@]}"} \
      ${WORKDIR_ARG[@]+"${WORKDIR_ARG[@]}"} \
      --outdir "$WYOUT" \
      --pvalues_exp   "$PVAL_EXP" \
      --pvalues_theta "$PVAL_THETA" \
      --save_resample_pvals

    echo "=== [3] plot FDR: sigma=$SIGMA theta=$THETA ==="
    python3 plot_fdr_v2.py \
      --fdr_exp   "$WYOUT/fdr_exp.csv" \
      --fdr_theta "$WYOUT/fdr_theta.csv" \
      --out "$PLOT_ROOT/tracerx_sigma${SIGMA}_theta${THETA}.pdf" \
      --title "TRACERx sigma=${SIGMA} theta=${THETA}"
  done
done
echo "=== discovery (TRACERx) complete -> $OUT_ROOT ==="
