#!/bin/bash
# Ablations on the real cohorts: five variants of the test, each breaking one
# design choice of MAESTRI (see ablation_calibration.py), calibrated on their
# own null and validated on per-patient null datasets, plus Q-Q plots and the
# real-data discoveries of every variant. Defaults are the paper configuration.
set -euo pipefail
cd "$(dirname "$0")/.."                     # -> <repo>/MAESTRI, the code directory

DATASETS=${DATASETS:-"breastCancer tracerx"}
SIGMA_LIST=${SIGMA_LIST:-"5"}
N_CAL=${N_CAL:-500}
N_VAL=${N_VAL:-500}
PAR=${PAR:-50}
MC_SAMPLES=${MC_SAMPLES:-20000}  # per-patient MC draws, 'perm' variants only
MC_CUTOFF=${MC_CUTOFF:-8}
EXACT_WHEN_CHEAPER=${EXACT_WHEN_CHEAPER:-1}  # exact above the cutoff when placements <= MC_SAMPLES
ALPHAS=${ALPHAS:-0.01,0.05,0.1}
OUT_BASE=${OUT_BASE:-results/ablation_v3}
WORKDIR=${WORKDIR:-}             # mining scratch; empty = <outdir>/work
SHARD=${SHARD:-}                 # I/N: this job's share of the null datasets (no analysis)

for DS in $DATASETS; do
  for SIGMA in $SIGMA_LIST; do
    INP=${OUT_BASE}/${DS}_inputs_sigma${SIGMA}
    if [ -f "$INP/inputs/graphs_all.txt" ]; then
      echo "=== [SKIP] inputs ${DS} sigma=${SIGMA} (exist) ==="
    else
      echo "=== [1] build inputs: ${DS} sigma=${SIGMA} ==="
      case "$DS" in
        breastCancer)
          python3 run_pipeline.py --npy ../data/breastCancer.npy --sigma "$SIGMA" --theta_list 1.0 \
            --seed 0 --outdir "$INP" ;;
        tracerx)
          python3 run_pipeline.py --graphs ../data/TRACERx/graphs_tracerx.txt \
            --owner ../data/TRACERx/owner_tracerx.txt --sigma "$SIGMA" --theta_list 1.0 \
            --seed 0 --outdir "$INP" ;;
        *) echo "unknown dataset $DS" >&2; exit 1 ;;
      esac
    fi

    OUT=${OUT_BASE}/${DS}_sigma${SIGMA}
    echo "=== [2] ablation: ${DS} sigma=${SIGMA} n_cal=${N_CAL} n_val=${N_VAL} ==="
    # Resumable: datasets whose p-values are already on disk are not recomputed
    python3 ablation_calibration.py \
      --graphs_all "$INP/inputs/graphs_all.txt" \
      -w           "$INP/inputs/weights_uniform.txt" \
      --owner      "$INP/inputs/owner.txt" \
      --sigma "$SIGMA" --n_cal "$N_CAL" --n_val "$N_VAL" --alphas "$ALPHAS" \
      --mc_cutoff "$MC_CUTOFF" --mc_samples "$MC_SAMPLES" --par "$PAR" \
      --outdir "$OUT" ${WORKDIR:+--workdir "$WORKDIR/${DS}_sigma${SIGMA}"} \
      ${SHARD:+--shard "$SHARD"} \
      $( [ "$EXACT_WHEN_CHEAPER" = 1 ] && echo --exact_when_cheaper )

    [ -n "$SHARD" ] && continue
    python3 plot_ablation_qq.py --dir "$OUT"
    python3 plot_ablation_qq.py --dir "$OUT" --family mined
  done
done
echo "=== ablation complete -> ${OUT_BASE} ==="
