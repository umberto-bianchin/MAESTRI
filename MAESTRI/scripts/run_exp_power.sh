#!/bin/bash
# Recovery power: implant a known trajectory into N carrier patients, in a
# fraction f of each carrier's trees, and measure how often each test recovers
# it at a Westfall-Young threshold (implant_experiment_multitree.py). Fully
# synthetic. Defaults are the paper configuration.
set -euo pipefail
# Every path below is relative to the code directory (<repo>/MAESTRI), which
# holds this scripts/ folder; the datasets are in <repo>/data (../data).
cd "$(dirname "$0")/.."

PAR=${PAR:-4}
THETA_LIST=${THETA_LIST:-"0.5 1.0"}
NULL=${NULL:-perm}
N_PATIENTS=${N_PATIENTS:-60}
MTREES=${MTREES:-10}
K=${K:-4}
N_LIST=${N_LIST:-2,3,5,10,20,30}
F_LIST=${F_LIST:-0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
SIGMA=${SIGMA:-2.0}
ALPHA=${ALPHA:-0.05}
N_TRIALS=${N_TRIALS:-200}
WY_M=${WY_M:-500}
CALIBRATION=${CALIBRATION:-cell}  # global: one carrier-free cohort; cell: one implanted cohort per (N, f)
# recall.csv, the WY cache and the trial scratch all go here.
OUT_BASE=${OUT_BASE:-results/power_v3_cell}
MC_CUTOFF=${MC_CUTOFF:-10}       # M_i above which a patient uses MC
MC_SAMPLES=${MC_SAMPLES:-2000}   # per-patient MC draws

for THETA in $THETA_LIST; do
  echo "=== power sweep: theta=$THETA null=$NULL calibration=$CALIBRATION wy_M=$WY_M ==="
  python3 implant_experiment_multitree.py \
    --n_patients "$N_PATIENTS" --M "$MTREES" --k "$K" \
    --N_list "$N_LIST" --f_list "$F_LIST" \
    --theta "$THETA" --sigma "$SIGMA" --alpha "$ALPHA" \
    --null "$NULL" --test both \
    --n_trials "$N_TRIALS" --wy_M "$WY_M" --calibration "$CALIBRATION" \
    --mc_cutoff "$MC_CUTOFF" --mc_samples "$MC_SAMPLES" --par "$PAR" \
    --outdir "${OUT_BASE}/implant_theta${THETA}"
done
echo "=== power experiment complete -> ${OUT_BASE} ==="
