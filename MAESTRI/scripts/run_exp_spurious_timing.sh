#!/bin/bash
# Mining overhead on breastCancer: how many itemsets LCM enumerates that are
# not relation-complete trajectories, and the time of each stage
# (analyze_spurious_itemsets.py). Timings depend on the hardware; record the
# machine with any number taken from here.
set -euo pipefail
# Every path below is relative to the code directory (<repo>/MAESTRI), which
# holds this scripts/ folder; the datasets are in <repo>/data (../data).
cd "$(dirname "$0")/.."

# 5e6 itemset cap and 180 s timeout per LCM invocation; weighted mode, seed 0.
SIGMAS=${SIGMAS:-2,3,4,5,6,8,10,15,20}
MAXITEMS=${MAXITEMS:-5000000}
TIMEOUT=${TIMEOUT:-180}

python3 analyze_spurious_itemsets.py \
  --npy ../data/breastCancer.npy \
  --sigmas "$SIGMAS" \
  --seed 0 \
  --mode weighted \
  --max_itemsets "$MAXITEMS" \
  --timeout "$TIMEOUT" \
  --outdir results/spurious_analysis

echo "=== spurious timing complete -> results/spurious_analysis ==="
