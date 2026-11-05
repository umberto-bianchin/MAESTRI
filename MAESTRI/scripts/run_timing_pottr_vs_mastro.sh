#!/bin/bash
# Wall-clock time of frequent-trajectory extraction, MAESTRI against POTTR, on
# matched breastCancer cohorts at several per-patient tree caps (or one prebuilt
# cohort, COHORT_DIR). No significance test on either side. The two methods run
# one after the other on the same trees; timeouts are recorded as data points.
# Writes $OUTDIR/timing.csv and environment.txt.
set -euo pipefail

# Every path below is relative to the code directory (<repo>/MAESTRI), which
# holds this scripts/ folder. COHORT_DIR is still resolved against the
# directory the script was invoked from, so remember it before moving.
INVOKED_FROM="$(pwd)"
cd "$(dirname "$0")/.."
CODE_DIR="$(pwd)"

NPY=${NPY:-../data/breastCancer.npy}
# Interpreter for every Python stage. Override when the PATH python3 is
# not the one carrying numpy / the POTTR deps.
PY=${PY:-python3}
# POTTR is cloned next to this repository: <parent>/MAESTRI and <parent>/POTTR
POTTR_REPO=${POTTR_REPO:-../../POTTR}
[ -f "${POTTR_REPO}/code/run_POTTR.py" ] || {
  echo "[err] no POTTR clone at ${POTTR_REPO} (clone it next to the repo, or set POTTR_REPO)."; exit 1; }
POTTR_REPO="$(cd "$POTTR_REPO" && pwd)"
echo "[info] POTTR repo: $POTTR_REPO"
POTTR_ENV=${POTTR_ENV:-pottr_env}
OUTDIR=${OUTDIR:-results/timing}

# Cohorts to sweep. An entry is either
#   <N>   cap: every patient keeps at most N distinct trees, heavy patients
#         truncated, cohort labeled cap<N>. 0 means no cap at all, i.e. the
#         full 37809-tree cohort
#   d<N>  drop: patients with more than N distinct trees are removed from the
#         cohort entirely, labeled drop<N>.
CAP_LIST=${CAP_LIST:-"1 2 10 d5"}
# A cohort that already exists, instead of one built from breastCancer.npy.
# Point it at a directory holding inputs/ (graphs_all.txt, owner.txt,
# weights_uniform.txt) and dags/ (one file per tree, named <patient>-<tree>).
# CAP_LIST is then ignored: there is one cohort and it is taken as given.
# Both sides must describe the same trees.
COHORT_DIR=${COHORT_DIR:-}
COHORT_NAME=${COHORT_NAME:-prebuilt}

# COHORT_DIR is resolved against the directory the script was invoked from, not
# the code directory it moves into.
if [ -n "$COHORT_DIR" ]; then
  case "$COHORT_DIR" in
    /*) ;;
    *) if [ -d "${INVOKED_FROM}/${COHORT_DIR}" ]; then
         COHORT_DIR="${INVOKED_FROM}/${COHORT_DIR}"
       fi ;;
  esac
  for sub in inputs dags; do
    [ -d "${COHORT_DIR}/${sub}" ] || {
      echo "[err] COHORT_DIR=${COHORT_DIR} has no ${sub}/ subdirectory."
      echo "      It must hold inputs/ (graphs_all.txt, owner.txt, weights_uniform.txt)"
      echo "      and dags/ (one file per tree, named <patient>-<tree>)."
      exit 1; }
  done
  COHORT_DIR="$(cd "$COHORT_DIR" && pwd)"
  # so that the run is labeled by the cohort actually used
  CAP_LIST="$COHORT_NAME"
  echo "[info] prebuilt cohort: $COHORT_DIR"
fi

# POTTR solves one ILP per recurrence level k, while MAESTRI returns the whole
# family for a given sigma in a single pass.
K_LIST=${K_LIST:-"2 5 10 20 50"}

SIGMA_LIST=${SIGMA_LIST:-"2 5"}
CORES=${CORES:-30}
TIMEOUT=${TIMEOUT:-14400}
SEED=${SEED:-0}
# Run only one half. SKIP_POTTR=1 also skips the Gurobi dependency check, so
# the MAESTRI side can be exercised on a machine without a Gurobi license.
SKIP_POTTR=${SKIP_POTTR:-0}
SKIP_MASTRO=${SKIP_MASTRO:-0}
# Abandon the k-sweep of a cohort after this many consecutive timeouts; 0 = never
MAX_CONSEC_TIMEOUTS=${MAX_CONSEC_TIMEOUTS:-0}
# Thread settings to time POTTR at (MAESTRI's mining is serial)
POTTR_THREADS=${POTTR_THREADS:-"$CORES"}
# Gurobi solution pool size; 0 (run_POTTR's default) returns one optimal solution
POOL=${POOL:-0}
# Pass -v to POTTR so Gurobi prints its MIP log
POTTR_VERBOSE=${POTTR_VERBOSE:-1}

mkdir -p "$OUTDIR"
# mastro_timing_run.py calls LCM directly: build it now if missing or foreign
"$PY" -c 'from pathlib import Path; import utils; utils.ensure_lcm(Path("lcm53"))'
# timing.csv is rebuilt from the per-run JSON records, so a resume never
# duplicates rows.
CSV="${OUTDIR}/timing.csv"

# ---- POTTR environment ------------------------------------------------------
if [ -n "$POTTR_ENV" ] && command -v conda >/dev/null 2>&1; then
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"
  if conda env list | awk '{print $1}' | grep -qx "$POTTR_ENV"; then
    conda activate "$POTTR_ENV"
  else
    echo "[warn] conda env '$POTTR_ENV' not found; using current python3"
  fi
fi
if [ "$SKIP_POTTR" != "1" ]; then
  ( cd "${POTTR_REPO}/code" && "$PY" -c "import gurobipy, networkx, pandas" ) \
    || { echo "[err] POTTR deps missing (gurobipy/networkx/pandas)."; \
         echo "      Set POTTR_ENV=<conda env> or install them for python3."; exit 1; }
fi

# Record the machine the numbers were measured on.
{
  echo "date          : $(date -Iseconds)"
  echo "host          : $(hostname)"
  echo "kernel        : $(uname -srm)"
  echo "cpu model     : $(grep -m1 'model name' /proc/cpuinfo 2>/dev/null | cut -d: -f2- | xargs || sysctl -n machdep.cpu.brand_string 2>/dev/null || echo unknown)"
  echo "cpu cores     : $(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo unknown)"
  echo "memory        : $(free -h 2>/dev/null | awk '/^Mem:/{print $2}' || echo unknown)"
  echo "python        : $("$PY" --version 2>&1)"
  echo "gurobi        : $("$PY" -c 'import gurobipy; print(gurobipy.gurobi.version())' 2>/dev/null || echo unknown)"
  echo "cores knob    : $CORES"
  echo "timeout       : ${TIMEOUT}s"
  if [ -n "$COHORT_DIR" ]; then
    echo "cohort        : $COHORT_NAME (prebuilt, $COHORT_DIR)"
  else
    echo "cohort list   : $CAP_LIST"
  fi
  echo "k list        : $(echo $K_LIST | tr ' ' '\n' | head -1)..$(echo $K_LIST | tr ' ' '\n' | tail -1)"
  echo "sigma list    : $SIGMA_LIST"
} > "${OUTDIR}/environment.txt"
cat "${OUTDIR}/environment.txt"

# ---- sweep ------------------------------------------------------------------
for CAP in $CAP_LIST; do
  echo
  if [ -n "$COHORT_DIR" ]; then
    echo "############ cohort: ${COHORT_NAME} (prebuilt) ############"
    COHORT_LABEL="${COHORT_NAME}"
    CDIR="${OUTDIR}/${COHORT_NAME}"
    ENS_DIR="$(cd "$COHORT_DIR/inputs" && pwd)"
    DAGS_DIR="$(cd "$COHORT_DIR/dags" && pwd)"
    mkdir -p "$CDIR"
  else
    # d<N> selects the drop arm, a bare number the cap arm
    case "$CAP" in
      d*) DROP_ABOVE="${CAP#d}"; MT=1000000
          COHORT_LABEL="drop${DROP_ABOVE}"
          echo "############ cohort: drop patients with > ${DROP_ABOVE} trees ############" ;;
      *)  DROP_ABOVE=0; MT="$CAP"
          # max_trees 0 means "no cap"; the converter takes a large number
          [ "$CAP" = "0" ] && MT=1000000
          COHORT_LABEL="cap${CAP}"
          echo "############ cohort: max_trees=${CAP} ############" ;;
    esac
    CDIR="${OUTDIR}/${COHORT_LABEL}"
    ENS_DIR="${CDIR}/inputs"
    DAGS_DIR="${CDIR}/dags"
    mkdir -p "$CDIR"
  fi

  # --- (A) matched cohort: same trees to both methods ---
  if [ -n "$COHORT_DIR" ]; then
    echo "=== (A) using prebuilt cohort at $COHORT_DIR ==="
  elif [ -f "${ENS_DIR}/manifest.csv" ]; then
    echo "[SKIP] cohort ${COHORT_LABEL} already built"
  else
    echo "=== (A) build matched cohort: ${COHORT_LABEL} ==="
    "$PY" breastcancer_to_pottr.py \
      --npy "$NPY" --out "$DAGS_DIR" --multitree_out "$ENS_DIR" \
      --max_trees "$MT" --drop_above "$DROP_ABOVE" \
      --n_patients 0 --seed "$SEED"
  fi

  if [ -f "${ENS_DIR}/manifest.csv" ]; then
    N_PAT=$(awk -F, 'NR>1{print $1}' "${ENS_DIR}/manifest.csv" | sort -u | wc -l | tr -d ' ')
  else
    N_PAT=$(sort -u "${ENS_DIR}/owner.txt" | grep -cve '^[[:space:]]*$' | tr -d ' ')
  fi
  N_TREE=$(grep -cve '^[[:space:]]*$' "${ENS_DIR}/graphs_all.txt" | tr -d ' ')
  echo "    cohort: ${N_PAT} patients, ${N_TREE} trees"

  # --- (B) MAESTRI, one pass per sigma, serial ---
  for S in $SIGMA_LIST; do
    [ "$SKIP_MASTRO" = "1" ] && break
    J="${CDIR}/mastro_sigma${S}.json"
    if [ -f "$J" ]; then
      echo "[SKIP] mastro sigma=${S}"
    else
      echo "=== (B) MAESTRI mining, sigma=${S} ==="
      rm -rf "${CDIR}/mastro_work_s${S}"
      if timeout "$TIMEOUT" "$PY" mastro_timing_run.py \
           --inputs "$ENS_DIR" --sigma "$S" \
           --workdir "${CDIR}/mastro_work_s${S}" \
           --label "${COHORT_LABEL}" --out "$J"; then
        :
      else
        echo "{\"method\":\"multi-mastro\",\"sigma\":${S},\"completed\":false,\"status\":\"timeout_or_error\",\"total_s\":${TIMEOUT}}" > "$J"
        echo "    [TIMEOUT/ERR] mastro sigma=${S}"
      fi
    fi
    # Refresh the CSV after every cell, so it is current even if the run is killed.
    "$PY" "${CODE_DIR}/rebuild_timing_csv.py" --outdir "$OUTDIR" >/dev/null
  done

  # --- (C) POTTR, one ILP per k, at every thread setting of POTTR_THREADS ---
  for TH in $POTTR_THREADS; do
    [ "$SKIP_POTTR" = "1" ] && break
    CONSEC=0
    PFLAG=""; [ "$TH" != "1" ] && PFLAG="--parallel"
    for K in $K_LIST; do
      J="${CDIR}/pottr_t${TH}_k${K}.json"
      if [ -f "$J" ]; then
        echo "[SKIP] pottr threads=${TH} k=${K}"
      else
        echo "=== (C) POTTR k=${K}, threads=${TH} ==="
        KOUT="${CDIR}/pottr_t${TH}/k${K}"
        LOG="${CDIR}/pottr_t${TH}_k${K}.log"
        mkdir -p "$KOUT"
        VFLAG=""; [ "$POTTR_VERBOSE" = "1" ] && VFLAG="--pottr_verbose"
        if timeout "$TIMEOUT" "$PY" "${CODE_DIR}/pottr_timing_run.py" \
             --pottr_repo "$POTTR_REPO" --dags "$(cd "$DAGS_DIR" && pwd)" \
             -k "$K" -c "$TH" --pool "$POOL" $PFLAG $VFLAG \
             --output_path "$KOUT" --label "${COHORT_LABEL}" --out "$J" 2>&1 | tee "$LOG"; then
          :
        else
          # pottr_timing_run.py writes a partial record on SIGTERM; write a
          # bare stub only when it left nothing.
          if ! "$PY" -c "import json,sys; d=json.load(open('$J')); sys.exit(0 if d.get('stages') else 1)" 2>/dev/null; then
            echo "{\"method\":\"pottr\",\"k\":${K},\"threads\":${TH},\"completed\":false,\"status\":\"timeout_no_partial\",\"total_s\":${TIMEOUT},\"stages\":[]}" > "$J"
          fi
          echo "    [TIMEOUT/ERR] pottr k=${K} threads=${TH} (log: $LOG)"
        fi
      fi
      "$PY" "${CODE_DIR}/rebuild_timing_csv.py" --outdir "$OUTDIR" >/dev/null
      if "$PY" -c "import json,sys; sys.exit(0 if json.load(open('$J')).get('completed') else 1)"; then
        CONSEC=0
      else
        CONSEC=$((CONSEC + 1))
        if [ "$MAX_CONSEC_TIMEOUTS" != "0" ] && [ "$CONSEC" -ge "$MAX_CONSEC_TIMEOUTS" ]; then
          echo "    [STOP] ${CONSEC} consecutive timeouts; abandoning the k-sweep (${COHORT_LABEL}, threads=${TH})"
          for KR in $K_LIST; do
            [ "$KR" -le "$K" ] && continue
            [ -f "${CDIR}/pottr_t${TH}_k${KR}.json" ] && continue
            echo "{\"method\":\"pottr\",\"k\":${KR},\"threads\":${TH},\"completed\":false,\"status\":\"not_attempted_after_timeouts\",\"total_s\":null,\"stages\":[]}" \
              > "${CDIR}/pottr_t${TH}_k${KR}.json"
          done
          break
        fi
      fi
    done
  done
done

echo
echo "=== done ==="
"$PY" "${CODE_DIR}/rebuild_timing_csv.py" --outdir "$OUTDIR"
echo "csv         : $CSV"
echo "environment : ${OUTDIR}/environment.txt"
# Pretty-print when column(1) is available; the csv is complete without it.
if command -v column >/dev/null 2>&1; then
  column -s, -t "$CSV" | cut -c1-140 | head -40
fi
