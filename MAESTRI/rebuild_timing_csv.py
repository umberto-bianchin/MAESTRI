"""Rebuild timing.csv from the per-run JSON records.

The CSV is a function of the records on disk, so rebuilding it is idempotent
and stays correct after a resumed run (scripts/run_timing_pottr_vs_mastro.sh
calls this after every cell).

Usage:
    python3 rebuild_timing_csv.py --outdir results/timing
"""

import argparse
import csv
import json

from pottr_timing_run import read_trajectories
from pathlib import Path

from utils import SCRIPT_DIR

HEADER = ["cohort", "n_patients", "n_trees", "method", "threads",
          "param", "value", "total_s", "cpu_total_s", "cpu_per_wall", "status",
          "n_trajectories", "pottr_support", "pottr_patients", "detail"]


def pottr_output(cohort_dir: Path, threads, k):
    """What POTTR returned for one run, parsed from its output file with
    pottr_timing_run.read_trajectories (not taken from the JSON record).

    POTTR's support counts trees (occurrence ids are "P<patient>-<tree>"), so
    the distinct-patient count is reported next to it.
    """
    f = cohort_dir / f"pottr_t{threads}" / f"k{k}" / "converted_graphs.txt"
    found = read_trajectories(f)
    if not found["n_trajectories"]:
        return "", "", ""
    tr = found["trajectories"]
    return (found["n_trajectories"],
            ";".join(str(x["support"]) for x in tr),
            ";".join(str(x["n_patients"]) for x in tr))


def cohort_size(cohort_dir: Path):
    """(patients, trees) of a cohort.

    Prefers the transaction files (not tracked in the repository), then
    manifest.csv (one row per patient with its tree count).
    """
    owner = cohort_dir / "inputs" / "owner.txt"
    graphs = cohort_dir / "inputs" / "graphs_all.txt"
    if owner.exists() and graphs.exists():
        pats = len({l.strip() for l in open(owner) if l.strip()})
        trees = sum(1 for l in open(graphs) if l.strip())
        return pats, trees

    manifest = cohort_dir / "inputs" / "manifest.csv"
    if manifest.exists():
        rows = list(csv.DictReader(open(manifest)))
        return len(rows), sum(int(r["n_trees"]) for r in rows)

    # A prebuilt cohort keeps its inputs outside the run directory: fall back
    # to the MAESTRI record, which measured both while mining.
    for j in sorted(cohort_dir.glob("mastro_*.json")):
        try:
            rec = json.load(open(j))
        except (OSError, json.JSONDecodeError):
            continue
        if rec.get("n_patients") and rec.get("n_trees"):
            return rec["n_patients"], rec["n_trees"]
    return "", ""


def row_from(rec, cap, pats, trees, cohort_dir):
    stages = rec.get("stages", [])
    detail = ";".join(
        f"{s['stage']}={s['elapsed_s']}"
        + (f"@{s['cpu_per_wall']}x" if s.get("cpu_per_wall") else "")
        for s in stages)
    if rec.get("interrupted_in"):
        detail += f";KILLED_IN={rec['interrupted_in']}@{rec.get('interrupted_after_s')}s"
    status = rec.get("status") or ("ok" if rec.get("completed") else "timeout")
    # Interrupted records carry cpu_total_s but no ratio: derive it.
    cpw = rec.get("cpu_per_wall")
    if cpw in (None, "") and rec.get("cpu_total_s") and rec.get("total_s"):
        cpw = round(rec["cpu_total_s"] / rec["total_s"], 2)
    if rec["method"] == "pottr":
        n, sup, tpats = pottr_output(cohort_dir, rec.get("threads", ""), rec.get("k", ""))
        return [cap, pats, trees, "pottr", rec.get("threads", ""), "k",
                rec.get("k", ""), rec.get("total_s", ""),
                rec.get("cpu_total_s", ""), cpw if cpw is not None else "", status,
                n, sup, tpats, detail]
    fam = rec.get("n_trajectories", {})
    return [cap, pats, trees, "multi-mastro", rec.get("threads", 1), "sigma",
            rec.get("sigma", ""), rec.get("total_s", ""),
            rec.get("cpu_total_s", ""), cpw if cpw is not None else "", status,
            fam.get("exp", ""), "", "", detail]


def sort_key(r):
    # Numeric cohort labels are tree caps, sorted by size with 0 (no cap) last;
    # named cohorts sort after the numeric ones.
    label = str(r[0])
    if label.isdigit():
        cap = int(label)
        first = (0, cap if cap else 10 ** 9, "")
    else:
        first = (1, 0, label)
    return (first, r[3], int(r[4] or 0), float(r[6] or 0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=str(SCRIPT_DIR / "results/timing"))
    args = ap.parse_args()
    out = Path(args.outdir)

    rows = []
    # Every subdirectory is a cohort: cap<N> for a tree cap, otherwise the
    # name of a prebuilt cohort.
    for cdir in sorted(d for d in out.iterdir() if d.is_dir()):
        if not cdir.is_dir():
            continue
        cap = cdir.name[3:] if cdir.name.startswith("cap") else cdir.name
        pats, trees = cohort_size(cdir)
        for j in sorted(cdir.glob("*.json")):
            try:
                rec = json.load(open(j))
            except json.JSONDecodeError:
                print(f"[warn] unreadable, skipped: {j}")
                continue
            if "method" not in rec:
                continue
            rows.append(row_from(rec, cap, pats, trees, cdir))

    rows.sort(key=sort_key)
    csv_path = out / "timing.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(HEADER)
        w.writerows(rows)
    print(f"[OK] {csv_path}: {len(rows)} rows")
    for r in rows:
        print("  " + " | ".join(str(x) for x in r[:10]))


if __name__ == "__main__":
    main()
