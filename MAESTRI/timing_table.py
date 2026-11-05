"""Collect per-cohort timing runs (one subdirectory per cohort, each with the
timing.csv and environment.txt of a timing run) into one CSV, a Markdown
summary and LaTeX table rows (MAESTRI in seconds, POTTR in minutes).

Usage:
    python3 timing_table.py --indir results/timing_dei
"""

import argparse
import csv
from pathlib import Path

from utils import SCRIPT_DIR

COHORTS = ["cap1", "cap2", "cap10", "drop5", "tracerx"]
LABELS = {"cap1": "cap $1$", "cap2": "cap $2$", "cap10": "cap $10$",
          "drop5": "drop $5$", "tracerx": "TRACERx"}
SIGMAS = ["2.0", "5.0"]
KS = ["2", "5", "10", "20", "50"]


def secs(x):
    return f"{x:.3f}" if x < 1 else f"{x:.1f}"


def load(indir):
    rows = {}
    for c in COHORTS:
        f = indir / c / "timing.csv"
        if f.exists():
            rows[c] = [dict(r, cohort=c) for r in csv.DictReader(open(f))]
    return rows


def cell(rows, method, threads, param, value):
    for r in rows:
        if (r["method"], r["threads"], r["param"], r["value"]) == (method, threads, param, value):
            return r
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", default=str(SCRIPT_DIR / "results/timing_dei"))
    ap.add_argument("--pottr_threads", default="30")
    args = ap.parse_args()
    indir = Path(args.indir)
    data = load(indir)
    th = args.pottr_threads

    merged = [r for c in COHORTS for r in data.get(c, [])]
    with open(indir / "timing_all.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(merged[0]))
        w.writeheader()
        w.writerows(merged)

    md = [f"# Timing, {indir.name}", "",
          f"MAESTRI serial (s); POTTR at {th} threads (min), no solution pool.", "",
          "| cohort | patients | trees | MAESTRI s=2 | MAESTRI s=5 | " + " | ".join(f"k={k}" for k in KS) + " |",
          "|---" * (5 + len(KS)) + "|"]
    tex = []
    for c in COHORTS:
        rows = data.get(c)
        if not rows:
            continue
        mm = [cell(rows, "multi-mastro", "1", "sigma", s) for s in SIGMAS]
        po = [cell(rows, "pottr", th, "k", k) for k in KS]

        def fm(r):
            return "--" if r is None else secs(float(r["total_s"])) + ("" if r["status"] == "ok" else "*")

        def fp(r):
            return "--" if r is None else f"{float(r['total_s']) / 60:.1f}" + ("" if r["status"] == "ok" else "*")
        n_pat, n_tree = rows[0]["n_patients"], rows[0]["n_trees"]
        md.append(f"| {c} | {n_pat} | {n_tree} | " + " | ".join(fm(r) for r in mm) + " | "
                  + " | ".join(fp(r) for r in po) + " |")
        tex.append(f"    {LABELS[c]:<9} & ${n_pat}$ & ${n_tree}$ & ${fm(mm[0])}$ s & "
                   + " & ".join(f"${fp(r)}$ m" for r in po) + r" \\")
        if c == "tracerx":
            one = [cell(rows, "pottr", "1", "k", k) for k in KS]
            if all(one):
                md.append("| tracerx, POTTR 1 thread | | | | | " + " | ".join(fp(r) for r in one) + " |")
    md += ["", "`*` not completed (timeout). Totals over k = "
           + ", ".join(KS) + ":", ""]
    for c in COHORTS:
        rows = data.get(c)
        if not rows:
            continue
        for t in sorted({r["threads"] for r in rows if r["method"] == "pottr"}, key=int):
            tot = sum(float(r["total_s"]) for r in rows
                      if r["method"] == "pottr" and r["threads"] == t and r["value"] in KS)
            md.append(f"- {c}: POTTR {t} thread(s) {tot:.1f} s")
    hosts = sorted({l.split(":", 1)[1].strip() for c in data
                    for l in open(indir / c / "environment.txt") if l.startswith("cpu model")})
    md += ["", f"CPU: {', '.join(hosts)}"]
    (indir / "summary.md").write_text("\n".join(md) + "\n")
    (indir / "timing_table.tex").write_text("\n".join(tex) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
