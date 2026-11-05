"""Where the p-values of an ablation-toy sweep cross alpha.

Reads one ablation_toys.py CSV (columns toy, x, variant, pattern, pval_exp,
pval_theta) and, for every toy, variant, pattern and test (expected support,
and the theta test where the CSV has it), reports whether the p-value is
monotone in x, every x at which it crosses alpha (with the values on both
sides), and the p-values at the requested x. Written for the sweeps of
ablation_toys.py --sweep, e.g. --toy competition --sweep 1:200, where x is M.

Usage:
    python3 summarize_toy_sweep.py --csv results/ablation_toys/competition_sweep.csv \\
        --at 10,50,100,200
"""

import argparse
import csv
from pathlib import Path

DISPLAY_NAMES = {"maestri": "MAESTRI", "count_trees": "tree counts",
                 "perm_tree": "per-tree null", "trees_as_patients": "pooled trees",
                 "count_perm_tree": "tree counts + per-tree null"}


def monotone(p, tol=1e-12):
    """Up to float noise: exact laws summed to 1 - 1e-15 must not read as a
    change of direction."""
    if all(a <= b + tol for a, b in zip(p, p[1:])):
        return "constant" if abs(p[0] - p[-1]) <= tol else "non-decreasing"
    if all(a + tol >= b for a, b in zip(p, p[1:])):
        return "non-increasing"
    return "not monotone"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", required=True)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--at", default="10,50,100,200",
                    help="x values whose p-values are listed")
    ap.add_argument("--out", default=None, help="Markdown summary (default: <csv>.md)")
    args = ap.parse_args()

    curves = {}
    for r in csv.DictReader(open(args.csv)):
        for test, col in (("exp", "pval_exp"), ("theta", "pval_theta")):
            if r.get(col):
                curves.setdefault((r["toy"], r["variant"], r["pattern"], test), []).append(
                    (int(r["x"]), float(r[col])))
    at = [int(x) for x in args.at.split(",") if x]
    a = args.alpha

    md = [f"# {Path(args.csv).name}: crossings of p = {a}", "",
          f"Significant means p <= {a} (uncorrected p-value; test exp = expected "
          "support, theta = theta-consensus).", "",
          "| toy | variant | pattern | test | x range | at first x | crossings | monotone | "
          + " | ".join(f"x={x}" for x in at) + " |",
          "|---" * (8 + len(at)) + "|"]
    for (toy, v, pat, test), pts in curves.items():
        pts.sort()
        xs, ps = [x for x, _ in pts], [p for _, p in pts]
        sig = [p <= a for p in ps]
        cross = []
        for i in range(1, len(xs)):
            if sig[i] != sig[i - 1]:
                cross.append(f"{'significant' if sig[i] else 'not significant'} from x={xs[i]}"
                             f" ({ps[i - 1]:.3g} at {xs[i - 1]}, {ps[i]:.3g} at {xs[i]})")
        vals = {x: p for x, p in pts}
        md.append(f"| {toy} | {DISPLAY_NAMES.get(v, v)} | {pat} | {test} | {xs[0]}-{xs[-1]} | "
                  f"{'significant' if sig[0] else 'not significant'} | "
                  f"{'; '.join(cross) or 'none'} | {monotone(ps)} | "
                  + " | ".join(f"{vals[x]:.2g}" if x in vals else "--" for x in at) + " |")
    out = Path(args.out) if args.out else Path(args.csv).with_suffix(".md")
    out.write_text("\n".join(md) + "\n")
    print("\n".join(md))
    print(f"\n[out] {out}")


if __name__ == "__main__":
    main()
