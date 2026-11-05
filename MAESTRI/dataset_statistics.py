"""Print and save summary statistics for a multi-tree dataset.

Loads a dataset, either a .npy (object array indexed
data[patient][tree][edge] = [parent, child], e.g. breastCancer) or one
transaction per tree plus its owner (--graphs/--owner, e.g. TRACERx), and
computes the statistics of utils.dataset_stats_from_mutsets (number of
patients, single- and multi-tree patients, trees-per-patient distribution,
mutation counts, fraction of patients whose trees all share the same mutation
set). It prints a short summary, optionally the trees-per-patient histogram,
and writes the full statistics to a JSON file.

    python3 dataset_statistics.py --npy ../data/breastCancer.npy --tree_dist --outdir OUT
    python3 dataset_statistics.py --graphs ../data/TRACERx/graphs_tracerx.txt \
        --owner ../data/TRACERx/owner_tracerx.txt --tree_dist --outdir OUT
"""

import argparse
import json
import re
from pathlib import Path
from datetime import datetime
from utils import (
    ensure_dir,
    load_npy,
    compute_dataset_stats,
    dataset_stats_from_mutsets,
)


def main():
    ap = argparse.ArgumentParser()
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--npy", help="Path to a .npy dataset (object array data[patient][tree][edge]=[parent,child])")
    src.add_argument("--graphs", help="One transaction per tree (e.g. graphs_tracerx.txt); needs --owner")
    ap.add_argument("--owner", help="Patient of each line of --graphs")
    ap.add_argument("--outdir", default=None, help="Output directory (default: results_<timestamp>)")
    ap.add_argument("--keep_gl", action="store_true", help="Keep GL in mutation sets")
    ap.add_argument("--tree_dist", action="store_true",
                    help="Also compute and print the distribution of #trees per patient")
    args = ap.parse_args()

    drop_gl = not args.keep_gl

    if args.outdir is None:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        outdir = Path(f"results_{ts}").resolve()
    else:
        outdir = Path(args.outdir).resolve()
    ensure_dir(outdir)

    if args.npy:
        patients_trees = list(load_npy(Path(args.npy).resolve()))
        stats = compute_dataset_stats(patients_trees, drop_gl=drop_gl)
    else:
        if not args.owner:
            ap.error("--graphs needs --owner")
        trees = {}
        # --graphs and --owner must have the same number of lines (zip truncates)
        for line, o in zip(open(args.graphs).read().splitlines(), open(args.owner).read().split()):
            # the separators of utils.EDGE_SEPARATORS
            muts = {x for it in line.split() for x in re.split(r"->-|-\?-|-/-", it)}
            if drop_gl:
                # a root labeled "0" is treated as germline too
                muts -= {"GL", "0"}
            trees.setdefault(o, []).append(muts)
        patients_trees = list(trees.values())
        stats = dataset_stats_from_mutsets(patients_trees, drop_gl)
    print(f"patients {stats['n_patients']}: single-tree {stats['single_tree_patients']}, "
          f"multi-tree {stats['multi_tree_patients']}; trees per patient median "
          f"{stats['trees_per_patient']['median']:g}, mean {stats['trees_per_patient']['mean']:.2f}, "
          f"max {stats['trees_per_patient']['max']}")

    if args.tree_dist:
        from collections import Counter
        dist = Counter(len(tlist) for tlist in patients_trees)
        stats["trees_per_patient_distribution"] = {
            str(k): v for k, v in sorted(dist.items())
        }
        print("\n--- Distribution: #trees -> #patients ---")
        for n_trees, n_patients in sorted(dist.items()):
            bar = "#" * min(n_patients, 60)
            print(f"  {n_trees:4d} trees: {n_patients:5d} patients  {bar}")

    (outdir / "dataset_stats.json").write_text(json.dumps(stats, indent=2))
    print("[OK] Dataset stats written to", outdir / "dataset_stats.json")

    print("\n[DONE] Outputs are under:", outdir)


if __name__ == "__main__":
    main()
