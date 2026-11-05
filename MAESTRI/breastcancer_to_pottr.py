"""
Prepare a matched POTTR + MAESTRI cohort from a multi-tree-per-patient .npy
dataset (e.g. breastCancer.npy), so the two methods run on exactly the same
trees.

Why a dedicated cohort is needed
--------------------------------
POTTR reads a directory in which each file holds one transitively-closed DAG
(one tree). It groups files into evolutionary processes (patients) by the file
name ``<patient>-<tree>...`` and its ILP then selects at most one tree per
patient. Two things matter on breastCancer:

1. Scale. POTTR builds a pairwise conflict graph for every pair of trees
   across patients, O(T^2) in the total number of trees T. breastCancer has
   37,809 candidate trees, about 6.6e8 cross-patient pairs: the conflict-graph
   stage does not finish.

2. Support counts trees. ``compute_support.py`` tallies support by graph name
   (one per tree, ``P<i>-<j>``) and rescans all candidate trees, not just the
   one the ILP picked, so on a multi-tree cohort a patient counts once per
   candidate tree that carries the trajectory. On a single-tree cohort
   (AML/NSCLC) name and patient coincide.

To keep the run tractable and the two inputs identical, this script:
  * deduplicates identical trees within a patient (POTTR would count them
    separately, see point 2),
  * caps the number of candidate trees per patient (``--max_trees``) or drops
    the patients that exceed it (``--drop_above``),
  * optionally restricts to multi-tree patients (``--multitree_only``), where
    tree and patient counts differ, and/or to a random subset of
    ``--n_patients`` patients,
  * writes the POTTR dags (``P<i>-<j>_bc.txt``) and the matched MAESTRI inputs
    (``graphs_all.txt``, ``weights_uniform.txt``, ``owner.txt``,
    ``graphs_sampled.txt``) built from the same selected trees, plus a
    ``manifest.csv`` recording the selection.

The MAESTRI inputs feed ``pottr_significance.py``; because both sides describe
the same trees, occurrence matching is exact.

Usage (as scripts/run_pottr_breastcancer.sh calls it):
    python3 breastcancer_to_pottr.py \\
        --npy ../data/breastCancer.npy \\
        --out results/pottr_cmp/dags \\
        --multitree_out results/pottr_cmp/bc_inputs \\
        --max_trees 2 --n_patients 0 --seed 0
"""

import argparse
import csv
from pathlib import Path

import numpy as np

from utils import (
    build_inputs,
    extract_edges_from_tree,
    tree_edges_to_transaction_items,
)


def _tree_key(tr, drop_gl):
    """Canonical hashable key for a tree = its set of pairwise-relation items."""
    edges = extract_edges_from_tree(tr)
    return frozenset(tree_edges_to_transaction_items(edges, drop_gl=drop_gl))


def _dedup_trees(trees, drop_gl):
    """Drop trees with an identical pairwise-relation set, and trees with no
    pair at all (fewer than two mutations), keeping first occurrence order."""
    seen = set()
    out = []
    for tr in trees:
        k = _tree_key(tr, drop_gl)
        if not k or k in seen:
            continue
        seen.add(k)
        out.append(tr)
    return out


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--npy", required=True, help="Multi-tree .npy dataset")
    ap.add_argument("--out", required=True, help="Output dir of POTTR DAG files")
    ap.add_argument("--multitree_out", default="",
                    help="If set, also write the matched MAESTRI inputs "
                         "(graphs_all/weights/owner/sampled) here")
    ap.add_argument("--max_trees", type=int, default=5,
                    help="Cap on distinct candidate trees per patient after "
                         "dedup (default 5). Total trees drive POTTR's O(T^2) cost")
    ap.add_argument("--min_trees", type=int, default=1,
                    help="Keep only patients with >= this many distinct trees")
    ap.add_argument("--drop_above", type=int, default=0,
                    help="Drop patients with more than this many distinct "
                         "trees instead of truncating them (0 = keep all). "
                         "Unlike --max_trees this removes the patient from "
                         "the cohort, so no patient is reported on a subset "
                         "of its own candidate trees")
    ap.add_argument("--multitree_only", action="store_true",
                    help="Keep only patients with >=2 distinct trees (where "
                         "tree and patient counts differ)")
    ap.add_argument("--n_patients", type=int, default=0,
                    help="Random subset of eligible patients (0 = all eligible)")
    ap.add_argument("--keep_gl", action="store_true",
                    help="Do not drop the germline root GL")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    drop_gl = not args.keep_gl
    data = np.load(args.npy, allow_pickle=True)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    # purge stale dags from a previous selection: same-named files get
    # overwritten, but files for trees no longer selected would linger and
    # POTTR would read them as extra trees. Clear them first.
    for f in out.glob("*_bc.txt"):
        f.unlink()

    min_trees = max(args.min_trees, 2 if args.multitree_only else 1)

    # ---- select cohort once: dedup, cap, eligibility, subset ---------------
    selected = []  # list of (orig_patient_idx, [trees])
    n_dropped_above = 0
    n_pairless_trees = 0
    n_pairless_patients = 0
    for i, tlist in enumerate(data):
        n_pairless = sum(1 for tr in tlist if not _tree_key(tr, drop_gl))
        n_pairless_trees += n_pairless
        if len(tlist) > 0 and n_pairless == len(tlist):
            n_pairless_patients += 1
        trees = _dedup_trees(list(tlist), drop_gl)
        if len(trees) < min_trees:
            continue
        # --max_trees keeps heavy patients on part of their trees; --drop_above
        # removes them, so every kept patient has all its candidate trees
        if args.drop_above and len(trees) > args.drop_above:
            n_dropped_above += 1
            continue
        if len(trees) > args.max_trees:
            trees = trees[:args.max_trees]  # deterministic: first distinct trees
        selected.append((i, trees))

    if args.n_patients and args.n_patients < len(selected):
        idx = sorted(rng.choice(len(selected), size=args.n_patients,
                                replace=False).tolist())
        selected = [selected[j] for j in idx]

    # ---- write POTTR dags: one file per tree, named P<i>-<j>_bc.txt ---------
    n_files = 0
    tree_counts = []
    for pi, trees in selected:
        tree_counts.append(len(trees))
        for j, tr in enumerate(trees):
            edges = extract_edges_from_tree(tr)
            items = tree_edges_to_transaction_items(edges, drop_gl=drop_gl)
            (out / f"P{pi}-{j}_bc.txt").write_text(" ".join(items) + "\n")
            n_files += 1

    # ---- write the matched MAESTRI inputs from the same selected trees -----
    if args.multitree_out:
        ens = Path(args.multitree_out)
        ens.mkdir(parents=True, exist_ok=True)
        patients_trees = [trees for _, trees in selected]
        build_inputs(patients_trees, ens, seed=args.seed, drop_gl=drop_gl)
        with (ens / "manifest.csv").open("w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["ensemble_owner_id", "orig_patient_idx", "n_trees"])
            for pos, (pi, trees) in enumerate(selected):
                w.writerow([pos, pi, len(trees)])

    tc = np.array(tree_counts) if tree_counts else np.array([0])
    print(f"[bc->POTTR] patients selected = {len(selected)}, tree files = {n_files}")
    print(f"[bc->POTTR] dropped for having no mutation pair (< 2 mutations): "
          f"{n_pairless_trees} trees, {n_pairless_patients} patients left with no tree")
    if args.drop_above:
        print(f"[bc->POTTR] dropped for > {args.drop_above} distinct trees: "
              f"{n_dropped_above} patients")
    print(f"[bc->POTTR] distinct trees/patient: min={tc.min()} "
          f"median={int(np.median(tc))} max={tc.max()} "
          f"(multi-tree patients = {(tc > 1).sum()})")
    print(f"[bc->POTTR] total trees T = {int(tc.sum())} "
          f"-> POTTR pairwise cost ~ {int(tc.sum())**2 // 2} conflict graphs")
    print(f"[bc->POTTR] POTTR dags: {out}")
    if args.multitree_out:
        print(f"[bc->POTTR] matched MAESTRI inputs: {args.multitree_out}")


if __name__ == "__main__":
    main()
