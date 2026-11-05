"""
Re-evaluate POTTR trajectories under the MAESTRI significance tests.

POTTR selects, via ILP, the maximum recurrent trajectory (an incomplete poset
in general, with clusters possibly resolved and one tree chosen per patient)
and assesses its significance with its adaptation of MASTRO's single-tree test.
This script re-scores each POTTR trajectory under the multi-tree tests, with
every candidate tree of a patient weighted 1/M_i (expected-support and
theta-consensus), on occurrences recomputed from the MAESTRI transactions.

For each POTTR trajectory we keep only its ordered and incomparable relations
(``a->-b`` and, if present, ``a-/-b``), dropping the germline root ``0`` and
the same-node ``-?-`` pairs. With ``--keep_cluster`` the ``-?-`` pairs are
kept instead: in the TRACERx input ``-?-`` means "same node"
(data/TRACERx/create_tracerx_trees.py), and a POTTR cluster pair it did not
resolve is still a same-node claim on every tree that supports the trajectory.
We then recompute the trajectory's occurrences directly in the MAESTRI
transaction database and run
``compute_significance_multitree.py --allow-incomplete`` on the resulting
(incomplete) pattern, so only the relations POTTR actually commits to constrain
the null.

Output: one row per POTTR trajectory with its size, recomputed MAESTRI support,
expected-support and theta-consensus p-values, and, when a POTTR
``significance_output.txt`` is available, POTTR's own permutation p-value for a
side-by-side comparison.

Usage:
    python3 pottr_significance.py \\
        --pottr_dir results/pottr_cmp/pottr_tracerx --k_range 2,35 \\
        --graphs_all results/pottr_cmp/tracerx_cohort/inputs/graphs_all.txt \\
        -w results/pottr_cmp/tracerx_cohort/inputs/weights_uniform.txt \\
        --owner results/pottr_cmp/tracerx_cohort/inputs/owner.txt \\
        --theta 1.0 --null perm --mc_samples 20000 --exact_when_cheaper \\
        --keep_cluster \\
        --out results/pottr_cmp/pottr_tracerx_significance.csv
"""

import argparse
import csv
import re
import subprocess
import sys
from pathlib import Path

from utils import (
    SCRIPT_DIR,
    ensure_dir,
    read_result_pairs,
    parse_items_from_pattern_line,
    load_weights_and_owner,
)

ROOT_LABEL = "0"          # POTTR germline / root node (dropped, like MASTRO's GL)


def endpoints(item):
    """Return (a, sep, b) for an edge item, or None."""
    for sep in ("->-", "-/-", "-?-"):
        if sep in item:
            a, b = item.split(sep, 1)
            return a, sep, b
    return None


def resolved_relations(items, keep_cluster=False):
    """Keep only the ordered/incomparable relations among non-root nodes;
    drop '-?-' (same node) and anything touching the root '0'.
    keep_cluster keeps '-?-' as a same-node relation, endpoints sorted as in
    graphs_all.txt (POTTR writes the pair in either order, sometimes both)."""
    keep = []
    for it in items:
        ep = endpoints(it)
        if ep is None:
            continue
        a, sep, b = ep
        if a == ROOT_LABEL or b == ROOT_LABEL:
            continue
        if sep == "-?-":
            if not keep_cluster:
                continue
            a, b = sorted((a, b))
            it = f"{a}-?-{b}"
            if it in keep:
                continue
        keep.append(it)
    return keep


def load_transaction_sets(graphs_all):
    """Each transaction (tree) -> frozenset of its edge-items."""
    sets = []
    for line in open(graphs_all):
        sets.append(frozenset(line.split()))
    return sets


def match_occurrences(pattern_items, tx_sets):
    """Indices of transactions that contain all of the pattern's items."""
    pat = frozenset(pattern_items)
    return [t for t, s in enumerate(tx_sets) if pat <= s]


def forward_closure(items):
    """Canonical key for a trajectory: the transitive closure of its non-root
    ``a->-b`` relations.

    ``converted_graphs.txt`` stores the transitive *reduction* of the
    trajectory (``run_POTTR.py`` calls ``nx.transitive_reduction``), while
    ``significance_output.txt`` stores the *closure* (its ``load_graph`` calls
    ``nx.transitive_closure_dag``). Comparing the two verbatim only ever
    matches single-edge trajectories, so both sides are closed here.
    """
    succ = {}
    for it in items:
        ep = endpoints(it)
        if ep is None or ep[1] != "->-":
            continue
        a, _, b = ep
        if a == ROOT_LABEL or b == ROOT_LABEL:
            continue
        succ.setdefault(a, set()).add(b)
        succ.setdefault(b, set())

    closed = set()
    for a in succ:
        stack, seen = list(succ[a]), set()
        while stack:
            b = stack.pop()
            if b in seen or b == a:
                continue
            seen.add(b)
            closed.add(f"{a}->-{b}")
            stack.extend(succ.get(b, ()))
    return frozenset(closed)


def pottr_sig_key(items):
    """Key matching a converted_graphs.txt trajectory to its row in
    significance_output.txt: the non-root alterations plus the closure of the
    ``->-`` relations. The closure alone is not enough: POTTR writes '-?-'
    pairs and unordered alterations as plain children of the root, so every
    trajectory without an ordered pair (e.g. KRAS-?-STK11_N) has an empty
    closure and would collide with the others. POTTR keeps one trajectory per
    alteration set at each k, so within a k the key is unique."""
    nodes = set()
    for it in items:
        ep = endpoints(it)
        if ep:
            nodes.update((ep[0], ep[2]))
    nodes.discard(ROOT_LABEL)
    return frozenset(nodes), forward_closure(items)


def load_pottr_pvals(sig_path):
    """Map pottr_sig_key -> POTTR permutation p-value, from a POTTR
    significance_output.txt (';'-separated; edges in the first column)."""
    out = {}
    if not sig_path.exists():
        return out
    with open(sig_path) as f:
        header = f.readline().rstrip("\n").split(";")
        try:
            i_edges = header.index("edges_traj")
            i_pperm = header.index("pval_perm")
        except ValueError:
            return out
        for line in f:
            cols = line.rstrip("\n").split(";")
            if len(cols) <= max(i_edges, i_pperm):
                continue
            edges = re.findall(r"[^\[\],\s]+->-[^\[\],\s]+", cols[i_edges])
            out[pottr_sig_key(edges)] = cols[i_pperm]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pottr_dir", required=True,
                    help="Directory with per-k subdirs (kN/converted_graphs.txt)")
    ap.add_argument("--k_range", default="2,50", help="min,max k (inclusive)")
    ap.add_argument("--graphs_all", required=True)
    ap.add_argument("-w", "--weights", required=True)
    ap.add_argument("--owner", required=True)
    ap.add_argument("--theta", type=float, default=1.0)
    ap.add_argument("--null", choices=["perm", "indep"], default="perm")
    ap.add_argument("--mc_cutoff", type=int, default=8)
    ap.add_argument("--mc_samples", type=int, default=2000)
    ap.add_argument("--exact_when_cheaper", action="store_true",
                    help="exact above the cutoff when placements <= mc_samples")
    ap.add_argument("--n_jobs", type=int, default=1)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--pottr_sig_name", default="significance_output.txt",
                    help="POTTR significance file inside each kN dir (optional)")
    ap.add_argument("--pottr_sig_global", default="",
                    help="extra POTTR significance file covering several k at "
                         "once (written by pottr_force_significance.py); used "
                         "for trajectories no per-k file scored")
    ap.add_argument("--k_dir", action="append", default=[], metavar="K=DIR",
                    help="take k from DIR instead of <pottr_dir>/k<K> (e.g. a k "
                         "run at another thread count); repeatable")
    ap.add_argument("--keep_cluster", action="store_true",
                    help="keep '-?-' pairs as same-node relations (TRACERx encoding)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    kmin, kmax = (int(x) for x in args.k_range.split(","))
    pottr_dir = Path(args.pottr_dir)
    k_dirs = {int(k): Path(d) for k, d in (x.split("=", 1) for x in args.k_dir)}
    global_pvals = load_pottr_pvals(Path(args.pottr_sig_global)) if args.pottr_sig_global else {}
    tx_sets = load_transaction_sets(args.graphs_all)
    _, owner0, _, _ = load_weights_and_owner(args.weights, args.owner)

    # --- gather POTTR trajectories across k, dedup by resolved-relation set ---
    # meta[pattern_key] = dict(n_edges_full, n_rel, n_nodes_rel, occ, n_patients, ks, fwd, pottr_pval)
    meta = {}
    for k in range(kmin, kmax + 1):
        kdir = k_dirs.get(k, pottr_dir / f"k{k}")
        cg = kdir / "converted_graphs.txt"
        if not cg.exists():
            continue
        pottr_pvals = load_pottr_pvals(kdir / args.pottr_sig_name)
        for pl, _ol in read_result_pairs(cg):
            items = parse_items_from_pattern_line(pl)
            if not items:
                continue
            rel = resolved_relations(items, args.keep_cluster)
            if not rel:
                continue
            key = " ".join(sorted(rel))
            if key not in meta:
                occ = match_occurrences(rel, tx_sets)
                fwd = pottr_sig_key(items)
                meta[key] = dict(
                    n_edges_full=len(items),
                    n_rel=len(rel),
                    n_nodes_rel=len({n for it in rel for n in (endpoints(it)[0], endpoints(it)[2])}),
                    occ=occ,
                    n_patients=len({owner0[t] for t in occ if 0 <= t < len(owner0)}),
                    ks=set(),
                    fwd=fwd,
                    pottr_pval="",
                )
            m = meta[key]
            if not m["pottr_pval"]:
                m["pottr_pval"] = pottr_pvals.get(m["fwd"], "")
            m["ks"].add(k)

    if global_pvals:
        for m in meta.values():
            if not m["pottr_pval"]:
                m["pottr_pval"] = global_pvals.get(m["fwd"], "")

    n_pv = sum(1 for m in meta.values() if m["pottr_pval"])
    print(f"[POTTR] {len(meta)} distinct POTTR trajectories over k in "
          f"[{kmin},{kmax}]; {n_pv} with a POTTR p-value", flush=True)

    # --- write a MASTRO-format filtered file and run the MAESTRI test -------
    out_path = Path(args.out)
    ensure_dir(out_path.parent)
    filt = out_path.with_suffix(".pottr_patterns.txt")
    with filt.open("w") as f:
        for key, m in meta.items():
            f.write(f"{key} ({len(m['occ'])})\n")
            f.write(" ".join(str(t) for t in m["occ"]) + "\n")

    sig_csv = out_path.with_suffix(".ensemble_sig.csv")
    subprocess.run([
        sys.executable, str(SCRIPT_DIR / "compute_significance_multitree.py"),
        "-i", str(filt), "-o", str(sig_csv),
        "-w", str(args.weights), "--owner", str(args.owner),
        "--graphs_all", str(args.graphs_all),
        "--test", "both", "--theta", str(args.theta), "--null", args.null,
        "--mc_cutoff", str(args.mc_cutoff), "--mc_samples", str(args.mc_samples),
        "--n_jobs", str(args.n_jobs), "--seed", str(args.seed),
        "--allow-incomplete",
    ] + (["--exact_when_cheaper"] if args.exact_when_cheaper else []), check=True)

    sig = {row["pattern"]: row for row in csv.DictReader(open(sig_csv))}

    # --- merge and write the final comparison table ---
    with out_path.open("w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["k_values", "n_nodes", "n_resolved_rel", "n_edges_full",
                     "support_trees", "support_patients", "s_exp", "s_theta",
                     "mastro_pval_exp", "mastro_pval_theta", "pottr_pval_perm"])
        for key, m in sorted(meta.items(), key=lambda kv: -kv[1]["n_nodes_rel"]):
            s = sig.get(key, {})
            wr.writerow([
                ",".join(str(x) for x in sorted(m["ks"])),
                m["n_nodes_rel"], m["n_rel"], m["n_edges_full"],
                len(m["occ"]), m["n_patients"], s.get("s_exp", ""), s.get("s_theta", ""),
                s.get("pval_exp", ""), s.get("pval_theta", ""), m["pottr_pval"],
            ])
    print(f"[POTTR] written -> {out_path}", flush=True)


if __name__ == "__main__":
    main()
