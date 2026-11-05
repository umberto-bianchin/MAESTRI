"""
POTTR vs MAESTRI comparison report

1. POTTR's own p-value against MAESTRI's on every trajectory POTTR outputs over
   the k-sweep (counts, quadrants; outputs.csv, plotted by plot_pottr_pvalues.py).
2. POTTR-only trajectories (POTTR p < 0.05, as in the POTTR paper, but not
   significant for MAESTRI after Westfall-Young): their frequencies and why they
   differ. Categories:
     (a)   cluster resolution: POTTR counts more patients than carry the
           exact order in at least one tree (a same-node pair counted as
           ordered);
     (b)   trees as patients: POTTR's support counts more trees than patients;
     (a+b) both;
     (c)   neither: same frequency, the difference is the null.
   (a) is inferred from the counts and then checked tree by tree: every tree
   POTTR lists is 'exact', 'clone' (ordered only by resolving a same-node
   pair) or 'other' (columns trees_exact/_clone/_other, summary.md).
3. MAESTRI-only trajectories (WY-significant, never output by POTTR): state at
   every swept k, read from the per-k outputs and, when present, the solution
   pool (pool_log.json):
     found       POTTR outputs it at this k;
     other       POTTR outputs the same alteration set with another relation
                 (dedup by node set keeps one order);
     infeasible  fewer than k patients have a tree compatible with it (the
                 exact pattern, or its ordered pairs in one node, which POTTR
                 may resolve). This bounds the support POTTR could reach;
                 the conflict graph can only lower it, so "absent" is optimistic;
     absent      the trajectory has the maximum size at this k and at least k
                 patients have a compatible tree, but it is not output
                 (with a full pool: not in the truncated enumeration);
     notmax      the maximum size at this k is a different number of alterations.
   Per trajectory, the status column summarizes these as "found", "other
   relation output", "maximum size but not output" or "never maximum size".

Inputs per dataset: <prefix>.csv and <prefix>.ensemble_sig.csv from
pottr_significance.py (rows aligned by its stable sort on node count), the
per-k POTTR dirs when available, the matched cohort (graphs + owner) and
MAESTRI's significance table + WY thresholds on the same cohort.

    python3 pottr_compare_report.py --datasets tracerx,breast_pool --out results/pottr_cmp/report
    cd results/pottr_cmp/report/tracerx && lualatex maestri_significant_trees.tex
"""

import argparse
from collections import Counter
import csv
import json
from pathlib import Path

import numpy as np

from pottr_significance import endpoints, resolved_relations, ROOT_LABEL
from utils import read_result_pairs, parse_items_from_pattern_line

HERE = Path(__file__).resolve().parent
ALPHA_POTTR = 0.05
FLOOR = 1e-300

# colorblind-safe categorical palette; gray for reference lines
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
GRAY = "#9a9a96"

CATS = [("a+b", "cluster resolution + trees as patients", BLUE),
        ("a", "cluster resolution", ORANGE),
        ("b", "trees as patients", AQUA),
        ("c", "same frequency, different null", YELLOW)]


def dataset_config(name, bc_inputs):
    R = "results/pottr_cmp"
    if name == "tracerx":
        return dict(label="TRACERx", prefix=f"{R}/pottr_tracerx_significance",
                    pottr_dir=f"{R}/pottr_tracerx", keep_cluster=True,
                    graphs="../data/TRACERx/graphs_tracerx.txt",
                    owner="../data/TRACERx/owner_tracerx.txt",
                    maestri="results/exact_tracerx/tracerx_sigma2/significance/"
                            "alg1_expected_uniform_pvalues_exp.csv",
                    wy="results/exact_tracerx/tracerx_wy_sigma2_theta1.0/wy_thresholds.txt")
    common = dict(keep_cluster=False,
                  graphs=f"{bc_inputs}/graphs_all.txt", owner=f"{bc_inputs}/owner.txt",
                  maestri=f"{R}/wy_sigma2/significance/alg1_expected_uniform_pvalues_exp.csv",
                  wy=f"{R}/wy_sigma2_out/wy_thresholds.txt")
    if name == "breast_nopool":
        return dict(label="Breast cancer (no pool)", prefix=f"{R}/pottr_bc_significance",
                    pottr_dir=f"{R}/pottr_bc", **common)
    if name == "breast_pool":
        return dict(label="Breast cancer (pool)", prefix=f"{R}/pottr_bc_pool_t6_significance",
                    pottr_dir=f"{R}/pottr_bc_pool_t6", **common)
    raise SystemExit(f"unknown dataset {name}")


# ---------------------------------------------------------------- helpers
def canon(item):
    """Same-node pairs are written in either order; sort their endpoints."""
    ep = endpoints(item)
    if ep and ep[1] == "-?-":
        a, b = sorted((ep[0], ep[2]))
        return f"{a}-?-{b}"
    return item


def nodes_of(items):
    out = set()
    for it in items:
        ep = endpoints(it)
        if ep:
            out.update((ep[0], ep[2]))
    out.discard(ROOT_LABEL)
    return frozenset(out)


def patient_of(tree_name):
    return tree_name.rsplit("-", 1)[0]


def wy_threshold(path, alpha="0.05"):
    for line in open(path):
        if line.startswith(f"alpha={alpha}"):
            return float(line.split("threshold_exp=")[1].split()[0])
    raise SystemExit(f"no alpha={alpha} in {path}")


def fnum(x):
    return float(x) if x not in ("", None) else float("nan")


def load_cohort(graphs, owner):
    trees = [frozenset(canon(it) for it in line.split()) for line in open(graphs)]
    own = [line.strip() for line in open(owner)]
    assert len(trees) == len(own), (graphs, owner)
    alts = [nodes_of(t) for t in trees]
    return trees, own, alts


def patients_with(pattern_items, cohort):
    """(patients with the exact pattern in >= 1 tree, patients with >= 1 tree
    compatible with it). Compatible: every ordered pair of the pattern is
    ordered the same way or placed in one node, which POTTR may resolve to
    that order; every other relation is exact. The compatible count bounds the
    support POTTR can reach for the pattern."""
    trees, own, _ = cohort
    pat = [canon(it) for it in pattern_items]
    alts = []
    for it in pat:
        a, sep, b = endpoints(it)
        alts.append((it, canon(f"{a}-?-{b}")) if sep == "->-" else (it,))
    exact = {own[t] for t, s in enumerate(trees) if all(it in s for it in pat)}
    compatible = {own[t] for t, s in enumerate(trees) if all(any(x in s for x in alt) for alt in alts)}
    return len(exact), len(compatible)


def tree_index(cohort):
    """POTTR's tree name -> row of the cohort. A dag is named <patient>-<j>,
    j counting the patient's trees in cohort order; the breast converter
    prefixes the numeric owner with 'P'."""
    _, own, _ = cohort
    seen, idx = Counter(), {}
    for t, o in enumerate(own):
        idx[f"{o}-{seen[o]}"] = idx[f"P{o}-{seen[o]}"] = t
        seen[o] += 1
    return idx


def tree_state(rel, tree):
    """How one tree of POTTR's occurrence list holds the trajectory: 'exact'
    (every relation is in the tree), 'clone' (only once a pair the tree puts
    in one node is counted as ordered, as POTTR resolves it) or 'other'
    (some relation is in the tree in neither form)."""
    clone = False
    for it in rel:
        if it in tree:
            continue
        a, sep, b = endpoints(it)
        if sep == "->-" and canon(f"{a}-?-{b}") in tree:
            clone = True
            continue
        return "other"
    return "clone" if clone else "exact"


# ---------------------------------------------------------------- loading
def load_pottr_table(cfg):
    rows = list(csv.DictReader(open(HERE / f"{cfg['prefix']}.csv")))
    sig = list(csv.DictReader(open(HERE / f"{cfg['prefix']}.ensemble_sig.csv")))
    # pottr_significance.py writes the .csv sorted by -n_nodes (stable) over
    # the insertion order of the .ensemble_sig.csv
    sig = sorted(sig, key=lambda r: -int(r["n_nodes"]))
    assert len(rows) == len(sig)
    out = []
    for r, s in zip(rows, sig):
        assert int(r["n_nodes"]) == int(s["n_nodes"]) and \
            abs(fnum(r["s_exp"]) - fnum(s["s_exp"])) < 1e-6, "row alignment broken"
        items = s["pattern"].split()
        out.append(dict(
            key=" ".join(sorted(canon(it) for it in items)), items=items,
            nodes=nodes_of(items), n_alt=int(r["n_nodes"]),
            ks=sorted(int(k) for k in r["k_values"].split(",")),
            exact_trees=int(r["support_trees"]), exact_patients=int(r["support_patients"]),
            s_exp=fnum(r["s_exp"]), s_theta=fnum(r["s_theta"]),
            p_maestri=fnum(r["mastro_pval_exp"]), p_pottr=fnum(r["pottr_pval_perm"]),
        ))
    return out


def load_pottr_runs(cfg, table, cohort):
    """Per k: POTTR's output keys (with its own support) and, when logged, the
    pool's max size, fullness and node sets. Also, per key, the state of every
    tree POTTR lists for it at any k (tree_state)."""
    pdir = HERE / cfg["pottr_dir"]
    runs = {}
    support = {}
    states = {}
    idx = tree_index(cohort)
    cohort_trees = cohort[0]
    kdirs = {int(d.name[1:]): d for d in pdir.glob("k*") if d.name[1:].isdigit()}
    kdirs.update({k: HERE / d for k, d in cfg.get("k_dirs", {}).items()})
    for k, kdir in sorted(kdirs.items()):
        run = {}
        cg = kdir / "converted_graphs.txt"
        if cg.exists():
            keys = set()
            for pl, ol in read_result_pairs(cg):
                items = parse_items_from_pattern_line(pl)
                rel = resolved_relations(items, cfg["keep_cluster"])
                if not rel:
                    continue
                key = " ".join(sorted(canon(it) for it in rel))
                keys.add(key)
                occ = ol.split()
                rel_c = [canon(it) for it in rel]
                st = states.setdefault(key, {})
                for t in occ:
                    st[t] = tree_state(rel_c, cohort_trees[idx[t]])
                sup = support.setdefault(key, (0, 0))
                support[key] = max(sup, (len(occ), len({patient_of(t) for t in occ})))
            run["keys"] = keys
        pl_path = kdir / "pool_log.json"
        if pl_path.exists():
            log = json.load(open(pl_path))
            run["max_alts"] = log["max_size"] - 1        # max_size counts the root
            run["pool_full"] = log["pool_full"]
            run["pool_nodes"] = {frozenset(t["nodes"]) - {ROOT_LABEL}
                                 for t in log["trajectories"]}
        if run:
            runs[k] = run
    # keys and max size from the rescored table where per-k files are missing
    for row in table:
        for k in row["ks"]:
            run = runs.setdefault(k, {})
            run.setdefault("table_keys", set()).add(row["key"])
            run["table_max"] = max(run.get("table_max", 0), row["n_alt"])
    for k, run in runs.items():
        run.setdefault("keys", run.get("table_keys", set()))
        run.setdefault("max_alts", run.get("table_max"))
    return runs, support, states


# ---------------------------------------------------------------- analyses
def analysis_pvalues(table, thr):
    both = [r for r in table if np.isfinite(r["p_maestri"]) and np.isfinite(r["p_pottr"])]
    pm = np.array([max(r["p_maestri"], FLOOR) for r in both])
    pp = np.array([max(r["p_pottr"], FLOOR) for r in both])
    gap = np.log10(pm) - np.log10(pp)
    sig_p, sig_m = pp < ALPHA_POTTR, pm <= thr
    return dict(
        n=len(both), n_pottr_smaller=int((pp < pm).sum()),
        median_orders=float(np.median(gap)), max_orders=float(gap.max()),
        n_ge1_order=int((gap >= 1).sum()),
        pottr_sig=int(sig_p.sum()), maestri_sig=int(sig_m.sum()),
        maestri_unc=int((pm < 0.05).sum()),
        both=int((sig_p & sig_m).sum()), pottr_only=int((sig_p & ~sig_m).sum()),
        maestri_only=int((~sig_p & sig_m).sum()), neither=int((~sig_p & ~sig_m).sum()),
    )


def frequency_category(table, support, states):
    """Every output with POTTR's tree/patient support and its category:
    (a) POTTR counts more patients than carry it exactly, (b) more trees than
    patients, (c) neither. The trees_* columns check (a) tree by tree: how the
    trees POTTR lists (over all k) hold the trajectory, see tree_state."""
    out = []
    for r in table:
        trees, pats = support.get(r["key"], (r["exact_trees"], r["exact_patients"]))
        a = pats > r["exact_patients"]
        b = trees > pats
        cat = "a+b" if a and b else "a" if a else "b" if b else "c"
        st = Counter(states.get(r["key"], {}).values())
        out.append(dict(r, pottr_trees=trees, pottr_patients=pats, cat=cat,
                        trees_exact=st["exact"], trees_clone=st["clone"],
                        trees_other=st["other"]))
    return out


def classify(table, support, states, thr):
    return [r for r in frequency_category(table, support, states)
            if r["p_pottr"] < ALPHA_POTTR and not r["p_maestri"] <= thr]


def maestri_status(cfg, runs, cohort, thr, table):
    rows = [r for r in csv.DictReader(open(HERE / cfg["maestri"]))
            if fnum(r["pval_exp"]) <= thr]
    ks = sorted(runs)
    out = []
    for r in rows:
        items = r["pattern"].split()
        key = " ".join(sorted(canon(it) for it in items))
        nodes, n_alt = nodes_of(items), len(nodes_of(items))
        exact, compatible = patients_with(items, cohort)
        states = {}
        for k in ks:
            run = runs[k]
            out_nodes = {t["nodes"]: t["key"] for t in table if k in t["ks"]}
            if key in run["keys"]:
                states[k] = "found"
            elif nodes in out_nodes or nodes in run.get("pool_nodes", set()):
                states[k] = "other"
            elif k > compatible:
                states[k] = "infeasible"
            elif run["max_alts"] == n_alt:
                states[k] = "absent"
            else:
                states[k] = "notmax"
        found = [k for k in ks if states[k] == "found"]
        other = [k for k in ks if states[k] == "other"]
        absent = [k for k in ks if states[k] == "absent"]
        if found:
            status = "found"
        elif other:
            status = "other relation output"
        elif absent:
            full = [k for k in absent if runs[k].get("pool_full")]
            status = ("maximum size but not output" +
                      (f" (pool full at {len(full)}/{len(absent)} such k)"
                       if any("pool_full" in runs[k] for k in absent) else ""))
        else:
            status = "never maximum size"
        other_keys = sorted({t["key"] for t in table if t["nodes"] == nodes and t["key"] != key})
        out.append(dict(pattern=r["pattern"], key=key, n_alt=n_alt, exact_patients=exact,
                        compatible=compatible, s_exp=fnum(r["s_exp"]), p_maestri=fnum(r["pval_exp"]),
                        status=status, found=found, other=other, absent=absent,
                        other_keys=other_keys, states=states))
    return out


# ---------------------------------------------------------------- output
def kranges(ks):
    if not ks:
        return ""
    ks, out, start = sorted(ks), [], None
    for i, k in enumerate(ks):
        if start is None:
            start = k
        if i + 1 == len(ks) or ks[i + 1] != k + 1:
            out.append(f"{start}" if start == k else f"{start}-{k}")
            start = None
    return ",".join(out)


def pretty(pattern):
    return (pattern.replace("->-", " < ").replace("-?-", " ~ ")
            .replace("-/-", " | "))


def tex_escape(s):
    return s.replace("_", r"\_").replace("~", r"$\sim$").replace("<", r"$\prec$").replace("|", r"$\mid$")


def sci(x):
    if not np.isfinite(x):
        return "--"
    m, e = f"{x:.1e}".split("e")
    return f"${m} \\times 10^{{{int(e)}}}$"


def write_csv(path, rows, cols):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in rows:
            w.writerow([r[c] if not isinstance(r[c], (list, set)) else kranges(r[c]) for c in cols])


# ---------------------------------------------------------------- trees (TikZ)
# Same style as plot_trajectory_trees.py (MAESTRI's figure): germline circle,
# one box per node, ancestor -> descendant left to right. Alterations POTTR or
# MAESTRI put in one node ('-?-') share a box; '-/-' pairs and pairs without
# a relation (POTTR's incomplete posets) get no edge.
TREE_UNITS, TREE_UNIT_W, TREE_GAP_X, TREE_ROW_GAP = 3, 5.2, 0.25, 0.34
TREE_HEAD_H, TREE_STEP_X, TREE_ROOT_X, TREE_FIRST_X = 0.70, 1.6, 0.16, 1.0


def clone_dag(items):
    """(clones, Hasse edges, levels): clones are tuples of alterations joined
    by '-?-', edges the immediate ->- relations between clones."""
    rel = [endpoints(canon(it)) for it in items]
    genes = sorted({g for a, _, b in rel for g in (a, b)} - {ROOT_LABEL})
    root = {g: g for g in genes}

    def find(g):
        while root[g] != g:
            g = root[g]
        return g
    for a, sep, b in rel:
        if sep == "-?-":
            root[find(a)] = find(b)
    clone = {g: tuple(sorted(h for h in genes if find(h) == find(g))) for g in genes}
    clones = sorted(set(clone.values()))
    succ = {c: set() for c in clones}
    for a, sep, b in rel:
        if sep == "->-" and ROOT_LABEL not in (a, b) and clone[a] != clone[b]:
            succ[clone[a]].add(clone[b])
    # transitive closure, then keep the immediate edges only
    closed = {c: set() for c in clones}
    for c in clones:
        stack = list(succ[c])
        while stack:
            d = stack.pop()
            if d not in closed[c]:
                closed[c].add(d)
                stack.extend(succ[d])
    edges = [(c, d) for c in clones for d in closed[c]
             if not any(d in closed[e] for e in closed[c])]
    preds = {c: [u for u, v in edges if v == c] for c in clones}
    level = {}

    def lvl(c):
        if c not in level:
            level[c] = 1 + max((lvl(u) for u in preds[c]), default=0)
        return level[c]
    for c in clones:
        lvl(c)
    return clones, edges, level


def tree_layout(clones, edges, level):
    """{clone: (level, slot)}, slots centered per level, children ordered by
    the mean slot of their parents."""
    depth = max(level.values())
    cols = {l: [c for c in clones if level[c] == l] for l in range(1, depth + 1)}
    height = max(len(v) for v in cols.values())
    pos = {}
    for l in range(1, depth + 1):
        par = lambda c: [pos[u][1] for u, v in edges if v == c and u in pos] or [0]
        col = sorted(cols[l], key=lambda c: (sum(par(c)) / len(par(c)), c))
        off = (height - len(col)) / 2
        for i, c in enumerate(col):
            pos[c] = (l, off + i)
    return pos, height, depth


def tex_gene(g):
    return g.replace("_", r"\_")


def draw_tree_entry(x0, y0, head1, head2, items, crimson):
    clones, edges, level = clone_dag(items)
    pos, height, depth = tree_layout(clones, edges, level)
    big = max(len(c) for c in clones)
    slot_h = 0.58 + 0.30 * (big - 1)
    L = [rf"\node[head] at ({x0:.3f},{y0 - 0.26:.3f}) {{{head1}}};",
         rf"\node[{'tagon' if crimson else 'sub'}] at ({x0 + 0.28:.3f},{y0 - 0.54:.3f}) {{{head2}}};"]
    nst, est = ("hlgene", "hled") if crimson else ("gene", "ed")
    top = y0 - TREE_HEAD_H - slot_h / 2
    xy = lambda c: (x0 + TREE_FIRST_X + (pos[c][0] - 1) * TREE_STEP_X, top - pos[c][1] * slot_h)
    ry = top - (height - 1) / 2 * slot_h
    name = {c: f"n{i}" for i, c in enumerate(clones)}
    L.append(rf"\node[root] (r) at ({x0 + TREE_ROOT_X:.3f},{ry:.3f}) {{}};")
    for c in clones:
        txt = r"\\".join(tex_gene(g) for g in c)
        L.append(r"\node[{}, align=center] ({}) at ({:.3f},{:.3f}) {{{}}};".format(nst, name[c], *xy(c), txt))
    for c in clones:
        if level[c] == 1:
            L.append(rf"\draw[{est}] (r.east) -- ({name[c]}.west);")
    for u, v in edges:
        L.append(rf"\draw[{est}] ({name[u]}.east) -- ({name[v]}.west);")
    width = TREE_FIRST_X + (depth - 1) * TREE_STEP_X + 0.8
    return L, width, TREE_HEAD_H + height * slot_h


def tree_figure(entries, path, rows_per_page=9):
    """entries: (head1, head2, items, crimson). Flow layout as in
    plot_trajectory_trees.figure, one standalone page every rows_per_page
    rows; compile with lualatex."""
    from plot_trajectory_trees import PREAMBLE
    L = [PREAMBLE.replace("plot_trajectory_trees.py", "pottr_compare_report.py")]
    y, col, row_h, rows = 0.0, 0, 0.0, 1
    for head1, head2, items, crimson in entries:
        _, w, h = draw_tree_entry(0, 0, head1, head2, items, crimson)
        span = max(1, -(-w // TREE_UNIT_W))
        if col + span > TREE_UNITS:
            y -= row_h + TREE_ROW_GAP
            col, row_h, rows = 0, 0.0, rows + 1
            if rows > rows_per_page:
                L += [r"\end{tikzpicture}", r"\begin{tikzpicture}"]
                y, rows = 0.0, 1
        lines, _, _ = draw_tree_entry(col * (TREE_UNIT_W + TREE_GAP_X), y, head1, head2, items, crimson)
        L += lines
        col += int(span)
        row_h = max(row_h, h)
    L += [r"\end{tikzpicture}", r"\end{document}"]
    Path(path).write_text("\n".join(L) + "\n")


def tk(ks):
    return "$k$ = " + kranges(ks).replace("-", "--").replace(",", ", ")


def tsci(p):
    """1.1 \\cdot 10^{-27}; plain decimals from 0.1 up."""
    if p >= 0.1:
        return f"{p:.2f}"
    m, e = f"{p:.1e}".split("e")
    return rf"{m} \cdot 10^{{{int(e)}}}"


def maestri_tree_entries(status, runs):
    out = []
    for i, r in enumerate(status, 1):
        head1 = (rf"{{\tiny {i})}}\enspace$p = {tsci(r['p_maestri'])}$\enspace\enspace"
                 rf"$s = {r['s_exp']:g}$\enspace\enspace {r['exact_patients']} patients")
        if r["found"]:
            found, absent = r["found"], r["absent"]
            full = absent and all(runs[k].get("pool_full") for k in absent)
            tail = " (pool full)" if full else ""
            # one line per entry: list the k while they fit, then count them
            # (the full lists are in summary.md and maestri_significant.csv)
            for f_txt, a_txt in ((tk(found), tk(absent)),
                                 (tk(found), f"{len(absent)} other $k$"),
                                 (f"{len(found)} $k$ in {min(found)}--{max(found)}",
                                  f"{len(absent)} other $k$")):
                head2 = f"POTTR: output at {f_txt}" + (f"; missing at {a_txt}{tail}" if absent else "")
                if len(head2.replace("$", "").replace("--", "-")) <= 72:
                    break
        elif r["other"]:
            head2 = f"POTTR: other relation output at {tk(r['other'])}"
        elif r["absent"]:
            head2 = f"POTTR: maximum size but never output ({tk(r['absent'])})"
        else:
            reach = sorted({runs[k]["max_alts"] for k in runs if k <= r["compatible"]})
            head2 = (f"POTTR: never maximum size ({r['compatible']} compatible patients; "
                     + (f"maximum at $k \\le {r['compatible']}$: {', '.join(map(str, reach))} alterations)"
                        if reach else f"$k \\le {r['compatible']}$ not run)"))
        out.append((head1, head2, r["pattern"].split(), not r["found"]))
    return out


def pottr_tree_entries(only):
    out = []
    for i, r in enumerate(sorted(only, key=lambda r: r["p_pottr"]), 1):
        head1 = (rf"{{\tiny {i})}}\enspace$p_{{\mathrm{{POTTR}}}} = {tsci(r['p_pottr'])}$\enspace"
                 rf"$p_{{\mathrm{{MAESTRI}}}} = {tsci(r['p_maestri'])}$")
        head2 = (f"{r['pottr_trees']} trees, {r['pottr_patients']} patients for POTTR; "
                 f"{r['exact_patients']} with the exact order\\enspace({r['cat']})")
        out.append((head1, head2, r["items"], False))
    return out


def write_tex_maestri(status, label, path):
    lines = [r"\begin{tabular}{lrrrll}", r"\toprule",
             r"Trajectory & Alt. & Patients & $p$ (MAESTRI) & POTTR status & $k$ \\", r"\midrule"]
    for r in status:
        ks = r["found"] or r["other"] or r["absent"]
        lines.append(f"{tex_escape(pretty(r['pattern']))} & {r['n_alt']} & {r['exact_patients']} & "
                     f"{sci(r['p_maestri'])} & {r['status']} & {kranges(ks)} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    Path(path).write_text(f"% {label}\n" + "\n".join(lines) + "\n")


def write_tex_only(only, label, path):
    lines = [r"\begin{tabular}{lrrrrrll}", r"\toprule",
             r"Trajectory & Trees & Patients & Exact & $s_{\mathrm{exp}}$ & $p$ (POTTR) & $p$ (MAESTRI) & Cat. \\",
             r"\midrule"]
    for r in sorted(only, key=lambda r: r["p_pottr"]):
        lines.append(f"{tex_escape(pretty(r['key']))} & {r['pottr_trees']} & {r['pottr_patients']} & "
                     f"{r['exact_patients']} & {r['s_exp']:g} & {sci(r['p_pottr'])} & "
                     f"{sci(r['p_maestri'])} & ({r['cat']}) \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    Path(path).write_text(f"% {label}\n" + "\n".join(lines) + "\n")


def run_dataset(name, cfg, out):
    od = out / name
    od.mkdir(parents=True, exist_ok=True)
    thr = wy_threshold(HERE / cfg["wy"])
    table = load_pottr_table(cfg)
    cohort = load_cohort(HERE / cfg["graphs"], HERE / cfg["owner"])
    runs, support, states = load_pottr_runs(cfg, table, cohort)
    missing = {r["key"] for r in table} - set(support)

    pv = analysis_pvalues(table, thr)
    allrows = frequency_category(table, support, states)
    only = classify(table, support, states, thr)
    status = maestri_status(cfg, runs, cohort, thr, table)

    tree_figure(maestri_tree_entries(status, runs), od / "maestri_significant_trees.tex")
    tree_figure(pottr_tree_entries(only), od / "pottr_only_trees.tex")
    write_csv(od / "pottr_only.csv", sorted(only, key=lambda r: r["p_pottr"]),
              ["key", "ks", "n_alt", "pottr_trees", "pottr_patients", "exact_patients",
               "s_exp", "s_theta", "p_pottr", "p_maestri", "cat",
               "trees_exact", "trees_clone", "trees_other"])
    write_csv(od / "outputs.csv", sorted(allrows, key=lambda r: r["p_pottr"]),
              ["key", "ks", "n_alt", "pottr_trees", "pottr_patients", "exact_patients",
               "s_exp", "p_pottr", "p_maestri", "cat",
               "trees_exact", "trees_clone", "trees_other"])
    write_csv(od / "maestri_significant.csv", status,
              ["pattern", "n_alt", "exact_patients", "compatible", "s_exp", "p_maestri",
               "status", "found", "other", "absent", "other_keys"])
    write_tex_only(only, cfg["label"], od / "pottr_only.tex")
    write_tex_maestri(status, cfg["label"], od / "maestri_significant.tex")

    cats = {c: sum(r["cat"] == c for r in only) for c, _, _ in CATS}
    ks = sorted(runs)
    md = [f"# {cfg['label']}", "",
          f"- WY threshold (alpha 0.05, expected support): {thr:.3e}",
          f"- k swept: {kranges(ks)}; unique POTTR outputs: {len(table)}"
          + (f"; POTTR support from converted_graphs for {len(table) - len(missing)}"
             if support else "; no converted_graphs here: POTTR support = exact occurrence"),
          "- max alterations per k: " + ", ".join(
              f"{v}@{kranges([k for k in ks if runs[k]['max_alts'] == v])}"
              for v in sorted({runs[k]['max_alts'] for k in ks}, reverse=True)),
          "", "## 1. p-values", "",
          f"- {pv['n']} outputs with both p-values; POTTR's p is smaller on {pv['n_pottr_smaller']}"
          f" ({100 * pv['n_pottr_smaller'] / pv['n']:.0f}%)",
          f"- log10(p_MAESTRI / p_POTTR): median {pv['median_orders']:.2f}, max {pv['max_orders']:.1f};"
          f" at least one order of magnitude on {pv['n_ge1_order']}",
          f"- POTTR p < 0.05: {pv['pottr_sig']}; MAESTRI p < 0.05 (uncorrected): {pv['maestri_unc']};"
          f" MAESTRI WY-significant among them: {pv['maestri_sig']} (of MAESTRI's"
          f" {len(status)}; POTTR never outputs the remaining {len(status) - pv['maestri_sig']})",
          f"- quadrants (POTTR < 0.05 / MAESTRI WY): both {pv['both']}, POTTR only {pv['pottr_only']},"
          f" MAESTRI only {pv['maestri_only']}, neither {pv['neither']}",
          "", "## 2. POTTR-only", "",
          f"- {len(only)} trajectories; " + ", ".join(f"({c}) {n}" for c, n in cats.items()),
          f"- of these, MAESTRI p < 0.05 uncorrected: {sum(r['p_maestri'] < 0.05 for r in only)}",
          f"- tree check, all outputs (distinct trees POTTR lists over all k): exact"
          f" {sum(r['trees_exact'] for r in allrows)}, ordered only by resolving a pair in one"
          f" clone {sum(r['trees_clone'] for r in allrows)}, other"
          f" {sum(r['trees_other'] for r in allrows)}",
          f"- (a) trajectories with no tree in 'other' (category confirmed): "
          f"{sum(r['cat'] in ('a', 'a+b') and not r['trees_other'] for r in allrows)} of"
          f" {sum(r['cat'] in ('a', 'a+b') for r in allrows)}"
          + "".join(f"; {pretty(r['key'])}: {r['trees_other']} tree(s) in 'other'"
                    for r in allrows if r["trees_other"]),
          "", "## 3. MAESTRI-significant trajectories", ""]
    for r in status:
        md.append(f"- {pretty(r['pattern'])}: {r['n_alt']} alt, {r['exact_patients']} patients exact,"
                  f" {r['compatible']} compatible, p {r['p_maestri']:.2e} -> {r['status']}"
                  + (f"; found k={kranges(r['found'])}" if r["found"] else "")
                  + (f"; other relation at k={kranges(r['other'])}: {', '.join(pretty(x) for x in r['other_keys'])}"
                     if r["other"] else "")
                  + (f"; absent at k={kranges(r['absent'])}" if r["absent"] else ""))
    (od / "summary.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--datasets", default="tracerx")
    ap.add_argument("--bc_inputs", default="results/pottr_cmp/bc_inputs",
                    help="capped breast cohort with graphs_all.txt / owner.txt")
    ap.add_argument("--out", default="results/pottr_cmp/report")
    args = ap.parse_args()
    out = HERE / args.out
    for name in args.datasets.split(","):
        run_dataset(name, dataset_config(name, args.bc_inputs), out)


if __name__ == "__main__":
    main()
