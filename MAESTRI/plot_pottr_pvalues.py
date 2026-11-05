"""Paper figure: POTTR's p-value against MAESTRI's on every trajectory POTTR
outputs, one log-log panel per cohort, colored by why POTTR's frequency
differs (pottr_compare_report.py categories).

Reads results/pottr_cmp/report/<dataset>/outputs.csv, written by
pottr_compare_report.py, and the WY threshold of each cohort.

Usage:
    python3 plot_pottr_pvalues.py --datasets tracerx,breast_pool \
        --out results/pottr_cmp/report/pvalues_paper.pdf
"""

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.lines import Line2D
import numpy as np

from pottr_compare_report import (ALPHA_POTTR, FLOOR, GRAY, BLUE, ORANGE, AQUA, YELLOW,
                                  dataset_config, wy_threshold)

HERE = Path(__file__).resolve().parent
TITLES = {"tracerx": "TRACERx", "breast_pool": "Breast cancer"}
# each category also gets its own marker, since the last two colors have low contrast
CATS = [("a+b", "resolved pair + trees as patients", BLUE, "o"),
        ("a", "resolved pair", ORANGE, "s"),
        ("b", "trees as patients", AQUA, "^"),
        ("c", "same frequency", YELLOW, "D")]


def load(name):
    rows = list(csv.DictReader(open(HERE / f"results/pottr_cmp/report/{name}/outputs.csv")))
    pp = np.array([max(float(r["p_pottr"]), FLOOR) for r in rows])
    pm = np.array([max(float(r["p_maestri"]), FLOOR) for r in rows])
    cat = np.array([r["cat"] for r in rows])
    return pp, pm, cat


def style(ax):
    ax.tick_params(labelsize=6.5)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def panel(ax_tail, ax, name, title, x_min):
    """One cohort on a broken x axis: the main log-log plot, and on its left a
    narrow strip with the trajectories whose POTTR p-value is below x_min,
    which would otherwise squeeze every other trajectory into the right edge."""
    pp, pm, cat = load(name)
    thr = wy_threshold(HERE / dataset_config(name, None)["wy"])
    tail = pp < x_min
    lo_y = 10 ** np.floor(np.log10(min(pm.min(), thr)) - 0.3)
    for a in (ax_tail, ax):
        a.set_yscale("log"); a.set_ylim(lo_y, 2)
        a.axhline(thr, ls=":", lw=0.8, color="0.35", zorder=1)
        a.yaxis.set_major_locator(matplotlib.ticker.LogLocator(numticks=6))
        a.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        style(a)
    ax.set_xscale("log"); ax.set_xlim(x_min, 2)
    d = max(x_min, lo_y)
    ax.plot([d, 2], [d, 2], ls="--", lw=0.7, color=GRAY, zorder=1)
    ax.axvline(ALPHA_POTTR, ls=":", lw=0.8, color="0.35", zorder=1)
    ax.text(ALPHA_POTTR, lo_y * 1.6, r" $0.05$", fontsize=6, color="0.35", ha="left", va="bottom")
    ax.text(x_min * 2, thr * 0.7, "WY threshold", fontsize=6, color="0.35", ha="left", va="top")
    for c, _, col, mk in CATS:
        for a, m, x in ((ax, (cat == c) & ~tail, pp), (ax_tail, (cat == c) & tail, np.zeros_like(pp))):
            a.scatter(x[m], pm[m], s=11, marker=mk, color=col, edgecolor="white",
                      linewidth=0.35, alpha=0.9, zorder=3)
    e = int(round(np.log10(x_min)))
    # no tick at x_min itself: it would sit next to the strip's label
    ax.set_xticks([10.0 ** t for t in range(e + 2, 1, 2)])
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", left=False, labelleft=False)
    ax.set_title(f"{title}, {len(pp)} trajectories", fontsize=8)
    ax.set_xlabel("POTTR $p$-value", fontsize=7.5)
    ax_tail.set_xlim(-1, 1)
    ax_tail.set_xticks([0], [f"$<10^{{{e}}}$"])
    # break marks on the x axis between the strip and the main plot
    kw = dict(color="k", lw=0.8, clip_on=False, transform=ax_tail.transAxes)
    ax_tail.plot([1.0, 1.25], [-0.02, 0.02], **kw)
    ax_tail.plot([1.25, 1.5], [-0.02, 0.02], **kw)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--datasets", default="tracerx,breast_pool")
    ap.add_argument("--out", default="results/pottr_cmp/report/pvalues_paper.pdf")
    ap.add_argument("--x_min", type=float, default=1e-12,
                    help="POTTR p-values below this go to the strip left of the break")
    args = ap.parse_args()
    names = args.datasets.split(",")
    fig = plt.figure(figsize=(2.55 * len(names), 3.3))
    # per cohort: tail strip, main plot, then a spacer before the next cohort
    ratios = [r for _ in names for r in (0.9, 7, 2.2)][:-1]
    gs = fig.add_gridspec(1, len(ratios), width_ratios=ratios, wspace=0.08,
                          left=0.1, right=0.95, top=0.92, bottom=0.3)
    for i, name in enumerate(names):
        ax_tail = fig.add_subplot(gs[0, 3 * i])
        ax = fig.add_subplot(gs[0, 3 * i + 1], sharey=ax_tail)
        panel(ax_tail, ax, name, TITLES.get(name, name), args.x_min)
        if i == 0:
            ax_tail.set_ylabel("MAESTRI $p$-value", fontsize=7.5)
    handles = [Line2D([], [], ls="", marker=mk, color=col, markeredgecolor="white",
                      markersize=4.5, label=lab) for _, lab, col, mk in CATS]
    handles += [Line2D([], [], ls="--", lw=0.7, color=GRAY, label="equal $p$-values"),
                Line2D([], [], ls=":", lw=0.8, color="0.35",
                       label="thresholds: POTTR 0.05, MAESTRI WY")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False,
               fontsize=6.5, handletextpad=0.3, columnspacing=1.2, handlelength=1.6)
    out = HERE / args.out
    fig.savefig(out)
    fig.savefig(out.with_suffix(".png"), dpi=200)
    print(f"written -> {out}")


if __name__ == "__main__":
    main()
