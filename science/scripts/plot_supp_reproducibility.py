"""Fig S4: sample-to-sample reproducibility of the five Protenix folds, cognate
vs composition-scramble, for one held-out A*11:01 clonotype.

Reads paper/data/per_sample_readouts.csv. Panels:
  (a) iptm_groove: stable, clean cognate/scramble separation across samples.
  (b) iptm_TCRpep_max: stable, small separation, both low.
  (c) ranking_score: stable, but cognate and scramble overlap.

No result-interpretation sentence lives on the figure; that goes in the caption.

Usage: python scripts/plot_supp_reproducibility.py
"""
from __future__ import annotations
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import matplotlib.pyplot as plt
from figstyle import PALETTE, SURFACE, DATA, apply, despine, panel_label, save

COLOR = {"cognate": PALETTE["blue"], "scramble": PALETTE["orange"]}
LABEL = {"cognate": "cognate", "scramble": "scramble"}

PANELS = [
    ("a", "iptm_groove", "iptm_groove", (0.70, 1.00)),
    ("b", "iptm_TCRpep_max", "iptm_TCRpep_max", (0.30, 0.46)),
    ("c", "ranking_score", "ranking_score", (0.88, 0.96)),
]


def _rows():
    with (DATA / "per_sample_readouts.csv").open() as f:
        return list(csv.DictReader(f))


def main():
    apply()
    rows = _rows()
    by_construct = {
        c: sorted([r for r in rows if r["construct"] == c], key=lambda r: int(r["sample"]))
        for c in ("cognate", "scramble")
    }

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))

    for ax, (letter, col, ylabel, ylim) in zip(axes, PANELS):
        despine(ax)
        for construct in ("cognate", "scramble"):
            recs = by_construct[construct]
            xs = [int(r["sample"]) for r in recs]
            ys = [float(r[col]) for r in recs]
            mean = sum(ys) / len(ys)
            color = COLOR[construct]
            ax.scatter(xs, ys, color=color, s=34, zorder=3,
                       label=LABEL[construct], edgecolor=SURFACE, linewidth=0.6)
            ax.hlines(mean, -0.4, 4.4, color=color, lw=1.2, ls="--", zorder=2, alpha=0.85)
        ax.set_xlim(-0.4, 4.4)
        ax.set_xticks([0, 1, 2, 3, 4])
        ax.set_xlabel("Protenix sample")
        ax.set_ylabel(ylabel)
        ax.set_ylim(*ylim)
        panel_label(ax, letter)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False,
               bbox_to_anchor=(0.5, -0.04), fontsize=9)

    fig.subplots_adjust(left=0.06, right=0.985, top=0.90, bottom=0.24, wspace=0.42)
    save(fig, "figS4_reproducibility")


if __name__ == "__main__":
    main()
