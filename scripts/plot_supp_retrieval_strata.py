"""Fig S1: full retrieval battery, all folds vs reconstructed-only.

Reads paper/data/tcr_retrieval_top1.csv. Two panels (A*02:01 discovery,
A*11:01 held-out); each readout gets an all-folds bar and a reconstructed-only
bar, so the effect of excluding poly-G stub folds is visible per readout. Chance
is per panel (0.25 vs 0.5). No interpretation on the axes.

Usage: python scripts/plot_supp_retrieval_strata.py
"""
from __future__ import annotations
import csv, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import matplotlib.pyplot as plt
from figstyle import PALETTE, SURFACE, DATA, apply, despine, panel_label, save

C_ALL, C_RECON = PALETTE["grey"], PALETTE["blue"]
PANELS = [("A*02:01  (discovery)", "panel1", 0.25, "a"),
          ("A*11:01  (held-out)", "hla_a1101", 0.5, "b")]


def _rows():
    with (DATA / "tcr_retrieval_top1.csv").open() as f:
        return list(csv.DictReader(f))


def main():
    apply()
    rows = _rows()
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.0))
    for ax, (title, run, chance, lab) in zip(axes, PANELS):
        allr = {r["readout"]: float(r["top1"]) for r in rows if r["run"] == run and r["panel"] == "all"}
        rec = {r["readout"]: float(r["top1"]) for r in rows if r["run"] == run and r["panel"] == "reconstructed"}
        order = sorted(allr, key=lambda m: rec.get(m, 0))
        h = 0.38
        for y, m in enumerate(order):
            ax.barh(y - h / 2, allr[m], height=h, color=C_ALL, edgecolor=SURFACE,
                    linewidth=1.1, label="all folds" if y == 0 else None, zorder=2)
            ax.barh(y + h / 2, rec[m], height=h, color=C_RECON, edgecolor=SURFACE,
                    linewidth=1.1, label="reconstructed only" if y == 0 else None, zorder=2)
        ax.axvline(chance, color=PALETTE["mute"], lw=1.3, ls="--", zorder=1)
        ax.annotate(f"chance {chance:g}", (chance, len(order) - 0.4), fontsize=8.5,
                    color=PALETTE["mute"], ha="center", va="bottom")
        ax.set_yticks(range(len(order)))
        ax.set_yticklabels(order, fontsize=8.5)
        ax.set_xlim(0, 0.78)
        ax.set_xlabel(f"Top-1 retrieval  ({title})")
        despine(ax)
        ax.legend(fontsize=8.5, frameon=False, loc="lower right")
        panel_label(ax, lab, x=-0.42)
    fig.tight_layout()
    save(fig, "figS1_retrieval_strata")


if __name__ == "__main__":
    main()
