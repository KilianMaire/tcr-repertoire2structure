"""Graphical abstract: structural confidence reads presentation, not recognition.

A single self-contained figure (no caption), so unlike the paper figures it does
carry the headline finding and the key numbers as text. Reuses committed PyMOL
renders and the canonical numbers. Left: repertoire to a predicted TCR-pMHC.
Right: two specificity questions, presentation answered (groove confidence
separates a binder from its scramble, AUROC up to 0.99) and recognition not
(pre-registered held-out retrieval 0.61, p 0.24). Bottom: the one-line takeaway.

Usage: python scripts/plot_graphical_abstract.py
"""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from figstyle import PALETTE, CHAIN, SURFACE, FIGS, apply, save

YES, NO = PALETTE["green"], PALETTE["orange"]


def _img(path, pad=10):
    im = plt.imread(FIGS / path)
    rgb = im[..., :3] if im.shape[-1] >= 3 else im
    mask = (rgb < 0.97).any(-1)
    ys, xs = np.where(mask)
    return im[max(ys.min() - pad, 0):ys.max() + pad, max(xs.min() - pad, 0):xs.max() + pad]


def _place(fig, path, rect):
    ax = fig.add_axes(rect)
    ax.imshow(_img(path))
    ax.axis("off")
    return ax


def main():
    apply()
    fig = plt.figure(figsize=(13, 6.6), facecolor=SURFACE)
    bg = fig.add_axes([0, 0, 1, 1]); bg.axis("off")
    bg.set_xlim(0, 1); bg.set_ylim(0, 1)

    # ---- input: repertoire -> predicted complex ------------------------------
    _place(fig, "_struct_view1.png", [0.01, 0.26, 0.26, 0.60])
    bg.text(0.14, 0.90, "10x TCR repertoire", ha="center", fontsize=12,
            fontweight="bold", color=PALETTE["ink"])
    bg.text(0.14, 0.215, "predicted TCR-pMHC\n(Protenix, 5 seeds)", ha="center",
            va="top", fontsize=10, color=PALETTE["mute"])
    for c, lab, x in ((CHAIN["A"], "TCR", 0.075), (CHAIN["C"], "MHC", 0.15),
                      (CHAIN["E"], "peptide", 0.235)):
        bg.text(x, 0.15, lab, ha="center", fontsize=8.5, color=c, fontweight="bold")

    bg.add_patch(FancyArrowPatch((0.275, 0.55), (0.33, 0.55), arrowstyle="-|>",
                                 mutation_scale=22, lw=2, color=PALETTE["ink"]))

    # ---- vertical divider between the two lanes ------------------------------
    bg.plot([0.34, 1.0], [0.53, 0.53], color=PALETTE["gridline"], lw=1.2)

    # ---- lane A: presentation (answered yes) ---------------------------------
    bg.text(0.36, 0.94, "MHC-peptide presentation", ha="left", fontsize=12.5,
            fontweight="bold", color=PALETTE["ink"])
    bg.text(0.36, 0.895, "does the peptide bind the HLA?", ha="left", fontsize=9.5,
            color=PALETTE["mute"])
    _place(fig, "_groove_conf_cognate.png", [0.36, 0.58, 0.20, 0.27])
    _place(fig, "_groove_conf_scramble.png", [0.55, 0.58, 0.20, 0.27])
    bg.text(0.46, 0.565, "cognate (high pLDDT)", ha="center", fontsize=8.8,
            color=PALETTE["blue"], fontweight="bold")
    bg.text(0.65, 0.565, "scramble (low pLDDT)", ha="center", fontsize=8.8,
            color=NO, fontweight="bold")
    bg.text(0.905, 0.79, "✓", ha="center", va="center", fontsize=34,
            color=YES, fontweight="bold")
    bg.text(0.905, 0.685, "groove confidence\nseparates binder\nfrom scramble\nAUROC up to 0.99",
            ha="center", va="top", fontsize=10, color=PALETTE["ink"])

    # ---- lane B: recognition (not confirmed) ---------------------------------
    bg.text(0.36, 0.47, "TCR-peptide recognition", ha="left", fontsize=12.5,
            fontweight="bold", color=PALETTE["ink"])
    bg.text(0.36, 0.425, "does this TCR read the peptide?", ha="left", fontsize=9.5,
            color=PALETTE["mute"])
    _place(fig, "_interface_flu.png", [0.37, 0.13, 0.26, 0.27])
    bg.text(0.905, 0.34, "✗", ha="center", va="center", fontsize=34,
            color=NO, fontweight="bold")
    bg.text(0.905, 0.235, "pre-registered held-out\nretrieval 0.61 (11/18),\np = 0.24: not confirmed",
            ha="center", va="top", fontsize=10, color=PALETTE["ink"])

    # ---- takeaway banner -----------------------------------------------------
    bg.add_patch(FancyBboxPatch((0.06, 0.015), 0.88, 0.075,
                                boxstyle="round,pad=0.006,rounding_size=0.02",
                                linewidth=0, facecolor="#EEF3F5", zorder=0))
    bg.text(0.5, 0.052, "Structural confidence reads presentation, not recognition",
            ha="center", va="center", fontsize=15, fontweight="bold", color=PALETTE["ink"])

    save(fig, "graphical_abstract")


if __name__ == "__main__":
    main()
