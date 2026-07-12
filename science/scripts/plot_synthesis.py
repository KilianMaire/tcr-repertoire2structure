"""Fig 6: the two-axis map of where structural confidence is informative.

Two orthogonal specificity questions define a plane: the horizontal axis is
MHC-peptide presentation (does this peptide bind this HLA), the vertical axis is
TCR-peptide recognition (does this TCR read this peptide). Structural confidence
sits lower-right: strong on presentation, blind on recognition. The supporting
numbers live in the caption, not on the figure.

Usage: python scripts/plot_synthesis.py
"""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from figstyle import PALETTE, apply, save

C_YES, C_NO, C_SEQ = PALETTE["green"], PALETTE["orange"], PALETTE["grey"]


def main():
    apply()
    fig, ax = plt.subplots(figsize=(8.0, 7.0))

    ax.axhline(0.5, color=PALETTE["gridline"], lw=1.2, zorder=0)
    ax.axvline(0.5, color=PALETTE["gridline"], lw=1.2, zorder=0)
    for arrow in [FancyArrowPatch((0.02, 0.02), (0.98, 0.02), arrowstyle="-|>",
                                  mutation_scale=18, color=PALETTE["ink"], lw=1.5),
                  FancyArrowPatch((0.02, 0.02), (0.02, 0.98), arrowstyle="-|>",
                                  mutation_scale=18, color=PALETTE["ink"], lw=1.5)]:
        ax.add_patch(arrow)
    ax.text(0.5, -0.05, "MHC-peptide presentation  (does the peptide bind the HLA?)",
            ha="center", fontsize=10.5, color=PALETTE["ink"])
    ax.text(-0.06, 0.5, "TCR-peptide recognition  (does this TCR read the peptide?)",
            va="center", rotation=90, fontsize=10.5, color=PALETTE["ink"])

    # structural confidence: high presentation, low recognition
    ax.scatter([0.82], [0.26], s=620, color=C_NO, edgecolor="white", linewidth=2, zorder=3)
    ax.annotate("structural confidence\n(ipTM)", (0.82, 0.26), (0.82, 0.40), ha="center",
                fontsize=10, color=C_NO, fontweight="bold",
                arrowprops=dict(arrowstyle="-", color=C_NO, lw=1.2))

    # sequence (TCRdist): abstains on novel TCRs, off the plane
    ax.scatter([0.14], [0.14], s=420, color=C_SEQ, edgecolor="white", linewidth=2, zorder=3)
    ax.annotate("sequence\n(TCRdist)", (0.14, 0.14), (0.30, 0.30), ha="center", fontsize=9.5,
                color=PALETTE["mute"], fontweight="bold",
                arrowprops=dict(arrowstyle="-", color=PALETTE["gridline"]))

    # the empty target: a specificity oracle would sit top-right
    ax.scatter([0.85], [0.85], s=420, facecolor="none", edgecolor="#bbb",
               linewidth=2, linestyle="--", zorder=3)
    ax.annotate("specificity oracle\n(unoccupied)", (0.85, 0.85), (0.85, 0.70), ha="center",
                va="top", fontsize=9, color=PALETTE["mute"],
                arrowprops=dict(arrowstyle="-", color=PALETTE["gridline"]))

    ax.set_xlim(-0.12, 1.05)
    ax.set_ylim(-0.12, 1.05)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    fig.tight_layout()
    save(fig, "fig6_two_axis_map")


if __name__ == "__main__":
    main()
