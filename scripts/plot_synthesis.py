"""Fig 6: the two-axis map. Where structural confidence is informative.

A synthesis. Two orthogonal specificity questions define a plane: the horizontal
axis is MHC-peptide presentation (does this peptide bind this HLA), the vertical
axis is TCR-peptide recognition (does this TCR read this peptide). The paper's
finding is that a Protenix or AlphaFold confidence readout lives in the lower
right: strong on presentation, blind on recognition. Positions are anchored to the
measured numbers (presentation AUROC up to 0.99; held-out recognition p=0.09, not
confirmed).

Usage: python scripts/plot_synthesis.py     # writes paper/figures/fig6_two_axis_map.png
"""
from __future__ import annotations
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

SURFACE = "#FCFCFB"
C_YES, C_NO, C_SEQ = "#009E73", "#D55E00", "#888888"


def main(out="paper/figures/fig6_two_axis_map.png"):
    fig, ax = plt.subplots(figsize=(8.5, 7.2))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    # quadrant guides
    ax.axhline(0.5, color="#dddddd", lw=1.2, zorder=0)
    ax.axvline(0.5, color="#dddddd", lw=1.2, zorder=0)

    # axes as arrows
    for arrow in [FancyArrowPatch((0.02, 0.02), (0.98, 0.02), arrowstyle="-|>",
                                  mutation_scale=18, color="#333", lw=1.5),
                  FancyArrowPatch((0.02, 0.02), (0.02, 0.98), arrowstyle="-|>",
                                  mutation_scale=18, color="#333", lw=1.5)]:
        ax.add_patch(arrow)
    ax.text(0.5, -0.05, "MHC-peptide presentation  (does the peptide bind the HLA?)",
            ha="center", fontsize=10.5, color="#333")
    ax.text(-0.06, 0.5, "TCR-peptide recognition  (does this TCR read the peptide?)",
            va="center", rotation=90, fontsize=10.5, color="#333")

    # structural confidence: high presentation, low recognition
    ax.scatter([0.82], [0.26], s=620, color=C_NO, edgecolor="white", linewidth=2, zorder=3)
    ax.annotate("structural confidence (ipTM)", (0.82, 0.26), (0.82, 0.40), ha="center",
                fontsize=10, color=C_NO, fontweight="bold",
                arrowprops=dict(arrowstyle="-", color=C_NO, lw=1.2))
    ax.annotate("presentation: AUROC up to 0.99\nrecognition: held-out p=0.09, not confirmed",
                (0.82, 0.26), (0.82, 0.13), ha="center", fontsize=8.6, color="#555",
                arrowprops=dict(arrowstyle="-", color="#ccc"))

    # sequence (TCRdist): places novel TCRs nowhere on this plane (abstains)
    ax.scatter([0.14], [0.14], s=420, color=C_SEQ, edgecolor="white", linewidth=2, zorder=3)
    ax.annotate("sequence (TCRdist)", (0.14, 0.14), (0.30, 0.30), fontsize=9.5, color="#555",
                fontweight="bold", arrowprops=dict(arrowstyle="-", color="#ccc"))
    ax.annotate("abstains on novel TCRs (Honesty Rule 1)", (0.30, 0.255), fontsize=8.3, color="#777")

    # the empty target: a true specificity oracle would sit top-right
    ax.scatter([0.85], [0.85], s=420, facecolor="none", edgecolor="#bbb",
               linewidth=2, linestyle="--", zorder=3)
    ax.annotate("a true specificity oracle\nwould sit here (nothing does yet)",
                (0.85, 0.85), (0.85, 0.70), ha="center", va="top", fontsize=8.5, color="#999",
                arrowprops=dict(arrowstyle="-", color="#ddd"))

    ax.set_xlim(-0.12, 1.05)
    ax.set_ylim(-0.12, 1.05)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title("Structural confidence reads presentation, not recognition",
                 fontsize=13, fontweight="bold", loc="left", y=1.02)
    fig.tight_layout()
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    print(f"wrote {out}")


if __name__ == "__main__":
    main(*(sys.argv[1:2] or []))
