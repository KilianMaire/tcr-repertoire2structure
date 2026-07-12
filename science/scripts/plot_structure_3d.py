"""Fig 1: "the system" - what a TCR-pMHC complex is and how the pipeline builds one.

Five-panel orientation figure. (a) is a matplotlib pipeline schematic (no
rendered structure needed). (b)-(d) embed the PyMOL ray-traced renders from
scripts/render_structure_pymol.py (never re-run PyMOL here), autocropped and
annotated with short text + leader lines pointing at named features. (e) is a
compact schematic of the two specificity axes the rest of the paper measures
(presentation, recognition). No title, no interpretation: the reading lives
in the caption.

Usage: python scripts/plot_structure_3d.py
"""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Ellipse
from figstyle import PALETTE, CHAIN, SURFACE, FIGS, apply, save

V1 = FIGS / "_struct_view1.png"
V2 = FIGS / "_struct_view2.png"
VI = FIGS / "_interface_flu.png"


def _autocrop(img, pad=12):
    rgb = img[..., :3] if img.shape[-1] == 4 else img
    mask = (rgb < 0.97).any(-1)
    if not mask.any():
        return img
    ys, xs = np.where(mask)
    y0, y1 = max(ys.min() - pad, 0), min(ys.max() + pad, img.shape[0])
    x0, x1 = max(xs.min() - pad, 0), min(xs.max() + pad, img.shape[1])
    return img[y0:y1, x0:x1]


def _show_with_margin(ax, path, left=0.06, right=0.06, top=0.06, bottom=0.06):
    """Autocrop, then reserve blank (SURFACE-colored) margin around the image
    on each side so annotation text has somewhere to sit that is not on top
    of ribbon. Returns (w, h) of the cropped image in data/pixel coords, with
    y increasing downward (matplotlib image convention)."""
    img = _autocrop(plt.imread(path))
    h, w = img.shape[0], img.shape[1]
    ax.imshow(img, extent=[0, w, h, 0], zorder=2)
    ax.set_xlim(-left * w, w * (1 + right))
    ax.set_ylim(h * (1 + bottom), -top * h)
    # keep the axes' physical box == the full gridspec cell (so panel labels
    # and neighbouring panels line up), and let the *data limits* absorb the
    # aspect correction instead of shrinking/recentering the box.
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_facecolor(SURFACE)
    ax.axis("off")
    return w, h


def _leader(ax, text, xy, xytext, color=PALETTE["ink"], ha="center", va="center"):
    ax.annotate(
        text, xy=xy, xycoords="data", xytext=xytext, textcoords="data",
        fontsize=7.8, color=color, ha=ha, va=va, zorder=3,
        arrowprops=dict(arrowstyle="-", lw=0.8, color=color, shrinkA=0, shrinkB=3),
    )


# ---------------------------------------------------------------- (a) pipeline
STAGES = [
    "10x TCR\nrepertoire",
    "annotate\n(TCRdist +\nleakage guard)",
    "build\nTCR-pMHC",
    "fold\n(Protenix,\n5 seeds)",
    "QC\n(scramble\nmargin)",
]


def _panel_a(ax):
    n = len(STAGES)
    xs = np.linspace(0.5, n - 0.5, n)
    box_w, box_h = 0.86, 0.82
    y = 0.5
    for x, label in zip(xs, STAGES):
        box = FancyBboxPatch(
            (x - box_w / 2, y - box_h / 2), box_w, box_h,
            boxstyle="round,pad=0.02,rounding_size=0.08",
            linewidth=1.3, edgecolor=PALETTE["blue"], facecolor=SURFACE, zorder=2,
        )
        ax.add_patch(box)
        ax.text(x, y, label, ha="center", va="center", fontsize=8.2,
                 color=PALETTE["ink"], zorder=3)
    for x0, x1 in zip(xs[:-1], xs[1:]):
        arr = FancyArrowPatch(
            (x0 + box_w / 2, y), (x1 - box_w / 2, y),
            arrowstyle="-|>", mutation_scale=11, lw=1.3,
            color=PALETTE["mute"], zorder=1,
        )
        ax.add_patch(arr)
    ax.set_xlim(0, n)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(-0.012, 1.05, "(a)", transform=ax.transAxes, fontsize=12,
             fontweight="bold", va="bottom", ha="right", color=PALETTE["ink"])


# ------------------------------------------------------------- (b) whole complex
def _panel_b(ax):
    w, h = _show_with_margin(ax, V1, left=0.11, right=0.10, top=0.05, bottom=0.04)
    # feature coordinates read off the cropped image (x right, y down)
    _leader(ax, "TCR Vα", (0.62 * w, 0.22 * h), (-0.09 * w, 0.02 * h), color=CHAIN["A"], ha="right", va="bottom")
    _leader(ax, "TCR Vβ", (0.92 * w, 0.42 * h), (1.08 * w, 0.30 * h), color=CHAIN["B"], ha="left", va="center")
    _leader(ax, "MHC class I", (0.30 * w, 0.68 * h), (-0.09 * w, 0.80 * h), color=CHAIN["C"], ha="right", va="center")
    _leader(ax, "β₂m", (0.24 * w, 0.42 * h), (-0.09 * w, 0.42 * h), color=CHAIN["D"], ha="right", va="center")
    _leader(ax, "peptide", (0.55 * w, 0.55 * h), (0.55 * w, 1.05 * h), color=CHAIN["E"], ha="center", va="top")
    ax.text(-0.012, 1.02, "(b)", transform=ax.transAxes, fontsize=12,
             fontweight="bold", va="bottom", ha="right", color=PALETTE["ink"])


# --------------------------------------------------------- (c) down the groove
def _panel_c(ax):
    w, h = _show_with_margin(ax, V2, left=0.06, right=0.06, top=0.14, bottom=0.14)
    _leader(ax, "α1 helix", (0.30 * w, 0.16 * h), (0.20 * w, -0.10 * h), color=CHAIN["C"], ha="center", va="bottom")
    _leader(ax, "α2 helix", (0.68 * w, 0.86 * h), (0.78 * w, 1.10 * h), color=CHAIN["C"], ha="center", va="top")
    _leader(ax, "β-sheet\nfloor", (0.42 * w, 0.30 * h), (0.06 * w, 0.20 * h), color=PALETTE["ink"], ha="right", va="center")
    _leader(ax, "peptide", (0.55 * w, 0.55 * h), (0.94 * w, 0.55 * h), color=CHAIN["E"], ha="left", va="center")
    ax.text(-0.012, 1.02, "(c)", transform=ax.transAxes, fontsize=12,
             fontweight="bold", va="bottom", ha="right", color=PALETTE["ink"])


# ------------------------------------------------------------ (d) interface
def _panel_d(ax):
    w, h = _show_with_margin(ax, VI, left=0.06, right=0.06, top=0.12, bottom=0.06)
    _leader(ax, "peptide", (0.55 * w, 0.30 * h), (0.55 * w, -0.06 * h), color=CHAIN["E"], ha="center", va="bottom")
    _leader(ax, "CDR loops\n(Vα)", (0.68 * w, 0.75 * h), (0.98 * w, 0.90 * h), color=CHAIN["A"], ha="left", va="center")
    _leader(ax, "CDR loops\n(Vβ)", (0.32 * w, 0.72 * h), (0.02 * w, 0.90 * h), color=CHAIN["B"], ha="right", va="center")
    ax.text(-0.012, 1.02, "(d)", transform=ax.transAxes, fontsize=12,
             fontweight="bold", va="bottom", ha="right", color=PALETTE["ink"])


# ---------------------------------------------------------- (e) two axes
def _axis_row(ax, y, label, left_color, question):
    ax.add_patch(Ellipse((0.20, y), 0.16, 0.16, color=left_color, zorder=2))
    ax.add_patch(Ellipse((0.40, y), 0.10, 0.10, color=PALETTE["orange"], zorder=3))
    ax.annotate("", xy=(0.62, y), xytext=(0.47, y),
                arrowprops=dict(arrowstyle="-|>", lw=1.2, color=PALETTE["mute"]))
    ax.text(0.66, y, "?", fontsize=11, fontweight="bold", ha="left", va="center", color=PALETTE["ink"])
    ax.text(0.02, y + 0.17, label, fontsize=8.4, fontweight="bold", ha="left", va="bottom", color=PALETTE["ink"])
    ax.text(0.02, y - 0.17, question, fontsize=7.6, ha="left", va="top", color=PALETTE["mute"])


def _panel_e(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    _axis_row(ax, 0.74, "presentation", CHAIN["C"], "does the peptide bind HLA?")
    _axis_row(ax, 0.26, "recognition", CHAIN["A"], "does this TCR read the peptide?")
    ax.text(-0.03, 1.05, "(e)", transform=ax.transAxes, fontsize=12,
             fontweight="bold", va="bottom", ha="right", color=PALETTE["ink"])


def main():
    apply()
    fig = plt.figure(figsize=(13, 8), facecolor=SURFACE)
    gs = fig.add_gridspec(2, 4, height_ratios=[0.32, 1.0], hspace=0.22, wspace=0.08,
                           left=0.03, right=0.99, top=0.95, bottom=0.03)

    _panel_a(fig.add_subplot(gs[0, :]))
    _panel_b(fig.add_subplot(gs[1, 0]))
    _panel_c(fig.add_subplot(gs[1, 1]))
    _panel_d(fig.add_subplot(gs[1, 2]))
    _panel_e(fig.add_subplot(gs[1, 3]))

    save(fig, "fig1_structure")


if __name__ == "__main__":
    main()
