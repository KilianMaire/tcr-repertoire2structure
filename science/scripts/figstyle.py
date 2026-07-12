"""Shared figure style for the paper. One source of truth for palette, type, and
output, so every panel is visually uniform.

Conventions enforced here (they answer the review requirements directly):
  * no figure titles: panels carry a bold "(a)" tag only, never a claim;
  * no result text baked into a figure: interpretation lives in the caption;
  * one Okabe-Ito colorblind-safe palette, assigned by role, never cycled;
  * vector output: every figure is written as PDF (for the manuscript) and PNG
    (for quick view), from the committed CSVs, no GPU.

Import and use:
    from figstyle import PALETTE, SURFACE, panel_label, save, apply
    apply()
    ...
    save(fig, "fig3_retrieval")
"""
from __future__ import annotations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "paper/data"
FIGS = ROOT / "paper/figures"

# Okabe-Ito, validated CVD-safe (scripts checked with the dataviz palette checker).
PALETTE = {
    "blue": "#0072B2",     # TCR alpha / raw / primary series
    "green": "#009E73",    # TCR beta / peptide-identity / "yes" axis
    "orange": "#D55E00",   # de-leaked / cognate / held-out / "the finding"
    "amber": "#E69F00",    # MHC heavy chain
    "pink": "#CC79A7",     # beta-2 microglobulin
    "grey": "#9A9A9A",     # residual / control / sequence baseline
    "gridline": "#DDDDDD",
    "ink": "#333333",
    "mute": "#777777",
}
SURFACE = "#FCFCFB"

# Structure chain colors (match the PyMOL renders and the Fig 1 legend).
CHAIN = {"A": "#0072B2", "B": "#009E73", "C": "#E69F00", "D": "#CC79A7", "E": "#D55E00"}


def apply() -> None:
    """Global rcParams so type and lines are identical across figures."""
    plt.rcParams.update({
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "font.family": "DejaVu Sans",
        "font.size": 9.5,
        "axes.titlesize": 10.5,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "axes.linewidth": 0.9,
        "axes.edgecolor": PALETTE["ink"],
        "xtick.color": PALETTE["ink"],
        "ytick.color": PALETTE["ink"],
        "axes.labelcolor": PALETTE["ink"],
        "text.color": PALETTE["ink"],
        "pdf.fonttype": 42,   # embed real TrueType glyphs, not Type-3 outlines
        "ps.fonttype": 42,
    })


def despine(ax) -> None:
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def panel_label(ax, letter: str, x: float = -0.02, y: float = 1.04) -> None:
    """A bold '(a)' tag in axes coordinates. No claim, no description."""
    ax.text(x, y, f"({letter})", transform=ax.transAxes, fontsize=12,
            fontweight="bold", va="bottom", ha="right", color=PALETTE["ink"])


def save(fig, stem: str) -> None:
    """Write both PDF (manuscript, vector) and PNG (quick view) into paper/figures."""
    FIGS.mkdir(parents=True, exist_ok=True)
    for ext, dpi in ((".pdf", None), (".png", 200)):
        fig.savefig(FIGS / f"{stem}{ext}", dpi=dpi, facecolor=SURFACE, bbox_inches="tight")
    print(f"wrote {FIGS / stem}.pdf + .png")
