"""Fig S5: gallery of predicted TCR-pMHC complexes across epitopes.

Composes four already-rendered PyMOL complexes (same canonical TCR-up orientation)
into a 2x2 gallery with one shared chain-colour legend. Shows the pipeline produces
consistent well-folded complexes across epitopes and both HLAs. No title, no result
text. Renders come from scripts/render_structure_pymol.py (single mode); this only
composes them, no GPU.

Usage: python scripts/plot_supp_gallery.py
"""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, ConnectionPatch
from figstyle import PALETTE, CHAIN, FIGS, apply, save

PANELS = [
    ("a", "_gallery1_GILGFVFTL.png", "GILGFVFTL / HLA-A*02:01"),
    ("b", "_gallery2_ELAGIGILTV.png", "ELAGIGILTV / HLA-A*02:01"),
    ("c", "_gallery3_FLYALALLL.png", "FLYALALLL / HLA-A*02:01"),
    ("d", "_gallery4_IVTDFSVIK.png", "IVTDFSVIK / HLA-A*11:01"),
]
LEGEND = [(CHAIN["A"], "TCR α chain"), (CHAIN["B"], "TCR β chain"),
          (CHAIN["C"], "MHC class I heavy chain"), (CHAIN["D"], "β₂-microglobulin"),
          (CHAIN["E"], "peptide antigen")]


def _autocrop(img, pad=12):
    rgb = img[..., :3] if img.shape[-1] == 4 else img
    mask = (rgb < 0.97).any(-1)
    if not mask.any():
        return img
    ys, xs = np.where(mask)
    y0, y1 = max(ys.min() - pad, 0), min(ys.max() + pad, img.shape[0])
    x0, x1 = max(xs.min() - pad, 0), min(xs.max() + pad, img.shape[1])
    return img[y0:y1, x0:x1]


def _peptide_xy(img):
    """Pixel centroid (col, row) of the peptide sticks, found by the peptide orange
    (high red, low green, near-zero blue), which is separable from the amber MHC by a
    larger red-minus-green gap. Returns None if no peptide pixels are found."""
    rgb = img[..., :3] if img.shape[-1] == 4 else img
    r, g, b = rgb[..., 0], rgb[..., 1], rgb[..., 2]
    mask = (r > 0.5) & (b < 0.25) & ((r - g) > 0.38)
    if not mask.any():
        return None
    ys, xs = np.where(mask)
    return xs.mean(), ys.mean()


def main():
    apply()
    fig = plt.figure(figsize=(11, 10.5))
    gs = fig.add_gridspec(2, 2, hspace=0.10, wspace=0.04,
                          left=0.02, right=0.98, top=0.98, bottom=0.10)
    for i, (lab, fname, caption) in enumerate(PANELS):
        ax = fig.add_subplot(gs[i // 2, i % 2])
        full = _autocrop(plt.imread(FIGS / fname))
        ax.imshow(full)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.text(0.0, 1.0, f"({lab})", transform=ax.transAxes, fontsize=13,
                fontweight="bold", va="top", ha="left", color=PALETTE["ink"])
        ax.text(0.5, -0.02, caption, transform=ax.transAxes, fontsize=9,
                ha="center", va="top", color=PALETTE["ink"])
        # inset: an enlarged zoom on the peptide in the groove, placed in the empty
        # lower corner toward the figure centre, with a leader line from the peptide
        # on the full complex.
        left_col = (i % 2 == 0)
        bounds = [0.57, 0.03, 0.42, 0.32] if left_col else [0.01, 0.03, 0.42, 0.32]
        iax = ax.inset_axes(bounds)
        iax.imshow(_autocrop(plt.imread(FIGS / fname.replace(".png", "_inset.png"))))
        iax.set_aspect("equal")
        iax.set_xticks([]); iax.set_yticks([])
        for sp in iax.spines.values():
            sp.set_visible(True); sp.set_edgecolor(PALETTE["mute"]); sp.set_linewidth(1.0)
        iax.text(0.5, 1.02, "peptide", transform=iax.transAxes, fontsize=8,
                 ha="center", va="bottom", color=PALETTE["mute"])
        pep = _peptide_xy(full)
        if pep is not None:
            anchor = (0.12, 1.0) if left_col else (0.88, 1.0)  # inset corner facing the peptide
            con = ConnectionPatch(xyA=anchor, coordsA=iax.transAxes,
                                  xyB=pep, coordsB=ax.transData,
                                  lw=1.1, ls=(0, (4, 3)), color=PALETTE["ink"], zorder=5)
            ax.add_artist(con)
    fig.legend(handles=[Patch(color=c, label=l) for c, l in LEGEND],
               loc="lower center", ncol=5, frameon=False, fontsize=9,
               bbox_to_anchor=(0.5, 0.005))
    save(fig, "figS5_complex_gallery")


if __name__ == "__main__":
    main()
