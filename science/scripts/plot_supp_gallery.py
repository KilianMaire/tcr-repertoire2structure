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
    fig = plt.figure(figsize=(12, 10.5))
    # Explicit layout: the tall complexes are pushed to the outer edges and the
    # enlarged peptide zooms fill the wide empty band down the middle. Positions are
    # figure fractions (x, y, w, h). For each panel: complex axes, inset axes toward
    # the centre, panel-side anchor for the leader.
    ROW = {"top": 0.545, "bot": 0.075}
    CH, CW = 0.43, 0.23            # complex axes height, width
    IW, IH = 0.235, 0.245          # inset axes size
    LAYOUT = [
        # complex_xy,          inset_xy,          left_col
        ((0.005, ROW["top"]), (0.250, ROW["top"] + 0.09), True),   # (a)
        ((0.765, ROW["top"]), (0.515, ROW["top"] + 0.09), False),  # (b)
        ((0.005, ROW["bot"]), (0.250, ROW["bot"] + 0.09), True),   # (c)
        ((0.765, ROW["bot"]), (0.515, ROW["bot"] + 0.09), False),  # (d)
    ]
    for (lab, fname, caption), ((cx, cy), (ix, iy), left_col) in zip(PANELS, LAYOUT):
        ax = fig.add_axes([cx, cy, CW, CH])
        full = _autocrop(plt.imread(FIGS / fname))
        ax.imshow(full)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.text(0.0, 1.0, f"({lab})", transform=ax.transAxes, fontsize=13,
                fontweight="bold", va="top", ha="left", color=PALETTE["ink"])
        ax.text(0.5, -0.02, caption, transform=ax.transAxes, fontsize=9,
                ha="center", va="top", color=PALETTE["ink"])

        iax = fig.add_axes([ix, iy, IW, IH])
        iax.imshow(_autocrop(plt.imread(FIGS / fname.replace(".png", "_inset.png"))))
        iax.set_aspect("equal")
        iax.set_xticks([]); iax.set_yticks([])
        for sp in iax.spines.values():
            sp.set_visible(True); sp.set_edgecolor(PALETTE["mute"]); sp.set_linewidth(1.0)
        iax.text(0.5, 1.02, "peptide", transform=iax.transAxes, fontsize=8.5,
                 ha="center", va="bottom", color=PALETTE["mute"])

        pep = _peptide_xy(full)
        if pep is not None:
            anchor = (0.0, 0.5) if left_col else (1.0, 0.5)   # inset edge facing the complex
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
