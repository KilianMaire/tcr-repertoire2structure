"""Fig S2: per-residue confidence in the groove, cognate vs composition-scramble.

Composes the two PyMOL groove_conf renders (peptide sticks colored by per-residue
pLDDT, residues labelled) for one A*11:01 clonotype: its cognate peptide and a
composition-scramble of the same residues. A shared pLDDT colorbar keys both.
The point (caption): the same amino acids are placed with high confidence when
they are the true epitope and low confidence when shuffled.

Renders come from:
  micromamba run -n pymol python scripts/render_structure_pymol.py groove_conf <cif> <png>

Usage: python scripts/plot_supp_groove_conf.py
"""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from figstyle import PALETTE, FIGS, apply, save

PLDDT_CMAP = LinearSegmentedColormap.from_list("plddt_rwb", ["#D7191C", "#FFFFFF", "#2C7BB6"])
PMIN, PMAX = 50, 95


def _crop(img, pad=6):
    rgb = img[..., :3]
    mask = (rgb < 0.97).any(-1)
    ys, xs = np.where(mask)
    return img[max(ys.min() - pad, 0):ys.max() + pad, max(xs.min() - pad, 0):xs.max() + pad]


def main():
    apply()
    cog = FIGS / "_groove_conf_cognate.png"
    scr = FIGS / "_groove_conf_scramble.png"
    if not (cog.exists() and scr.exists()):
        raise SystemExit("render the groove_conf PNGs first (see module docstring)")

    fig, (axa, axb) = plt.subplots(1, 2, figsize=(11, 4.6))
    for ax, path, lab in ((axa, cog, "a"), (axb, scr, "b")):
        ax.imshow(_crop(plt.imread(path)))
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        ax.text(0.0, 1.0, f"({lab})", transform=ax.transAxes, fontsize=13,
                fontweight="bold", va="top", ha="left", color=PALETTE["ink"])

    sm = ScalarMappable(norm=Normalize(PMIN, PMAX), cmap=PLDDT_CMAP)
    cbar = fig.colorbar(sm, ax=(axa, axb), orientation="horizontal",
                        fraction=0.05, pad=0.06, aspect=40)
    cbar.set_label("per-residue pLDDT", fontsize=9)
    cbar.set_ticks([50, 60, 70, 80, 90])
    save(fig, "figS2_groove_confidence")


if __name__ == "__main__":
    main()
