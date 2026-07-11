"""Fig 1 reader aid: annotated cartoon (ribbon) views of a predicted TCR-pMHC complex.

Composes the two PyMOL ray-traced renders (produced by
scripts/render_structure_pymol.py in the pymol micromamba env) into a single
labelled figure with a lay-reader caption and a colour legend. If the render PNGs
are missing it invokes PyMOL to make them.

Usage: python scripts/plot_structure_3d.py     # writes paper/figures/fig1_structure.png
"""
from __future__ import annotations
import subprocess, sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parent.parent
FIGS = ROOT / "paper/figures"
SURFACE = "#FCFCFB"
DEFAULT_CIF = ("runs/panel1/folds/b97fd808da3f__GILGFVFTL/b97fd808da3f__GILGFVFTL"
               "/seed_101/predictions/b97fd808da3f__GILGFVFTL_sample_1.cif")
LEGEND = [("#0072B2", "TCR α chain (recognises)"), ("#009E73", "TCR β chain (recognises)"),
          ("#E69F00", "MHC class I (presents)"), ("#CC79A7", "β₂-microglobulin"),
          ("#D55E00", "peptide antigen")]


def _autocrop(img, pad=12):
    """Trim near-white margins from a ray-traced render."""
    rgb = img[..., :3] if img.shape[-1] == 4 else img
    mask = (rgb < 0.97).any(-1)
    if not mask.any():
        return img
    ys, xs = np.where(mask)
    y0, y1 = max(ys.min() - pad, 0), min(ys.max() + pad, img.shape[0])
    x0, x1 = max(xs.min() - pad, 0), min(xs.max() + pad, img.shape[1])
    return img[y0:y1, x0:x1]


def _ensure_renders(cif):
    v1, v2 = FIGS / "_struct_view1.png", FIGS / "_struct_view2.png"
    if v1.exists() and v2.exists():
        return v1, v2
    mm = Path.home() / ".local/bin/micromamba"
    subprocess.run([str(mm), "run", "-n", "pymol", "python",
                    "scripts/render_structure_pymol.py", cif, str(FIGS)],
                   cwd=ROOT, check=True, env={"MAMBA_ROOT_PREFIX": str(Path.home() / "micromamba"),
                                              "PATH": "/usr/bin:/bin:/usr/sbin:/sbin"})
    return v1, v2


def main(cif=DEFAULT_CIF, out="paper/figures/fig1_structure.png"):
    v1, v2 = _ensure_renders(cif)
    fig = plt.figure(figsize=(12, 6.2))
    fig.patch.set_facecolor(SURFACE)
    for i, (vp, title) in enumerate([(v1, "the whole complex"),
                                     (v2, "looking down onto the groove")]):
        ax = fig.add_subplot(1, 2, i + 1)
        ax.imshow(_autocrop(plt.imread(vp)))
        ax.set_title(title, fontsize=10.5, color="#333")
        ax.axis("off")
    fig.legend(handles=[Patch(color=c, label=l) for c, l in LEGEND],
               loc="lower center", ncol=5, frameon=False, fontsize=9, bbox_to_anchor=(0.5, 0.0))
    fig.suptitle("A predicted TCR-pMHC complex (flu GILGFVFTL on HLA-A*02:01)",
                 fontsize=13, fontweight="bold", y=0.99)
    fig.text(0.5, 0.925, "the T cell receptor (blue, green) docks on top of the peptide (red) held in "
             "the groove of the MHC (orange); this is how a T cell reads an antigen",
             ha="center", fontsize=8.8, color="#666")
    fig.tight_layout(rect=[0, 0.06, 1, 0.9])
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    print(f"wrote {out}")


if __name__ == "__main__":
    main(*(sys.argv[1:]))
