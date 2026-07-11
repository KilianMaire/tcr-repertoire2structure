"""Figure: structural confidence judges MHC-peptide presentation.

Read-only. Reuses analyze_mhc_scramble.py. Grouped horizontal bars of the
binder-vs-scramble AUROC per metric, one series per HLA, with a chance line at
0.5 and TCR-bootstrap CIs. Okabe-Ito 2-color palette (validated CVD-safe).

Usage: python scripts/plot_mhc_scramble.py           # writes docs/mhc_peptide_presentation.png
"""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze_mhc_scramble import METRICS, summarize

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C_A0201, C_A1101 = "#0072B2", "#D55E00"
SURFACE = "#FCFCFB"
RUNS = [("A*02:01 (n=29)", "runs/panel1", C_A0201),
        ("A*11:01 (n=18)", "runs/hla_a1101", C_A1101)]


def main(out_path="docs/mhc_peptide_presentation.png"):
    scored = {label: {m: summarize(run_dir, fn, reconstructed_only=True)
                      for m, fn in METRICS.items()}
              for label, run_dir, _ in RUNS}
    # order metrics by mean AUROC across the two runs, best at top
    order = sorted(METRICS, key=lambda m: -sum(scored[l][m]["auroc"] for l, _, _ in RUNS))

    fig, ax = plt.subplots(figsize=(9, 5.2))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    h = 0.36
    ys = list(range(len(order)))[::-1]  # top metric highest on the axis
    for si, (label, _run, color) in enumerate(RUNS):
        offset = (si - 0.5) * h
        for y, m in zip(ys, order):
            s = scored[label][m]
            au, (lo, hi) = s["auroc"], s["auroc_ci"]
            ax.barh(y + offset, au - 0.5, left=0.5, height=h, color=color,
                    edgecolor=SURFACE, linewidth=1.2, zorder=2,
                    label=label if y == ys[0] else None)
            ax.plot([lo, hi], [y + offset, y + offset], color="#333333", lw=1.4, zorder=3)
            ax.annotate(f"{au:.2f}", xy=(max(au, hi) + 0.006, y + offset), va="center",
                        ha="left", fontsize=7.5, color="#333333")
    ax.axvline(0.5, color="#888888", lw=1.2, zorder=1)
    ax.annotate("chance", xy=(0.5, len(order) - 0.4), xytext=(0.503, len(order) - 0.35),
                fontsize=8, color="#888888", va="top")
    ax.set_yticks(ys)
    ax.set_yticklabels(order, fontsize=9)
    ax.set_xlim(0.45, 1.05)
    ax.set_xlabel("AUROC: genuine binder (cognate + decoys) vs composition-scramble")
    ax.set_title("the mirror image of the TCR result: the groove interface confidence "
                 "separates a real ligand from its scramble",
                 fontsize=8.8, color="#555555", loc="left", pad=8)
    fig.suptitle("Structural confidence judges MHC-peptide presentation",
                 fontsize=12.5, fontweight="bold", x=0.02, ha="left", y=0.98)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(fontsize=9, loc="lower right", frameon=False)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, facecolor=SURFACE)
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main(*(sys.argv[1:2] or []))
