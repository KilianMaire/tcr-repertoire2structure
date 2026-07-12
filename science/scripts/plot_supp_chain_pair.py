"""Fig S3: chain-pair ipTM matrix for one confident cognate A*11:01 complex.

Reads paper/data/chain_pair_iptm_example.csv (5x5 mean chain-pair interface
ipTM, chains TCRa/TCRb/MHC/b2m/peptide). Panel (a) is the annotated heatmap;
the two cells behind iptm_TCRpep (TCRa-peptide, TCRb-peptide) are outlined in
orange and the one cell behind iptm_groove (MHC-peptide) is outlined in
green, matching the readouts used in Fig 3 / Fig 5. No result text is drawn
on the figure; the mapping lives in a small legend and, optionally, panel (b).

Usage: python scripts/plot_supp_chain_pair.py
"""
from __future__ import annotations
import csv, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from figstyle import PALETTE, SURFACE, DATA, apply, panel_label, save

CHAINS = ["TCRa", "TCRb", "MHC", "b2m", "peptide"]

# cells used by paper readouts, as (row_chain, col_chain)
TCRPEP_CELLS = [("TCRa", "peptide"), ("TCRb", "peptide")]
GROOVE_CELLS = [("MHC", "peptide")]


def _load_matrix():
    m = np.zeros((len(CHAINS), len(CHAINS)))
    idx = {c: i for i, c in enumerate(CHAINS)}
    with (DATA / "chain_pair_iptm_example.csv").open() as f:
        for row in csv.DictReader(f):
            i, j = idx[row["row_chain"]], idx[row["col_chain"]]
            m[i, j] = float(row["iptm"])
    return m


def _outline(ax, i, j, color):
    ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                            edgecolor=color, linewidth=2.4, zorder=5))


def main():
    apply()
    m = _load_matrix()
    idx = {c: i for i, c in enumerate(CHAINS)}

    seq_cmap = LinearSegmentedColormap.from_list(
        "blue_seq", ["#F2F8FC", PALETTE["blue"]], N=256
    )

    fig, ax = plt.subplots(figsize=(7, 6))
    fig.subplots_adjust(left=0.16, right=0.86, top=0.90, bottom=0.14)

    im = ax.imshow(m, cmap=seq_cmap, vmin=0, vmax=1, aspect="equal", zorder=1)

    ax.set_xticks(range(len(CHAINS)))
    ax.set_yticks(range(len(CHAINS)))
    ax.set_xticklabels(CHAINS, fontsize=9.5)
    ax.set_yticklabels(CHAINS, fontsize=9.5)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)

    # annotate values, contrast-aware text color
    for i in range(len(CHAINS)):
        for j in range(len(CHAINS)):
            v = m[i, j]
            color = "white" if v > 0.55 else PALETTE["ink"]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                     fontsize=8.5, color=color, zorder=3)

    for r, c in TCRPEP_CELLS:
        _outline(ax, idx[r], idx[c], PALETTE["orange"])
    for r, c in GROOVE_CELLS:
        _outline(ax, idx[r], idx[c], PALETTE["green"])

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("chain-pair ipTM", fontsize=9.5)
    cbar.ax.tick_params(labelsize=8.5)

    legend_handles = [
        Line2D([0], [0], color=PALETTE["orange"], lw=2.4, label="iptm_TCRpep"),
        Line2D([0], [0], color=PALETTE["green"], lw=2.4, label="iptm_groove"),
    ]
    ax.legend(handles=legend_handles, loc="upper left", bbox_to_anchor=(1.32, 1.0),
               frameon=False, fontsize=8.5, handlelength=1.6, borderaxespad=0)

    panel_label(ax, "a", x=-0.16, y=1.05)

    save(fig, "figS3_chain_pair_iptm")


if __name__ == "__main__":
    main()
