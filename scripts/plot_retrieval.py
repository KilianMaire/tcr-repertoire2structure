"""Fig 3: structure vs sequence retrieval, discovery and pre-registered held-out.

Reads paper/data/tcr_retrieval_top1.csv (reconstructed panel). Two panels:
  (a) discovery A*02:01: Top-1 per confidence readout vs the sequence baseline
      (0.0 by construction) and naive per-panel chance (0.25), with exact-binomial
      significance.
  (b) the pre-registered held-out A*11:01 outcome: the primary metric at 0.61 vs
      chance 0.5, not clearing the pre-committed p, so the confirmation is not
      licensed.

Usage: python scripts/plot_retrieval.py     # writes paper/figures/fig3_retrieval.png
"""
from __future__ import annotations
import csv, math, sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "paper/data"
C_DISC, C_HELD, C_CTRL = "#0072B2", "#D55E00", "#BBBBBB"
SURFACE = "#FCFCFB"
CONTROL = "iptm_groove_ctrl"


def binom_upper_p(k, n, p):
    """One-sided P(X >= k) for X ~ Binomial(n, p)."""
    return sum(math.comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(k, n + 1))


def _rows():
    with (DATA / "tcr_retrieval_top1.csv").open() as f:
        return [r for r in csv.DictReader(f) if r["panel"] == "reconstructed"]


def main(out="paper/figures/fig3_retrieval.png"):
    rows = _rows()
    disc = {r["readout"]: r for r in rows if r["run"] == "panel1"}
    held = {r["readout"]: r for r in rows if r["run"] == "hla_a1101"}

    fig, (axa, axb) = plt.subplots(1, 2, figsize=(12.5, 5), gridspec_kw={"width_ratios": [1.5, 1]})
    fig.patch.set_facecolor(SURFACE)
    for ax in (axa, axb):
        ax.set_facecolor(SURFACE)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)

    # (a) discovery battery, chance 0.25
    order = sorted(disc, key=lambda r: float(disc[r]["top1"]))
    ys = range(len(order))
    for y, ro in zip(ys, order):
        n = int(disc[ro]["n_tcr"]); top1 = float(disc[ro]["top1"]); k = round(top1 * n)
        p = binom_upper_p(k, n, 0.25)
        color = C_CTRL if ro == CONTROL else C_DISC
        axa.barh(y, top1, color=color, edgecolor=SURFACE, linewidth=1.2, height=0.7, zorder=2)
        tag = f"{top1:.2f}" + (" *" if p < 0.05 and ro != CONTROL else (" (ctrl)" if ro == CONTROL else " ns"))
        axa.annotate(tag, (top1 + 0.008, y), va="center", fontsize=8, color="#333")
    axa.axvline(0.25, color="#888", lw=1.4, ls="--", zorder=1)
    axa.annotate("naive chance 0.25", (0.25, len(order) - 0.3), fontsize=8, color="#888", ha="center")
    axa.axvline(0.0, color="#333", lw=1)
    axa.annotate("sequence\nbaseline 0.00", (0.01, 0.6), fontsize=8, color="#333", va="center")
    axa.set_yticks(list(ys))
    axa.set_yticklabels(order, fontsize=8.5)
    axa.set_xlim(0, 0.72)
    axa.set_xlabel("Top-1 retrieval (discovery A*02:01, reconstructed n=29)")
    axa.set_title("(a) discovery: confidence beats chance and sequence", fontsize=10.5, loc="left")

    # (b) held-out primary vs chance 0.5
    ro = "iptm_TCRpep_max"
    n = int(held[ro]["n_tcr"]); top1 = float(held[ro]["top1"]); k = round(top1 * n)
    p = binom_upper_p(k, n, 0.5)
    axb.bar(0, top1, width=0.5, color=C_HELD, edgecolor=SURFACE, linewidth=1.5, zorder=2)
    axb.axhline(0.5, color="#888", lw=1.4, ls="--", zorder=1)
    axb.annotate("chance 0.5", (0.32, 0.5), fontsize=8.5, color="#888", va="bottom")
    axb.annotate(f"{top1:.2f}  ({k}/{n})", (0, top1 + 0.02), ha="center", fontsize=10, color="#333")
    axb.annotate(f"binomial p = {p:.2f}: does not clear the\npre-registered 0.05, so the confirmation\nis not licensed",
                 (0, 0.88), ha="center", va="top", fontsize=9, color=C_HELD, fontweight="bold")
    axb.set_xticks([0]); axb.set_xticklabels(["iptm_TCRpep_max"], fontsize=9)
    axb.set_ylim(0, 1.0)
    axb.set_ylabel("Top-1 retrieval")
    axb.set_title("(b) pre-registered held-out A*11:01 (n=18)", fontsize=10.5, loc="left")

    fig.suptitle("Structure vs sequence: a discovery signal that does not confirm",
                 fontsize=13, fontweight="bold", x=0.02, ha="left", y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    print(f"wrote {out}")


if __name__ == "__main__":
    main(*(sys.argv[1:2] or []))
