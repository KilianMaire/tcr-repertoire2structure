"""Fig 3: structure vs sequence retrieval - a discovery signal that does not confirm.

Reads paper/data/tcr_retrieval_top1.csv (reconstructed panel only). Panels:
  (a) discovery battery A*02:01: Top-1 per confidence readout, against the
      sequence baseline (0.0 by construction) and naive per-panel chance (0.25).
  (b) pre-registered held-out A*11:01: the primary metric against chance (0.5).
  (c) negative controls: groove control and CDR3b-peptide contact, both < chance.
  (d) the gradation for the primary metric alone: discovery vs held-out.
  (e) the held-out null distribution, Binomial(18, 0.5), with 11/18 marked.

Significance stars are the exact one-sided binomial tail against per-panel chance.
All interpretation lives in the figure caption, not on the axes.

Usage: python scripts/plot_retrieval.py
"""
from __future__ import annotations
import csv, math, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import matplotlib.pyplot as plt
from figstyle import PALETTE, SURFACE, DATA, apply, despine, panel_label, save

CONTROL = "iptm_groove_ctrl"
PRIMARY = "iptm_TCRpep_max"
CONTACT_TOP1 = 0.19  # CDR3b-peptide contact, refuted, not in the CSV (see CANONICAL_NUMBERS.md)


def binom_upper_p(k, n, p):
    return sum(math.comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(k, n + 1))


def binom_pmf(k, n, p):
    return math.comb(n, k) * p ** k * (1 - p) ** (n - k)


def _rows():
    with (DATA / "tcr_retrieval_top1.csv").open() as f:
        return [r for r in csv.DictReader(f) if r["panel"] == "reconstructed"]


def main():
    apply()
    rows = _rows()
    disc = {r["readout"]: r for r in rows if r["run"] == "panel1"}
    held = {r["readout"]: r for r in rows if r["run"] == "hla_a1101"}

    fig = plt.figure(figsize=(13, 8))
    gs = fig.add_gridspec(
        2, 3, width_ratios=[1.55, 1, 1], height_ratios=[1, 1],
        wspace=0.62, hspace=0.62, left=0.09, right=0.975, top=0.95, bottom=0.09,
    )
    axa = fig.add_subplot(gs[:, 0])
    axb = fig.add_subplot(gs[0, 1])
    axc = fig.add_subplot(gs[0, 2])
    axd = fig.add_subplot(gs[1, 1])
    axe = fig.add_subplot(gs[1, 2])
    for ax in (axa, axb, axc, axd, axe):
        despine(ax)

    # ---- (a) discovery battery against chance 0.25 --------------------------
    order = sorted(disc, key=lambda r: float(disc[r]["top1"]))
    for y, ro in enumerate(order):
        n = int(disc[ro]["n_tcr"]); top1 = float(disc[ro]["top1"]); k = round(top1 * n)
        p = binom_upper_p(k, n, 0.25)
        is_ctrl = ro == CONTROL
        color = PALETTE["grey"] if is_ctrl else PALETTE["blue"]
        axa.barh(y, top1, color=color, edgecolor=SURFACE, linewidth=1.2, height=0.72, zorder=2)
        star = " *" if (p < 0.05 and not is_ctrl) else ""
        axa.annotate(f"{top1:.2f}{star}", (top1 + 0.008, y), va="center",
                     fontsize=8.5, color=PALETTE["ink"])
    axa.axvline(0.25, color=PALETTE["mute"], lw=1.3, ls="--", zorder=1)
    axa.annotate("chance 0.25", (0.25, len(order) + 0.15), fontsize=8.5,
                 color=PALETTE["mute"], ha="center", va="bottom")
    axa.axvline(0.0, color=PALETTE["ink"], lw=1)
    axa.annotate("sequence 0.00", (0.008, -0.75), fontsize=8.5, color=PALETTE["ink"],
                 va="center", ha="left")
    axa.set_yticks(range(len(order)))
    axa.set_yticklabels(order, fontsize=9)
    axa.set_xlim(0, 0.76)
    axa.set_ylim(-1.1, len(order) + 0.55)
    axa.set_xlabel("Top-1 retrieval  (discovery A*02:01, n=29)")
    panel_label(axa, "a", x=-0.40)

    # ---- (b) held-out primary against chance 0.5 -----------------------------
    ro = PRIMARY
    n = int(held[ro]["n_tcr"]); top1 = float(held[ro]["top1"]); k = round(top1 * n)
    axb.bar(0, top1, width=0.5, color=PALETTE["orange"], edgecolor=SURFACE,
            linewidth=1.5, zorder=2)
    axb.axhline(0.5, color=PALETTE["mute"], lw=1.3, ls="--", zorder=1)
    axb.annotate("chance 0.5", (0.34, 0.5), fontsize=8, color=PALETTE["mute"], va="bottom")
    axb.annotate(f"{top1:.2f}  ({k}/{n})", (0, top1 + 0.03), ha="center",
                 fontsize=9.5, color=PALETTE["ink"])
    axb.set_xticks([0]); axb.set_xticklabels(["iptm_TCRpep_max"], fontsize=8)
    axb.set_xlim(-0.6, 0.6)
    axb.set_ylim(0, 1.0)
    axb.set_ylabel("Top-1  (held-out A*11:01, n=18)", fontsize=8.5)
    panel_label(axb, "b", x=-0.30)

    # ---- (c) negative controls, both below chance 0.25 -----------------------
    ctrl_top1 = float(disc[CONTROL]["top1"])
    labels_c = ["groove\ncontrol", "CDR3b-pep\ncontact"]
    vals_c = [ctrl_top1, CONTACT_TOP1]
    xs = [0, 1]
    axc.bar(xs, vals_c, width=0.55, color=PALETTE["grey"], edgecolor=SURFACE,
            linewidth=1.5, zorder=2)
    for x, v in zip(xs, vals_c):
        axc.annotate(f"{v:.2f}", (x, v + 0.02), ha="center", fontsize=9, color=PALETTE["ink"])
    axc.axhline(0.25, color=PALETTE["mute"], lw=1.3, ls="--", zorder=1)
    axc.annotate("chance 0.25", (1.58, 0.263), fontsize=8, color=PALETTE["mute"],
                 va="bottom", ha="right")
    axc.set_xticks(xs); axc.set_xticklabels(labels_c, fontsize=8)
    axc.set_xlim(-0.6, 1.6)
    axc.set_ylim(0, 0.45)
    axc.set_ylabel("Top-1  (discovery A*02:01)", fontsize=8.5)
    panel_label(axc, "c", x=-0.30)

    # ---- (d) gradation for the primary metric alone ---------------------------
    n_d = int(disc[PRIMARY]["n_tcr"]); top1_d = float(disc[PRIMARY]["top1"]); k_d = round(top1_d * n_d)
    p_d = binom_upper_p(k_d, n_d, 0.25)
    n_h = int(held[PRIMARY]["n_tcr"]); top1_h = float(held[PRIMARY]["top1"]); k_h = round(top1_h * n_h)
    p_h = binom_upper_p(k_h, n_h, 0.5)

    xs_d = [0, 1]
    axd.bar(xs_d, [top1_d, top1_h], width=0.55,
            color=[PALETTE["blue"], PALETTE["orange"]], edgecolor=SURFACE,
            linewidth=1.5, zorder=2)
    axd.plot([-0.28, 0.28], [0.25, 0.25], color=PALETTE["mute"], lw=1.3, ls="--", zorder=1)
    axd.plot([0.72, 1.28], [0.5, 0.5], color=PALETTE["mute"], lw=1.3, ls="--", zorder=1)
    axd.annotate(f"{top1_d:.2f}", (0, top1_d + 0.055), ha="center", fontsize=9, color=PALETTE["ink"])
    axd.annotate(f"p={p_d:.1e}", (0, top1_d + 0.015), ha="center", fontsize=7.5, color=PALETTE["mute"])
    axd.annotate(f"{top1_h:.2f}", (1, top1_h + 0.055), ha="center", fontsize=9, color=PALETTE["ink"])
    axd.annotate(f"p={p_h:.2f}", (1, top1_h + 0.015), ha="center", fontsize=7.5, color=PALETTE["mute"])
    axd.set_xticks(xs_d)
    axd.set_xticklabels([f"discovery\nchance 0.25\nn={n_d}", f"held-out\nchance 0.5\nn={n_h}"], fontsize=7.7)
    axd.set_xlim(-0.6, 1.6)
    axd.set_ylim(0, 0.85)
    axd.set_ylabel("Top-1  (iptm_TCRpep_max)", fontsize=8.5)
    panel_label(axd, "d", x=-0.30)

    # ---- (e) held-out null distribution, Binomial(18, 0.5) --------------------
    n_e, p_e, obs = 18, 0.5, 11
    ks = list(range(n_e + 1))
    pm = [binom_pmf(k, n_e, p_e) for k in ks]
    colors_e = [PALETTE["orange"] if k >= obs else PALETTE["grey"] for k in ks]
    axe.bar(ks, pm, color=colors_e, edgecolor=SURFACE, linewidth=0.6, width=0.8, zorder=2)
    axe.axvline(obs, color=PALETTE["ink"], lw=1.2, zorder=3)
    axe.annotate(f"observed\n{obs}/18", (obs, max(pm) * 1.02), ha="center", va="bottom",
                 fontsize=7.5, color=PALETTE["ink"])
    axe.set_xlim(-0.8, n_e + 0.8)
    axe.set_ylim(0, max(pm) * 1.32)
    axe.set_xticks([0, 3, 6, 9, 12, 15, 18])
    axe.set_xlabel("correct out of 18", fontsize=8.5)
    axe.set_ylabel("P(k)", fontsize=8.5)
    panel_label(axe, "e", x=-0.30)

    save(fig, "fig3_retrieval")


if __name__ == "__main__":
    main()
