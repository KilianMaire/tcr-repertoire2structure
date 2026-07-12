"""Fig 4: where the structural-confidence variance goes.

Reads the committed paper/data/confidence_variance.csv (reconstructed TCRs only;
poly-G stub folds excluded upstream). Panels:
  (a) sequential variance decomposition (TCR identity, peptide identity, cognate
      status, residual) per condition.
  (b) the cognate-status effect (cognate minus panel mean, ipTM units) with its
      bootstrap 95% CI; significance is read off whether the CI clears zero.
  (c) ICC per condition: the between-TCR share of variance.
  (d) cognate variance fraction per condition: the recognition term is small
      everywhere and does not replicate held-out.

No title, no interpretation on the axes. Reproduces from committed data, no GPU.

Usage: python scripts/plot_confidence_variance.py
"""
from __future__ import annotations
import csv, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from figstyle import PALETTE, SURFACE, DATA, apply, despine, panel_label, save

C_TCR, C_PEP, C_COG, C_RESID = PALETTE["blue"], PALETTE["green"], PALETTE["orange"], PALETTE["grey"]

# (label, run, readout) in display order
CONDS = [
    ("A*02:01\nTCR-pep\n(discovery)", "panel1", "iptm_TCRpep_max"),
    ("A*11:01\nTCR-pep\n(held-out)", "hla_a1101", "iptm_TCRpep_max"),
    ("A*02:01\ngroove\n(control)", "panel1", "iptm_groove_ctrl"),
    ("A*11:01\ngroove\n(control)", "hla_a1101", "iptm_groove_ctrl"),
]


def _load():
    with (DATA / "confidence_variance.csv").open() as f:
        idx = {(r["run"], r["readout"]): r for r in csv.DictReader(f)}
    out = []
    for label, run, ro in CONDS:
        r = idx[(run, ro)]
        vt, vp, vc = float(r["var_tcr"]), float(r["var_peptide"]), float(r["var_cognate"])
        out.append({
            "label": label, "tcr": vt, "pep": vp, "cog": vc,
            "resid": max(0.0, 1 - vt - vp - vc),
            "delta": float(r["cognate_delta"]),
            "lo": float(r["delta_ci_lo"]), "hi": float(r["delta_ci_hi"]),
            "icc": float(r["icc"]),
            "n": int(r["n_tcr"]),
        })
    return out


def main():
    apply()
    data = _load()
    x = range(len(data))
    labels = [d["label"] for d in data]

    fig = plt.figure(figsize=(12.5, 9))
    gs = fig.add_gridspec(2, 2, hspace=0.55, wspace=0.32,
                           left=0.07, right=0.98, top=0.95, bottom=0.13)
    axa = fig.add_subplot(gs[0, 0])
    axb = fig.add_subplot(gs[0, 1])
    axc = fig.add_subplot(gs[1, 0])
    axd = fig.add_subplot(gs[1, 1])

    # (a) stacked variance decomposition; cognate % on a clean row above the bars
    gap = 0.012
    for i, d in enumerate(data):
        bottom = 0.0
        for frac, color in ((d["tcr"], C_TCR), (d["pep"], C_PEP),
                            (d["cog"], C_COG), (d["resid"], C_RESID)):
            axa.bar(i, frac, bottom=bottom, width=0.62, color=color,
                    edgecolor=SURFACE, linewidth=1.5)
            bottom += frac + gap
        axa.annotate(f"cognate {d['cog']*100:.1f}%", (i, 1.10), ha="center", va="bottom",
                     fontsize=7.8, color=C_COG, fontweight="bold")
    despine(axa)
    axa.set_xticks(list(x))
    axa.set_xticklabels(labels, fontsize=8)
    axa.set_ylabel("fraction of readout variance")
    axa.set_ylim(0, 1.28)
    axa.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    panel_label(axa, "a", x=-0.15)

    # (b) cognate effect size with bootstrap CI
    axb.axhline(0, color=PALETTE["mute"], lw=1, zorder=1)
    for i, d in enumerate(data):
        signif = d["lo"] > 0
        color = C_COG if signif else PALETTE["grey"]
        axb.errorbar(i, d["delta"], yerr=[[d["delta"] - d["lo"]], [d["hi"] - d["delta"]]],
                     fmt="o", ms=8, color=color, ecolor=color, elinewidth=2, capsize=4, zorder=3)
    despine(axb)
    axb.set_xticks(list(x))
    axb.set_xticklabels(labels, fontsize=8)
    axb.set_ylabel("cognate minus panel mean (ipTM)")
    axb.set_xlim(-0.5, len(data) - 0.5)
    panel_label(axb, "b", x=-0.2)

    # (c) ICC per condition: the between-TCR share of variance
    for i, d in enumerate(data):
        color = C_TCR if d["icc"] > 0 else PALETTE["grey"]
        axc.bar(i, d["icc"], width=0.55, color=color, edgecolor=SURFACE, linewidth=1.2)
    axc.axhline(0, color=PALETTE["ink"], lw=0.9, zorder=1)
    despine(axc)
    axc.set_xticks(list(x))
    axc.set_xticklabels(labels, fontsize=8)
    axc.set_ylabel("ICC (between-TCR variance share)")
    axc.set_xlim(-0.5, len(data) - 0.5)
    panel_label(axc, "c", x=-0.2)

    # (d) cognate variance fraction per condition: small everywhere,
    # does not replicate held-out
    for i, d in enumerate(data):
        signif = d["lo"] > 0
        color = C_COG if signif else PALETTE["grey"]
        axd.bar(i, d["cog"] * 100, width=0.55, color=color, edgecolor=SURFACE, linewidth=1.2)
        axd.annotate(f"{d['cog']*100:.1f}%", (i, d["cog"] * 100), xytext=(0, 3),
                     textcoords="offset points", ha="center", va="bottom", fontsize=8)
    despine(axd)
    axd.set_xticks(list(x))
    axd.set_xticklabels(labels, fontsize=8)
    axd.set_ylabel("cognate-status variance fraction (%)")
    axd.set_ylim(0, max(d["cog"] for d in data) * 100 * 1.35)
    axd.set_xlim(-0.5, len(data) - 0.5)
    panel_label(axd, "d", x=-0.2)

    fig.legend(handles=[Patch(color=C_TCR, label="TCR identity"),
                        Patch(color=C_PEP, label="peptide identity"),
                        Patch(color=C_COG, label="cognate status"),
                        Patch(color=C_RESID, label="residual")],
               loc="lower center", ncol=4, frameon=False, fontsize=9.5,
               bbox_to_anchor=(0.5, 0.01))
    save(fig, "fig4_confidence_variance")


if __name__ == "__main__":
    main()
