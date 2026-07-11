"""Figure: structural confidence separates TCRs, not epitopes.

Read-only. Reuses the decomposition in analyze_confidence_variance.py and renders
a two panel figure:

  (a) stacked variance bars (TCR, peptide, cognate, residual) per condition,
      showing the TCR docking term dominates;
  (b) the cognate effect size (cognate minus panel mean, ipTM units) with a
      bootstrap CI, showing it is significant on discovery and null on held-out.

Palette is Okabe-Ito (colorblind safe; validated with the dataviz palette
checker). Secondary encoding: direct labels and 2px surface gaps between fills.

Usage: python scripts/plot_confidence_variance.py            # writes docs/confidence_variance.png
"""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze_confidence_variance import (READOUTS, load_panel,
                                         variance_decomposition, cognate_effect)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Okabe-Ito, validated CVD-safe for the meaningful segments
C_TCR, C_PEP, C_COG, C_RESID = "#0072B2", "#009E73", "#D55E00", "#BBBBBB"
SURFACE = "#FCFCFB"

# conditions: (label, run_dir, readout_name)
CONDS = [
    ("A*02:01\nTCR-pep\n(discovery)", "runs/panel1", "iptm_TCRpep_max"),
    ("A*11:01\nTCR-pep\n(held-out)", "runs/hla_a1101", "iptm_TCRpep_max"),
    ("A*02:01\ngroove\n(control)", "runs/panel1", "iptm_groove_ctrl"),
    ("A*11:01\ngroove\n(control)", "runs/hla_a1101", "iptm_groove_ctrl"),
]


def gather():
    rows_out = []
    for label, run_dir, readout in CONDS:
        rows = load_panel(run_dir, READOUTS[readout])
        vd = variance_decomposition(rows)
        ce = cognate_effect(rows)
        resid = max(0.0, 1 - vd["tcr"] - vd["cognate_within_tcr"] - vd["peptide_within_tcr"])
        rows_out.append({"label": label, "vd": vd, "resid": resid, "ce": ce})
    return rows_out


def main(out_path="docs/confidence_variance.png"):
    data = gather()
    x = list(range(len(data)))
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(11, 5), gridspec_kw={"width_ratios": [1.25, 1]})
    fig.patch.set_facecolor(SURFACE)

    # (a) stacked variance
    axa.set_facecolor(SURFACE)
    gap = 0.012  # 2px-ish surface gap between fills, in axis (fraction) units
    for i, d in enumerate(data):
        vd = d["vd"]
        segs = [(vd["tcr"], C_TCR, "TCR"), (vd["peptide_within_tcr"], C_PEP, "peptide"),
                (vd["cognate_within_tcr"], C_COG, "cognate"), (d["resid"], C_RESID, "residual")]
        bottom = 0.0
        for frac, color, _ in segs:
            axa.bar(i, frac, bottom=bottom, width=0.62, color=color,
                    edgecolor=SURFACE, linewidth=1.5)
            bottom += frac + gap
        # highlight the cognate fraction with a direct label
        cog = vd["cognate_within_tcr"]
        axa.annotate(f"cognate {cog*100:.1f}%", xy=(i, vd["tcr"] + vd["peptide_within_tcr"] + cog + 2*gap),
                     xytext=(i, min(0.97, vd["tcr"] + vd["peptide_within_tcr"] + cog + 0.10)),
                     ha="center", va="bottom", fontsize=8.5, color=C_COG, fontweight="bold",
                     arrowprops=dict(arrowstyle="-", color=C_COG, lw=1))
    axa.set_xticks(x)
    axa.set_xticklabels([d["label"] for d in data], fontsize=8.5)
    axa.set_ylabel("fraction of readout variance")
    axa.set_ylim(0, 1.12)
    axa.set_title("(a) where the confidence variance goes", fontsize=11, loc="left")
    for s in ("top", "right"):
        axa.spines[s].set_visible(False)
    # legend proxies
    from matplotlib.patches import Patch
    axa.legend(handles=[Patch(color=C_TCR, label="TCR identity"),
                        Patch(color=C_PEP, label="peptide identity"),
                        Patch(color=C_COG, label="cognate status"),
                        Patch(color=C_RESID, label="residual")],
               fontsize=8, loc="upper right", frameon=False, ncol=2)

    # (b) cognate effect size with bootstrap CI
    axb.set_facecolor(SURFACE)
    axb.axhline(0, color="#888888", lw=1, zorder=1)
    for i, d in enumerate(data):
        ce = d["ce"]
        lo, hi = ce["ci"]
        signif = lo > 0  # CI excludes zero
        color = C_COG if signif else "#999999"
        axb.errorbar(i, ce["mean_delta"], yerr=[[ce["mean_delta"] - lo], [hi - ce["mean_delta"]]],
                     fmt="o", ms=8, color=color, ecolor=color, elinewidth=2, capsize=4, zorder=3)
        tag = f"p={ce['perm_p']:.3g}" + ("  *" if signif else "  ns")
        axb.annotate(tag, xy=(i, hi), xytext=(i, hi + 0.006), ha="center", va="bottom",
                     fontsize=8, color=color, fontweight="bold" if signif else "normal")
    axb.set_xticks(x)
    axb.set_xticklabels([d["label"] for d in data], fontsize=8.5)
    axb.set_ylabel("cognate minus panel mean (ipTM)")
    axb.set_title("(b) peptide-specificity effect size", fontsize=11, loc="left")
    for s in ("top", "right"):
        axb.spines[s].set_visible(False)

    fig.suptitle("Structural confidence separates TCRs, not epitopes",
                 fontsize=13, fontweight="bold", x=0.5)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, facecolor=SURFACE)
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main(*(sys.argv[1:2] or []))
