"""Fig 5: confidence reads MHC-peptide presentation, independent of TCR.

Reads paper/data/mhc_presentation.csv, paper/data/peptide_plddt.csv, and
paper/data/scramble_anchor_permissiveness.csv (all reconstructed-panel, all
binder-vs-scramble). Panels:
  (a) binder-vs-scramble AUROC per confidence readout, grouped by HLA, with
      95 percent CI whiskers and a chance line at 0.5.
  (b)/(c) rendered grooves (cognate vs scramble) for one A*11:01 clonotype,
      peptide colored by per-residue pLDDT, sharing one colorbar.
  (d) per-residue pLDDT along the same peptide, cognate vs scramble.
  (e) anchor-residue permissiveness, cognate vs scramble, both HLAs.

All interpretation lives in the figure caption, not on the axes.

Usage: python scripts/plot_mhc_scramble.py
"""
from __future__ import annotations
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from figstyle import PALETTE, SURFACE, DATA, FIGS, apply, despine, panel_label, save

PLDDT_CMAP = LinearSegmentedColormap.from_list("plddt_rwb", ["#D7191C", "#FFFFFF", "#2C7BB6"])
PMIN, PMAX = 50, 95

HLA_COLOR = {"A*02:01": PALETTE["blue"], "A*11:01": PALETTE["orange"]}

METRIC_LABEL = {
    "neg_gpde_groove": "neg_gpde_groove",
    "iptm_b2m_pep": "iptm_b2m_pep",
    "iptm_groove": "iptm_groove",
    "pep_ptm": "pep_ptm",
    "pep_plddt": "pep_plddt",
    "pep_iptm": "pep_iptm",
    "ranking_score": "ranking_score",
}


def _crop(img, pad=6):
    rgb = img[..., :3]
    mask = (rgb < 0.97).any(-1)
    ys, xs = np.where(mask)
    return img[max(ys.min() - pad, 0):ys.max() + pad, max(xs.min() - pad, 0):xs.max() + pad]


def _read_presentation():
    with (DATA / "mhc_presentation.csv").open() as f:
        return list(csv.DictReader(f))


def _read_plddt():
    with (DATA / "peptide_plddt.csv").open() as f:
        rows = list(csv.DictReader(f))
    cog = sorted([r for r in rows if r["construct"] == "cognate"], key=lambda r: int(r["pos"]))
    scr = sorted([r for r in rows if r["construct"] == "scramble"], key=lambda r: int(r["pos"]))
    return cog, scr


def _read_anchor():
    with (DATA / "scramble_anchor_permissiveness.csv").open() as f:
        return list(csv.DictReader(f))


def main():
    apply()
    pres = _read_presentation()
    cog, scr = _read_plddt()
    anchor = _read_anchor()

    fig = plt.figure(figsize=(13, 9))
    gs = fig.add_gridspec(
        2, 3, width_ratios=[1.4, 1, 1], height_ratios=[1, 1.05],
        wspace=0.55, hspace=0.30, left=0.08, right=0.98, top=0.97, bottom=0.08,
    )
    axa = fig.add_subplot(gs[:, 0])
    axb = fig.add_subplot(gs[0, 1])
    axc = fig.add_subplot(gs[0, 2])
    axd = fig.add_subplot(gs[1, 1])
    axe = fig.add_subplot(gs[1, 2])
    despine(axa); despine(axd); despine(axe)

    # ---- (a) binder-vs-scramble AUROC, grouped bars per HLA -------------------
    metrics = list(METRIC_LABEL)

    def mean_auroc(m):
        vals = [float(r["auroc"]) for r in pres if r["metric"] == m]
        return sum(vals) / len(vals)

    # ascending order -> plotted bottom to top, so best metric lands at the top
    order = sorted(metrics, key=mean_auroc)

    hlas = ["A*02:01", "A*11:01"]
    bar_h = 0.36
    for y, m in enumerate(order):
        for j, hla in enumerate(hlas):
            row = next(r for r in pres if r["metric"] == m and r["hla"] == hla)
            auroc = float(row["auroc"])
            lo = float(row["auroc_ci_lo"])
            hi = float(row["auroc_ci_hi"])
            yy = y + (j - 0.5) * bar_h
            color = HLA_COLOR[hla]
            axa.barh(yy, auroc, height=bar_h * 0.92, color=color, edgecolor=SURFACE,
                     linewidth=1.0, zorder=2)
            axa.plot([lo, hi], [yy, yy], color=PALETTE["ink"], lw=1.1, zorder=3)
            axa.plot([lo, lo], [yy - 0.05, yy + 0.05], color=PALETTE["ink"], lw=1.1, zorder=3)
            axa.plot([hi, hi], [yy - 0.05, yy + 0.05], color=PALETTE["ink"], lw=1.1, zorder=3)
            axa.annotate(f"{auroc:.2f}", (hi + 0.012, yy), va="center", fontsize=7.6,
                         color=PALETTE["ink"])
    axa.axvline(0.5, color=PALETTE["mute"], lw=1.3, ls="--", zorder=1)
    axa.annotate("chance 0.5", (0.5, len(order) - 0.55), fontsize=8.5,
                 color=PALETTE["mute"], ha="center", va="bottom")
    axa.set_yticks(range(len(order)))
    axa.set_yticklabels([METRIC_LABEL[m] for m in order], fontsize=9)
    axa.set_xlim(0.35, 1.12)
    axa.set_ylim(-0.6, len(order) - 0.15)
    axa.set_xlabel("binder-vs-scramble AUROC")
    handles = [plt.Rectangle((0, 0), 1, 1, color=HLA_COLOR[h]) for h in hlas]
    axa.legend(handles, hlas, loc="lower right", frameon=False, fontsize=8.5,
               handlelength=1.1, handleheight=1.1)
    panel_label(axa, "a", x=-0.30)

    # ---- (b)/(c) rendered grooves, shared pLDDT colorbar -----------------------
    cog_png = FIGS / "_groove_conf_cognate.png"
    scr_png = FIGS / "_groove_conf_scramble.png"
    if not (cog_png.exists() and scr_png.exists()):
        raise SystemExit("missing rendered groove_conf PNGs in paper/figures/")
    for ax, path, tag, lab in ((axb, cog_png, "cognate", "b"), (axc, scr_png, "scramble", "c")):
        ax.imshow(_crop(plt.imread(path)))
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        panel_label(ax, lab, x=-0.02, y=1.02)
        ax.text(0.5, -0.02, tag, transform=ax.transAxes, fontsize=9,
                ha="center", va="top", color=PALETTE["ink"])

    sm = ScalarMappable(norm=Normalize(PMIN, PMAX), cmap=PLDDT_CMAP)
    cbar = fig.colorbar(sm, ax=(axb, axc), orientation="horizontal",
                        fraction=0.05, pad=0.10, aspect=42)
    cbar.set_label("per-residue pLDDT", fontsize=9)
    cbar.set_ticks([50, 60, 70, 80, 90])

    # ---- (d) per-residue pLDDT along the peptide -------------------------------
    xs = [int(r["pos"]) for r in cog]
    y_cog = [float(r["plddt"]) for r in cog]
    y_scr = [float(r["plddt"]) for r in scr]
    axd.plot(xs, y_cog, marker="o", ms=5, color=PALETTE["blue"], lw=1.8, label="cognate", zorder=3)
    axd.plot(xs, y_scr, marker="o", ms=5, color=PALETTE["orange"], lw=1.8, label="scramble", zorder=3)
    axd.set_xticks(xs)
    axd.set_xlabel("peptide position", fontsize=8.5)
    axd.set_ylabel("pLDDT", fontsize=8.5)
    axd.set_ylim(50, 100)
    axd.legend(loc="center left", bbox_to_anchor=(0.02, 0.42), frameon=False, fontsize=8.5)
    panel_label(axd, "d", x=-0.30)

    # ---- (e) anchor permissiveness, grouped bars -------------------------------
    hla_labels = {"HLA-A*02:01": "A*02:01", "HLA-A*11:01": "A*11:01"}
    xs_e = np.arange(len(anchor))
    w = 0.32
    cog_vals = [float(r["cognate_anchor_frac"]) for r in anchor]
    scr_vals = [float(r["scramble_anchor_frac"]) for r in anchor]
    axe.bar(xs_e - w / 2, cog_vals, width=w, color=PALETTE["blue"], edgecolor=SURFACE,
            linewidth=1.2, zorder=2, label="cognate")
    axe.bar(xs_e + w / 2, scr_vals, width=w, color=PALETTE["orange"], edgecolor=SURFACE,
            linewidth=1.2, zorder=2, label="scramble")
    for x, v in zip(xs_e - w / 2, cog_vals):
        axe.annotate(f"{v:.2f}", (x, v + 0.02), ha="center", fontsize=8, color=PALETTE["ink"])
    for x, v in zip(xs_e + w / 2, scr_vals):
        axe.annotate(f"{v:.2f}", (x, v + 0.02), ha="center", fontsize=8, color=PALETTE["ink"])
    axe.set_xticks(xs_e)
    axe.set_xticklabels([hla_labels[r["hla"]] for r in anchor], fontsize=9)
    axe.set_ylabel("anchor-residue permissiveness", fontsize=8.5)
    axe.set_ylim(0, 1.18)
    # legend sits in the empty central strip between the two HLA groups (no bar at x=0.5)
    axe.legend(loc="center", bbox_to_anchor=(0.5, 0.62), ncol=1,
               frameon=False, fontsize=8.5, handlelength=1.1)
    panel_label(axe, "e", x=-0.30)

    save(fig, "fig5_mhc_presentation")


if __name__ == "__main__":
    main()
