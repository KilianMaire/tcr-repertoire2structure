"""Fig 2: honest annotation and the leakage guard.

Reads the committed CSVs in paper/data/ and docs/validation_donor1_metrics.json.
Panels (2x2):
  (a) precision, recall, unannotatable rate, raw vs de-leaked.
  (b) precision and recall across the TCRdist cut, raw vs de-leaked.
  (c) TCRdist percentiles of correct calls (P10..P90).
  (d) abstention waterfall: labeled -> de-leaked scored -> predicted -> correct.

No titles, no interpretation on the axes: the caption carries the reading.

Usage: python scripts/plot_validation.py
"""
from __future__ import annotations
import csv, json, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import matplotlib.pyplot as plt
from figstyle import PALETTE, SURFACE, ROOT, DATA, apply, despine, panel_label, save

C_RAW, C_DELEAK = PALETTE["blue"], PALETTE["orange"]


def _rows(name):
    with (DATA / name).open() as f:
        return list(csv.DictReader(f))


def main():
    apply()
    ann = {r["set"]: r for r in _rows("validation_annotation.csv")}
    sweep = _rows("validation_threshold_sweep.csv")
    perc = json.loads((ROOT / "docs/validation_donor1_metrics.json").read_text())["correct_distance_percentiles"]

    fig = plt.figure(figsize=(12, 8))
    gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.30,
                           left=0.07, right=0.98, top=0.95, bottom=0.08)
    axa = fig.add_subplot(gs[0, 0])
    axb = fig.add_subplot(gs[0, 1])
    axc = fig.add_subplot(gs[1, 0])
    axd = fig.add_subplot(gs[1, 1])
    for ax in (axa, axb, axc, axd):
        despine(ax)

    # (a) grouped bars: precision, recall, unannotatable rate
    cats = [("precision", "precision"), ("recall", "recall"), ("unannotatable_rate", "unannotatable")]
    w = 0.36
    for i, (key, _) in enumerate(cats):
        axa.bar(i - w / 2, float(ann["raw"][key]), w, color=C_RAW, edgecolor=SURFACE,
                linewidth=1.5, label="raw" if i == 0 else None, zorder=2)
        axa.bar(i + w / 2, float(ann["deleaked"][key]), w, color=C_DELEAK, edgecolor=SURFACE,
                linewidth=1.5, label="de-leaked" if i == 0 else None, zorder=2)
        for off, s in ((-w / 2, "raw"), (w / 2, "deleaked")):
            v = float(ann[s][key])
            axa.annotate(f"{v:.2f}", (i + off, v + 0.025), ha="center", fontsize=8.5,
                         color=PALETTE["ink"])
    axa.set_xticks(range(len(cats)))
    axa.set_xticklabels([c[1] for c in cats], fontsize=9)
    axa.set_ylim(0, 1.12)
    axa.set_ylabel("rate")
    axa.legend(fontsize=8.5, frameon=False, loc="upper center", ncol=2)
    panel_label(axa, "a")

    # (b) precision / recall across the TCRdist cut
    cut = [float(r["tcrdist_cut"]) for r in sweep]
    axb.plot(cut, [float(r["precision_deleaked"]) for r in sweep], "-o", color=C_DELEAK,
             ms=4, lw=2, label="precision (de-leaked)", zorder=3)
    axb.plot(cut, [float(r["precision_raw"]) for r in sweep], "--o", color=C_RAW,
             ms=4, lw=1.6, label="precision (raw)", zorder=3)
    axb.plot(cut, [float(r["recall_deleaked"]) for r in sweep], "-s", color=PALETTE["amber"],
             ms=4, lw=1.4, label="recall (de-leaked)", zorder=3)
    for thr in (12, 24, 48):
        axb.axvline(thr, color=PALETTE["gridline"], lw=1, zorder=0)
    axb.set_xlabel("TCRdist cut")
    axb.set_ylabel("rate")
    axb.set_ylim(0, 1.12)
    axb.set_xlim(0, 95)
    axb.legend(fontsize=8, frameon=False, loc="upper right")
    panel_label(axb, "b")

    # (c) TCRdist percentiles of correct calls
    ps = sorted(perc, key=lambda k: int(k))
    vals = [perc[p] for p in ps]
    axc.bar([f"P{p}" for p in ps], vals, color=C_RAW, edgecolor=SURFACE, linewidth=1.5, width=0.65, zorder=2)
    for i, v in enumerate(vals):
        axc.annotate(f"{v:g}", (i, v + 0.12), ha="center", fontsize=8.5, color=PALETTE["ink"])
    axc.set_xlabel("percentile of correct calls")
    axc.set_ylabel("TCRdist of correct call")
    axc.set_ylim(0, max(vals) + 1.2)
    panel_label(axc, "c")

    # (d) abstention waterfall: labeled -> de-leaked scored -> predicted -> correct
    steps = [
        ("labeled\nclonotypes", 3325, PALETTE["grey"]),
        ("de-leaked\nscored", 1853, PALETTE["blue"]),
        ("predicted", 193, PALETTE["amber"]),
        ("correct", 151, PALETTE["orange"]),
    ]
    xs = range(len(steps))
    heights = [s[1] for s in steps]
    colors = [s[2] for s in steps]
    labels = [s[0] for s in steps]
    axd.bar(xs, heights, color=colors, edgecolor=SURFACE, linewidth=1.5, width=0.6, zorder=2)
    for i, v in enumerate(heights):
        axd.annotate(f"{v:,}", (i, v + 60), ha="center", fontsize=9, color=PALETTE["ink"])
    axd.set_xticks(list(xs))
    axd.set_xticklabels(labels, fontsize=8.5)
    axd.set_ylabel("clonotypes")
    axd.set_ylim(0, 3325 * 1.14)
    panel_label(axd, "d")

    save(fig, "fig2_validation")


if __name__ == "__main__":
    main()
