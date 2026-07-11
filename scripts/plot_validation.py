"""Fig 2: honest annotation and the leakage guard.

Reads the committed CSVs in paper/data/. Three panels:
  (a) precision, recall, unannotatable rate, raw vs de-leaked.
  (b) precision and recall across the TCRdist cut, raw vs de-leaked.
  (c) the leakage signature: distance percentiles of correct calls (median 0).

Usage: python scripts/plot_validation.py     # writes paper/figures/fig2_validation.png
"""
from __future__ import annotations
import csv, json, sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "paper/data"
C_RAW, C_DELEAK = "#0072B2", "#D55E00"   # validated CVD-safe pair
SURFACE = "#FCFCFB"


def _rows(name):
    with (DATA / name).open() as f:
        return list(csv.DictReader(f))


def main(out="paper/figures/fig2_validation.png"):
    ann = {r["set"]: r for r in _rows("validation_annotation.csv")}
    sweep = _rows("validation_threshold_sweep.csv")
    perc = json.loads((ROOT / "docs/validation_donor1_metrics.json").read_text())["correct_distance_percentiles"]

    fig, (axa, axb, axc) = plt.subplots(1, 3, figsize=(13, 4.3))
    fig.patch.set_facecolor(SURFACE)
    for ax in (axa, axb, axc):
        ax.set_facecolor(SURFACE)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)

    # (a) grouped bars
    cats = [("precision", "precision"), ("recall", "recall"), ("unannotatable_rate", "unannotatable")]
    x = range(len(cats))
    w = 0.38
    for i, (key, _) in enumerate(cats):
        axa.bar(i - w / 2, float(ann["raw"][key]), w, color=C_RAW, edgecolor=SURFACE,
                linewidth=1.5, label="raw" if i == 0 else None)
        axa.bar(i + w / 2, float(ann["deleaked"][key]), w, color=C_DELEAK, edgecolor=SURFACE,
                linewidth=1.5, label="de-leaked" if i == 0 else None)
        for off, s in ((-w / 2, "raw"), (w / 2, "deleaked")):
            v = float(ann[s][key])
            axa.annotate(f"{v:.2f}", (i + off, v + 0.02), ha="center", fontsize=8, color="#333")
    axa.set_xticks(list(x))
    axa.set_xticklabels([c[1] for c in cats], fontsize=9)
    axa.set_ylim(0, 1.05)
    axa.set_ylabel("rate")
    axa.set_title("(a) annotation, raw vs de-leaked", fontsize=10.5, loc="left")
    axa.legend(fontsize=8.5, frameon=False, loc="upper center")

    # (b) precision / recall across the cut
    cut = [float(r["tcrdist_cut"]) for r in sweep]
    axb.plot(cut, [float(r["precision_deleaked"]) for r in sweep], "-o", color=C_DELEAK,
             ms=4, lw=2, label="precision (de-leaked)")
    axb.plot(cut, [float(r["precision_raw"]) for r in sweep], "--o", color=C_RAW,
             ms=4, lw=1.6, label="precision (raw)")
    axb.plot(cut, [float(r["recall_deleaked"]) for r in sweep], "-s", color=C_DELEAK,
             ms=4, lw=1.4, alpha=0.55, label="recall (de-leaked)")
    axb.set_xlabel("TCRdist cut")
    axb.set_ylim(0, 1.05)
    axb.set_title("(b) precision holds, recall stays low", fontsize=10.5, loc="left")
    axb.legend(fontsize=8, frameon=False, loc="center right")
    for thr, name in ((12, "high"), (24, "med"), (48, "low")):
        axb.axvline(thr, color="#cccccc", lw=1, zorder=0)

    # (c) leakage signature: distance percentiles of correct calls
    ps = sorted(perc, key=lambda k: int(k))
    vals = [perc[p] for p in ps]
    axc.bar([f"P{p}" for p in ps], vals, color=C_RAW, edgecolor=SURFACE, linewidth=1.5, width=0.7)
    for i, v in enumerate(vals):
        axc.annotate(f"{v:g}", (i, v + 0.08), ha="center", fontsize=8.5, color="#333")
    axc.set_ylabel("TCRdist of correct call")
    axc.set_title("(c) leakage signature: correct calls sit at distance 0", fontsize=10.5, loc="left")
    axc.set_ylim(0, max(vals) + 0.8)

    fig.suptitle("Honest annotation: precise where it fires, abstains on novel TCRs, leakage-corrected",
                 fontsize=12.5, fontweight="bold", x=0.02, ha="left", y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    print(f"wrote {out}")


if __name__ == "__main__":
    main(*(sys.argv[1:2] or []))
