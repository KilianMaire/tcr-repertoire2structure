"""Which structural metric cleanly judges MHC-peptide specificity (presentation)?

Read-only over the folded panels (no refold). The scramble is a composition
preserving shuffle of the cognate, so it holds amino-acid content and length
fixed and only breaks the anchor order: a clean negative for peptide binding.
The cognate AND the same-HLA decoys are all genuine binders. So a metric that
scores presentation should separate binders (cognate + decoys) from the scramble
regardless of TCR recognition, which is the mirror image of the TCR result where
no metric separates cognate from decoy.

For each candidate metric we report, per run:
  frac cog>scr   fraction of TCRs whose cognate median beats its own scramble
  frac dec>scr   same for pooled decoys (confirms it is a binding signal)
  delta          mean(cognate - scramble) with a bootstrap CI over TCRs
  AUROC          binder medians vs scramble medians, pooled, with a TCR-bootstrap CI

Usage: python scripts/analyze_mhc_scramble.py runs/panel1 runs/hla_a1101
"""
from __future__ import annotations
import glob, json, random, statistics as st, sys
from pathlib import Path

# chain index in the construct: A=0 TCRa, B=1 TCRb, C=2 MHC, D=3 b2m, E=4 pep
METRICS = {
    "iptm_groove": lambda c: c["chain_pair_iptm"][2][4],       # MHC-peptide interface ipTM
    "neg_gpde_groove": lambda c: -c["chain_pair_gpde"][2][4],  # MHC-peptide PAE-analog (neg so higher=better)
    "iptm_b2m_pep": lambda c: c["chain_pair_iptm"][3][4],
    "pep_plddt": lambda c: c["chain_plddt"][4],
    "pep_ptm": lambda c: c["chain_ptm"][4],
    "pep_iptm": lambda c: c["chain_iptm"][4],
    "ranking_score": lambda c: c["ranking_score"],
}


def _median(cid, ep, folds_root, fn):
    vals = []
    for jp in glob.glob(str(Path(folds_root) / f"{cid}__{ep}" / "**"
                             / "*summary_confidence_sample_*.json"), recursive=True):
        try:
            vals.append(fn(json.loads(Path(jp).read_text())))
        except (ValueError, OSError, KeyError, IndexError, TypeError):
            continue
    return st.median(vals) if vals else None


def collect(run_dir, fn):
    """Per TCR: (cognate_median, scramble_median, [decoy_medians])."""
    manifest = json.loads((Path(run_dir) / "manifest.json").read_text())
    folds_root = Path(run_dir) / "folds"
    rows = []
    for cid, ent in manifest.items():
        cog = _median(cid, ent["cognate"], folds_root, fn)
        scr = _median(cid, "__scramble__", folds_root, fn)
        if cog is None or scr is None:
            continue
        dec = [d for d in (_median(cid, ep, folds_root, fn) for ep in ent["decoys"]) if d is not None]
        rows.append((cog, scr, dec))
    return rows


def _auroc(pos, neg):
    n = d = 0.0
    for p in pos:
        for q in neg:
            d += 1
            n += 1.0 if p > q else (0.5 if p == q else 0.0)
    return n / d if d else None


def auroc_binder_vs_scramble(rows):
    binders = [v for (cog, _, dec) in rows for v in [cog, *dec]]
    scram = [scr for (_, scr, _) in rows]
    return _auroc(binders, scram)


def _boot_ci(sample_fn, seed=0, nb=2000):
    rng = random.Random(seed)
    vals = sorted(v for v in (sample_fn(rng) for _ in range(nb)) if v is not None)
    if not vals:
        return None, None
    return vals[int(0.025 * len(vals))], vals[int(0.975 * len(vals))]


def summarize(run_dir, fn):
    rows = collect(run_dir, fn)
    deltas = [cog - scr for (cog, scr, _) in rows]
    dec_deltas = [d - scr for (_, scr, dec) in rows for d in dec]
    au = auroc_binder_vs_scramble(rows)

    def resample_delta(rng):
        return st.mean(deltas[rng.randrange(len(deltas))] for _ in deltas)

    def resample_auroc(rng):
        boot = [rows[rng.randrange(len(rows))] for _ in rows]
        return auroc_binder_vs_scramble(boot)

    dlo, dhi = _boot_ci(resample_delta)
    alo, ahi = _boot_ci(resample_auroc)
    return {
        "n": len(rows),
        "frac_cog_gt_scr": sum(d > 0 for d in deltas) / len(deltas) if deltas else None,
        "frac_dec_gt_scr": sum(d > 0 for d in dec_deltas) / len(dec_deltas) if dec_deltas else None,
        "delta": st.mean(deltas) if deltas else None, "delta_ci": [dlo, dhi],
        "auroc": au, "auroc_ci": [alo, ahi],
    }


def main(dirs):
    for run_dir in dirs:
        print(f"\n=== {Path(run_dir).name} ===")
        print(f"{'metric':16} {'cog>scr':>8} {'dec>scr':>8} {'delta':>9} "
              f"{'delta_CI':>18} {'AUROC':>7} {'AUROC_CI':>16}")
        scored = {m: summarize(run_dir, fn) for m, fn in METRICS.items()}
        for m in sorted(scored, key=lambda m: -(scored[m]["auroc"] or 0)):
            s = scored[m]
            print(f"{m:16} {s['frac_cog_gt_scr']:>8.2f} {s['frac_dec_gt_scr']:>8.2f} "
                  f"{s['delta']:>+9.3f} [{s['delta_ci'][0]:>+7.3f},{s['delta_ci'][1]:>+7.3f}] "
                  f"{s['auroc']:>7.3f} [{s['auroc_ci'][0]:>6.3f},{s['auroc_ci'][1]:>6.3f}]")


if __name__ == "__main__":
    main(sys.argv[1:] or ["runs/panel1", "runs/hla_a1101"])
