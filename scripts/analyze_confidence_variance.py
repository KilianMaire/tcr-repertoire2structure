"""Is the structural-confidence readout peptide-agnostic?

Read-only analysis over the already-folded panels (no GPU, no refold). It asks
the mechanistic question behind the failed pre-registered confirmation: does an
AlphaFold/Protenix confidence readout (iptm_TCRpep) carry information about the
TCR (its docking onto the shared HLA) but NOT about which peptide is the
cognate?

Design. For each TCR the panel is {cognate + same-HLA decoys}. We take the
per-(TCR, epitope) median readout across the 5 Protenix samples (the value the
benchmark ranks on). Then:

  1. ICC (intraclass correlation) across TCRs: how much of the readout variance
     is a between-TCR property vs within-TCR (across which peptide is shown).
     High ICC => the readout is set by the TCR, roughly constant over peptides.

  2. Sequential variance decomposition: eta^2 for TCR identity first, then, on
     the within-TCR residual, eta^2 for cognate status and for peptide identity.
     Prediction: TCR large, cognate ~ 0.

  3. Within-panel cognate effect: mean(cognate - panel mean) in ipTM units with
     a bootstrap CI, and a within-TCR label-permutation p (reassign which panel
     member is the pseudo-cognate). This is the retrieval signal expressed as an
     effect size rather than a single binary Top-1.

Usage: python scripts/analyze_confidence_variance.py runs/panel1 runs/hla_a1101
"""
from __future__ import annotations
import glob
import json
import random
import statistics as st
import sys
from pathlib import Path

READOUTS = {
    # index order in the construct: A=0 TCRa, B=1 TCRb, C=2 MHC, D=3 b2m, E=4 pep
    "iptm_TCRpep_max": lambda c: max(c["chain_pair_iptm"][0][4], c["chain_pair_iptm"][1][4]),
    "iptm_beta_pep": lambda c: c["chain_pair_iptm"][1][4],
    # negative-control geometry: MHC-peptide groove (should be TCR-agnostic AND
    # high for every peptide). Included to contrast its variance structure.
    "iptm_groove_ctrl": lambda c: c["chain_pair_iptm"][2][4],
}


def median_per_epitope(cid, ent, folds_root, fn):
    out = {}
    for ep in ent["epitopes"]:
        if ep == "__scramble__":
            continue
        vals = []
        for jp in glob.glob(str(Path(folds_root) / f"{cid}__{ep}" / "**"
                                 / "*summary_confidence_sample_*.json"), recursive=True):
            try:
                vals.append(fn(json.loads(Path(jp).read_text())))
            except (ValueError, OSError, KeyError, IndexError, TypeError):
                continue
        if vals:
            out[ep] = st.median(vals)
    return out


def load_panel(run_dir, fn):
    """Return list of dicts: {tcr, cognate, values: {ep: median}} for TCRs whose
    cognate was scored."""
    manifest = json.loads((Path(run_dir) / "manifest.json").read_text())
    folds_root = Path(run_dir) / "folds"
    rows = []
    for cid, ent in manifest.items():
        vals = median_per_epitope(cid, ent, folds_root, fn)
        if ent["cognate"] in vals and len(vals) >= 2:
            rows.append({"tcr": cid, "cognate": ent["cognate"], "values": vals})
    return rows


def icc(rows):
    """One-way ICC(1): between-TCR variance / total. Uses each TCR's panel of
    per-epitope medians as its group."""
    groups = [list(r["values"].values()) for r in rows]
    all_vals = [v for g in groups for v in g]
    grand = st.mean(all_vals)
    n_total = len(all_vals)
    k = n_total / len(groups)  # avg panel size
    ss_between = sum(len(g) * (st.mean(g) - grand) ** 2 for g in groups)
    ss_within = sum((v - st.mean(g)) ** 2 for g in groups for v in g)
    df_between = len(groups) - 1
    df_within = n_total - len(groups)
    ms_between = ss_between / df_between
    ms_within = ss_within / df_within if df_within else 0.0
    denom = ms_between + (k - 1) * ms_within
    return (ms_between - ms_within) / denom if denom else 0.0


def variance_decomposition(rows):
    """Sequential eta^2: TCR first, then cognate-status and peptide-identity on
    the within-TCR residual. Returns fractions of TOTAL sum of squares."""
    cells = []  # (tcr, ep, is_cognate, value)
    for r in rows:
        rowmean = st.mean(r["values"].values())
        for ep, v in r["values"].items():
            cells.append((r["tcr"], ep, ep == r["cognate"], v, rowmean))
    grand = st.mean([c[3] for c in cells])
    ss_total = sum((c[3] - grand) ** 2 for c in cells)
    if ss_total == 0:
        return {"tcr": 0.0, "cognate_within_tcr": 0.0, "peptide_within_tcr": 0.0}
    # TCR: row means vs grand
    ss_tcr = sum((c[4] - grand) ** 2 for c in cells)
    # residual after removing TCR mean
    resid = [(c[0], c[1], c[2], c[3] - c[4]) for c in cells]
    # cognate-status on residual: two group means (cognate vs not)
    def ss_factor(labels):
        groups = {}
        for lab, val in labels:
            groups.setdefault(lab, []).append(val)
        gm = st.mean([v for _, v in labels])
        return sum(len(g) * (st.mean(g) - gm) ** 2 for g in groups.values())
    ss_cog = ss_factor([(c[2], c[3]) for c in resid])
    ss_pep = ss_factor([(c[1], c[3]) for c in resid])
    return {"tcr": ss_tcr / ss_total,
            "cognate_within_tcr": ss_cog / ss_total,
            "peptide_within_tcr": ss_pep / ss_total}


def cognate_effect(rows, n_boot=5000, n_perm=10000, seed=0):
    """Within-panel: cognate value minus panel mean, in readout units. Bootstrap
    CI over TCRs, and a within-TCR permutation p (pseudo-cognate = random panel
    member) for 'cognate scores above a random panel member'."""
    deltas, panels = [], []
    for r in rows:
        vals = r["values"]
        pm = st.mean(vals.values())
        deltas.append(vals[r["cognate"]] - pm)
        panels.append((list(vals.values()), list(vals).index(r["cognate"])))
    obs = st.mean(deltas)
    rng = random.Random(seed)
    boots = sorted(st.mean([deltas[rng.randrange(len(deltas))] for _ in deltas])
                   for _ in range(n_boot))
    lo, hi = boots[int(0.025 * n_boot)], boots[int(0.975 * n_boot)]
    ge = 0
    for _ in range(n_perm):
        s = 0.0
        for vlist, _ in panels:
            pm = st.mean(vlist)
            s += vlist[rng.randrange(len(vlist))] - pm
        if s / len(panels) >= obs:
            ge += 1
    return {"mean_delta": obs, "ci": [lo, hi], "frac_pos": sum(d > 0 for d in deltas) / len(deltas),
            "perm_p": (ge + 1) / (n_perm + 1), "n": len(deltas)}


def report(run_dir, readout_name, fn):
    rows = load_panel(run_dir, fn)
    ic = icc(rows)
    vd = variance_decomposition(rows)
    ce = cognate_effect(rows)
    print(f"[{Path(run_dir).name}] {readout_name}  (n_TCR={len(rows)})")
    print(f"  ICC across TCRs           = {ic:.3f}   (1.0 = readout is a pure TCR property)")
    print(f"  variance explained:  TCR  = {vd['tcr']*100:5.1f}%")
    print(f"                    cognate = {vd['cognate_within_tcr']*100:5.1f}%   <- peptide-specificity signal")
    print(f"                    peptide = {vd['peptide_within_tcr']*100:5.1f}%")
    print(f"  cognate - panel mean      = {ce['mean_delta']:+.4f}  CI[{ce['ci'][0]:+.4f},{ce['ci'][1]:+.4f}]"
          f"  (won {ce['frac_pos']*100:.0f}%, perm p={ce['perm_p']:.4f})")
    print()


if __name__ == "__main__":
    dirs = sys.argv[1:] or ["runs/panel1", "runs/hla_a1101"]
    for d in dirs:
        for name, fn in READOUTS.items():
            report(d, name, fn)
