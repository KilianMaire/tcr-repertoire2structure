"""Data-quality QC: how poly-G TCR stubs contaminate the TCR-recognition analysis.

Surfaced by an audit. When a clonotype's V domain cannot be reconstructed,
build_tcr_seqs falls back to a poly-G stub (ten glycines + CDR3) for BOTH TCR
chains, and those constructs were still folded. A stubbed TCR is a glycine
backbone with a floating CDR3, so any TCR-peptide interface readout for it is
meaningless. This script quantifies the contamination and recomputes the
TCR-recognition numbers on genuine (reconstructed) TCRs only.

A clonotype is a stub iff its manifest chain_b_seq starts with ten glycines. This
is an objective, outcome-independent filter (a data-quality criterion, not a
result-driven one).

Note on the pre-registration: the frozen held-out result was computed on the FULL
panel, so the full-panel number remains the official pre-registered outcome. The
reconstructed-only recomputation here is a post-hoc data-quality reanalysis and is
reported as such; it cannot retroactively license a positive claim.

Usage: python scripts/analyze_stub_contamination.py runs/panel1 runs/hla_a1101
"""
from __future__ import annotations
import glob, json, statistics as st, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze_confidence_variance import variance_decomposition, cognate_effect

STUB_PREFIX = "G" * 10
_tcrpep = lambda c: max(c["chain_pair_iptm"][0][4], c["chain_pair_iptm"][1][4])
_groove = lambda c: c["chain_pair_iptm"][2][4]


def is_stub(chain_b_seq: str) -> bool:
    return chain_b_seq.startswith(STUB_PREFIX)


def _median_per_ep(cid, ent, folds_root, fn):
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


def _median(cid, ep, folds_root, fn):
    vals = []
    for jp in glob.glob(str(Path(folds_root) / f"{cid}__{ep}" / "**"
                             / "*summary_confidence_sample_*.json"), recursive=True):
        try:
            vals.append(fn(json.loads(Path(jp).read_text())))
        except (ValueError, OSError, KeyError, IndexError, TypeError):
            continue
    return st.median(vals) if vals else None


def _auroc(pos, neg):
    n = d = 0.0
    for p in pos:
        for q in neg:
            d += 1
            n += 1.0 if p > q else (0.5 if p == q else 0.0)
    return n / d if d else None


def analyze(run_dir):
    manifest = json.loads((Path(run_dir) / "manifest.json").read_text())
    folds = Path(run_dir) / "folds"
    n_total = len(manifest)
    n_stub = sum(is_stub(e["chain_b_seq"]) for e in manifest.values())
    out = {"run": Path(run_dir).name, "n_total": n_total, "n_stub": n_stub,
           "stub_frac": round(n_stub / n_total, 3)}
    for mode in ("all", "reconstructed"):
        rows, panels, binders, scrs = [], [], [], []
        for cid, ent in manifest.items():
            if mode == "reconstructed" and is_stub(ent["chain_b_seq"]):
                continue
            pe = _median_per_ep(cid, ent, folds, _tcrpep)
            if ent["cognate"] in pe:
                rows.append({"tcr": cid, "values": pe, "cognate": ent["cognate"]})
                panels.append((pe, ent["cognate"]))
            for ep in [ent["cognate"], *ent["decoys"]]:
                v = _median(cid, ep, folds, _groove)
                if v is not None:
                    binders.append(v)
            s = _median(cid, "__scramble__", folds, _groove)
            if s is not None:
                scrs.append(s)
        hits = sum(1 for p, cog in panels
                   if p[cog] == max(p.values()) and sum(v == max(p.values()) for v in p.values()) == 1)
        vd = variance_decomposition(rows)
        ce = cognate_effect(rows)
        out[mode] = {
            "n": len(rows), "top1": round(hits / len(rows), 3),
            "var_tcr": round(vd["tcr"], 3), "var_cognate": round(vd["cognate_within_tcr"], 3),
            "cognate_delta": round(ce["mean_delta"], 4),
            "cognate_p": round(ce["perm_p"], 4),
            "groove_auroc": round(_auroc(binders, scrs), 3),
        }
    return out


def main(dirs):
    rows_csv = [("run", "mode", "n", "tcrpep_top1", "var_tcr", "var_cognate",
                 "cognate_delta", "cognate_perm_p", "groove_auroc")]
    for run_dir in dirs:
        r = analyze(run_dir)
        print(f"\n=== {r['run']}: {r['n_stub']}/{r['n_total']} stubbed "
              f"({r['stub_frac']:.0%}) ===")
        print(f"{'mode':14} {'n':>3} {'TCRpep_top1':>12} {'var_TCR':>8} {'var_cog':>8} "
              f"{'cog_delta':>10} {'cog_p':>7} {'groove_AUROC':>13}")
        for mode in ("all", "reconstructed"):
            m = r[mode]
            print(f"{mode:14} {m['n']:>3} {m['top1']:>12.3f} {m['var_tcr']:>8.3f} "
                  f"{m['var_cognate']:>8.3f} {m['cognate_delta']:>+10.4f} {m['cognate_p']:>7.4f} "
                  f"{m['groove_auroc']:>13.3f}")
            rows_csv.append((r["run"], mode, m["n"], m["top1"], m["var_tcr"], m["var_cognate"],
                             m["cognate_delta"], m["cognate_p"], m["groove_auroc"]))
    return rows_csv


if __name__ == "__main__":
    main(sys.argv[1:] or ["runs/panel1", "runs/hla_a1101"])
