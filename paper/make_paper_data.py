"""Regenerate the paper's derived data tables from the raw folded panels.

Reproducibility layer. The raw Protenix folds live in runs/ (gitignored, need a
GPU to regenerate). This script turns them into small tidy CSVs under paper/data/
that every figure and table reads. Commit the CSVs; the figures are then
reproducible without a GPU. Re-run this only when the raw folds change.

Usage: python paper/make_paper_data.py
"""
from __future__ import annotations
import csv, glob, json, statistics as st, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from analyze_confidence_variance import load_panel, variance_decomposition, cognate_effect, icc
from analyze_mhc_scramble import METRICS as MHC_METRICS, summarize
from check_scramble_anchors import MOTIF, _pep_chain_E, _median_groove

OUT = Path(__file__).resolve().parent / "data"
RUNS = [("panel1", "A*02:01", ROOT / "runs/panel1"),
        ("hla_a1101", "A*11:01", ROOT / "runs/hla_a1101")]

# full confidence battery, chain order 0=TCRa 1=TCRb 2=MHC 3=b2m 4=pep
BATTERY = {
    "iptm_TCRpep_max": lambda c: max(c["chain_pair_iptm"][0][4], c["chain_pair_iptm"][1][4]),
    "iptm_TCRpep_mean": lambda c: (c["chain_pair_iptm"][0][4] + c["chain_pair_iptm"][1][4]) / 2,
    "iptm_beta_pep": lambda c: c["chain_pair_iptm"][1][4],
    "iptm_alpha_pep": lambda c: c["chain_pair_iptm"][0][4],
    "neg_gpde_beta_pep": lambda c: -c["chain_pair_gpde"][1][4],
    "iptm_global": lambda c: c["iptm"],
    "ptm_global": lambda c: c["ptm"],
    "ranking_score": lambda c: c["ranking_score"],
    "iptm_groove_ctrl": lambda c: c["chain_pair_iptm"][2][4],
}


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


def _top1(cid_to_panel):
    hits = won = 0
    for panel, cog in cid_to_panel:
        if cog not in panel:
            continue
        won += 1
        mx = max(panel.values())
        if panel[cog] == mx and sum(1 for v in panel.values() if v == mx) == 1:
            hits += 1
    return hits, won


def _write(name, header, rows):
    with (OUT / name).open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)
    print(f"wrote paper/data/{name}  ({len(rows)} rows)")


def validation_tables():
    d = json.loads((ROOT / "docs/validation_donor1_metrics.json").read_text())
    rows = []
    for kind in ("metrics_raw", "metrics_deleaked"):
        m = d[kind]
        rows.append([kind.replace("metrics_", ""), m["precision"], m["recall"],
                     m["unannotatable_rate"], m["n"], m["n_annotated"], m["n_correct"]])
    _write("validation_annotation.csv",
           ["set", "precision", "recall", "unannotatable_rate", "n", "n_annotated", "n_correct"], rows)
    sweep = [[s["thr"], s["n_pred"], s["precision"], s["recall"],
              s["precision_deleaked"], s["recall_deleaked"]] for s in d["threshold_sweep"]]
    _write("validation_threshold_sweep.csv",
           ["tcrdist_cut", "n_pred", "precision_raw", "recall_raw",
            "precision_deleaked", "recall_deleaked"], sweep)
    _write("validation_leakage.csv", ["n_labeled", "n_leakage_suspected", "frac_leaked"],
           [[d["n_labeled_clonotypes"], d["n_leakage_suspected"],
             d["n_leakage_suspected"] / d["n_labeled_clonotypes"]]])


def tcr_tables():
    ret_rows, var_rows = [], []
    for run, hla, run_dir in RUNS:
        manifest = json.loads((run_dir / "manifest.json").read_text())
        folds = run_dir / "folds"
        # retrieval Top-1 for the full battery
        for readout, fn in BATTERY.items():
            panels = []
            for cid, ent in manifest.items():
                pe = _median_per_ep(cid, ent, folds, fn)
                if ent["cognate"] in pe:
                    panels.append((pe, ent["cognate"]))
            hits, won = _top1(panels)
            ret_rows.append([run, hla, readout, round(hits / won, 4) if won else "", won])
        # variance decomposition for the three interpretable readouts
        for readout in ("iptm_TCRpep_max", "iptm_beta_pep", "iptm_groove_ctrl"):
            rows = load_panel(run_dir, BATTERY[readout])
            vd = variance_decomposition(rows)
            ce = cognate_effect(rows)
            var_rows.append([run, hla, readout, round(icc(rows), 4),
                             round(vd["tcr"], 4), round(vd["cognate_within_tcr"], 4),
                             round(vd["peptide_within_tcr"], 4), round(ce["mean_delta"], 4),
                             round(ce["ci"][0], 4), round(ce["ci"][1], 4), round(ce["perm_p"], 4)])
    _write("tcr_retrieval_top1.csv", ["run", "hla", "readout", "top1", "n_tcr"], ret_rows)
    _write("confidence_variance.csv",
           ["run", "hla", "readout", "icc", "var_tcr", "var_cognate", "var_peptide",
            "cognate_delta", "delta_ci_lo", "delta_ci_hi", "perm_p"], var_rows)


def mhc_tables():
    rows = []
    for run, hla, run_dir in RUNS:
        for metric, fn in MHC_METRICS.items():
            s = summarize(run_dir, fn)
            rows.append([run, hla, metric, round(s["frac_cog_gt_scr"], 3),
                         round(s["frac_dec_gt_scr"], 3), round(s["delta"], 4),
                         round(s["delta_ci"][0], 4), round(s["delta_ci"][1], 4),
                         round(s["auroc"], 4), round(s["auroc_ci"][0], 4), round(s["auroc_ci"][1], 4)])
    _write("mhc_presentation.csv",
           ["run", "hla", "metric", "frac_cog_gt_scr", "frac_dec_gt_scr", "delta",
            "delta_ci_lo", "delta_ci_hi", "auroc", "auroc_ci_lo", "auroc_ci_hi"], rows)


def anchor_table():
    groove = lambda c: c["chain_pair_iptm"][2][4]
    rows = []
    for run, hla, run_dir in RUNS:
        manifest = json.loads((run_dir / "manifest.json").read_text())
        motif = MOTIF.get(hla if hla.startswith("HLA-") else f"HLA-{hla}", lambda p: None)
        cog_ok = scr_ok = n = 0
        binders, scrambles = [], []
        for cid, ent in manifest.items():
            scr_pep = _pep_chain_E(run_dir / "constructs" / f"{cid}____scramble__.fasta")
            if scr_pep is None:
                continue
            n += 1
            cog_ok += bool(motif(ent["cognate"]))
            scr_ok += bool(motif(scr_pep))
            for ep in [ent["cognate"], *ent["decoys"]]:
                v = _median_groove(cid, ep, run_dir / "folds")
                if v is not None:
                    binders.append(v)
            sv = _median_groove(cid, "__scramble__", run_dir / "folds")
            if sv is not None:
                scrambles.append(sv)
        rows.append([run, f"HLA-{hla}", n, round(cog_ok / n, 3), round(scr_ok / n, 3),
                     round(st.median(binders), 4), round(st.median(scrambles), 4)])
    _write("scramble_anchor_permissiveness.csv",
           ["run", "hla", "n_tcr", "cognate_anchor_frac", "scramble_anchor_frac",
            "binder_groove_median", "scramble_groove_median"], rows)


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    validation_tables()
    tcr_tables()
    mhc_tables()
    anchor_table()
    print("done")
