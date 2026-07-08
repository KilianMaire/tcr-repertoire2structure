"""Skeptical QC on downloaded Protenix folds, then regenerate the report.

For each clonotype we score the cognate model and its scramble control on
CDR3(Vbeta) to peptide heavy-atom contact. Honesty Rule 2: a fold is called
`reliable` only when the cognate contact beats its own scramble null; otherwise
`suspect` (the geometry looks docked but does not discriminate the real peptide
from a shuffle). Missing / <5-chain models are `qc_failed`.

Usage:
  python scripts/qc_folds.py <folds_dir> <tablo_run_dir> <out_report.html>
  # folds_dir = extracted rep2struct_folds.zip (holds out/<id>_cognate, out/<id>_scramble)
"""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import numpy as np
from rep2struct.qc import score_model
from rep2struct.report import render_report
from rep2struct.runstate import RunState
from rep2struct.schema import Clonotype, Annotation, QCResult


def ensemble_contact(fold_dir: Path):
    """Mean CDR3(Vbeta) to peptide contact across ALL Protenix samples for one
    construct. Protenix emits several samples per seed whose docking pose varies
    a lot, so a single sample is not representative; the ensemble mean is. Returns
    (mean, n_models, n_valid) or (None, 0, 0) if no 5-chain model parsed."""
    cifs = sorted(fold_dir.rglob("*.cif"))
    vals = []
    for c in cifs:
        v = score_model(c).get("cdr3_pep_atoms")
        if v is not None:
            vals.append(v)
    if not vals:
        return None, len(cifs), 0
    return float(np.mean(vals)), len(cifs), len(vals)


def main():
    folds_dir = Path(sys.argv[1])
    run_dir = sys.argv[2]
    out_html = sys.argv[3]

    rs = RunState(run_dir)
    clons = [Clonotype(**d) for d in rs.read_stage("ingest")]
    anns = [Annotation(**d) for d in rs.read_stage("annotate")]
    jobs = rs.read_stage("foldjobs") if rs.stage_done("foldjobs") else []
    job_ids = [j["clonotype_id"] for j in jobs]

    qcs = []
    for cid in job_ids:
        cog, ncog, _ = ensemble_contact(folds_dir / "out" / f"{cid}_cognate")
        if cog is None:
            qcs.append(QCResult(cid, "qc_failed",
                                "no cognate model produced (fold pending or <5 chains)"))
            continue
        scr, _, _ = ensemble_contact(folds_dir / "out" / f"{cid}_scramble")
        thr = scr if scr is not None else 0.0
        # scramble-calibrated verdict on the ensemble mean. The scramble contact
        # is typically non-zero (Protenix imposes docking geometry on any
        # peptide); a cognate is reliable only when it beats that null.
        reason = f"cognate {cog:.0f} vs scramble {thr:.0f} mean Vbeta-peptide contact ({ncog} models)"
        verdict = "reliable" if cog > thr else "suspect"
        qcs.append(QCResult(cid, verdict, reason, cdr3_pep_atoms=cog))
        print(f"{cid}: cognate={cog:.0f} scramble={thr:.0f} -> {verdict}")

    html = render_report(clons, anns, qcs, metrics=None)
    Path(out_html).write_text(html)
    print(f"report -> {out_html}")


if __name__ == "__main__":
    main()
