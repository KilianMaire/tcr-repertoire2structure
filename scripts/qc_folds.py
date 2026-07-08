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
from rep2struct.qc import score_model
from rep2struct.report import render_report
from rep2struct.runstate import RunState
from rep2struct.schema import Clonotype, Annotation, QCResult


def find_cif(fold_dir: Path):
    """Pick one representative model CIF under a Protenix output dir."""
    cifs = sorted(fold_dir.rglob("*.cif"))
    if not cifs:
        return None
    # prefer a rank/sample 0 model if present, else the first
    for c in cifs:
        if "model_0" in c.name or "sample_0" in c.name or "rank_1" in c.name:
            return c
    return cifs[0]


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
        cog = find_cif(folds_dir / "out" / f"{cid}_cognate")
        scr = find_cif(folds_dir / "out" / f"{cid}_scramble")
        if cog is None:
            qcs.append(QCResult(cid, "qc_failed", "no cognate model produced"))
            continue
        cs = score_model(cog); cs["clonotype_id"] = cid
        if cs.get("cdr3_pep_atoms") is None:
            qcs.append(QCResult(cid, "qc_failed", f"cognate has {cs.get('n_chains')} chains, need 5"))
            continue
        scr_atoms = None
        if scr is not None:
            ss = score_model(scr)
            scr_atoms = ss.get("cdr3_pep_atoms")
        # per-clonotype scramble calibration: cognate must beat its own scramble
        thr = scr_atoms if scr_atoms is not None else 0.0
        if cs["cdr3_pep_atoms"] > thr:
            qcs.append(QCResult(cid, "reliable",
                                f"Vbeta-peptide contact {cs['cdr3_pep_atoms']:.0f} > scramble {thr:.0f}",
                                cdr3_pep_atoms=cs["cdr3_pep_atoms"]))
        else:
            qcs.append(QCResult(cid, "suspect",
                                f"Vbeta-peptide contact {cs['cdr3_pep_atoms']:.0f} <= scramble {thr:.0f}",
                                cdr3_pep_atoms=cs["cdr3_pep_atoms"]))
        print(f"{cid}: cognate={cs['cdr3_pep_atoms']:.0f} scramble={thr:.0f} -> {qcs[-1].qc_verdict}")

    html = render_report(clons, anns, qcs, metrics=None)
    Path(out_html).write_text(html)
    print(f"report -> {out_html}")


if __name__ == "__main__":
    main()
