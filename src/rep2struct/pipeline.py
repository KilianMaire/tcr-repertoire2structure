from __future__ import annotations
from pathlib import Path
from .runstate import RunState
from .ingest import parse_10x, standardize_alleles
from .annotate import annotate
from .foldprep import select_top, build_construct
from .fold import run_folds
from .qc import score_model, verdict
from .report import render_report
from .validate import annotation_metrics
from .schema import QCResult


def _tcr_seq_stub(clonotype):
    # placeholder chain sequences: real runs use reconstructed V domains from
    # TCR Explorer. For the construct we need any residues; use the CDR3s padded.
    return {"A": "G" * 10 + clonotype.cdr3a, "B": "G" * 10 + clonotype.cdr3b}


def run_pipeline(csv_path, run_dir, top_n, sim_fn=None, assign_fn=None, fold_fn=None,
                 tcr_seqs=None, mhc_seqs=None, scramble_threshold=1.0, labels=None):
    rs = RunState(run_dir)

    clons = standardize_alleles(parse_10x(csv_path), assign_fn=assign_fn)
    anns = annotate(clons, sim_fn=sim_fn)

    metrics = annotation_metrics(anns, labels) if labels else None

    top = select_top(clons, anns, n=top_n)
    # Guard: select_top can return clonotypes whose annotation is not
    # foldable (unannotatable, hla is None). build_construct would raise
    # KeyError on mhc_seqs[None] for those, so filter to foldable entries
    # only. Filtered-out clonotypes simply get no fold job and show up in
    # the report as "not folded".
    foldable = [(c, a) for c, a in top if a.annotatable and a.hla]

    seqs = tcr_seqs or {c.id: _tcr_seq_stub(c) for c, _ in foldable}
    jobs = [build_construct(c, a, seqs, mhc_seqs) for c, a in foldable]
    jobs = run_folds(jobs, fold_fn, rs)

    qcs = []
    for job in jobs:
        if job.status != "done" or not job.model_paths:
            qcs.append(QCResult(job.clonotype_id, "qc_failed", "no model produced"))
            continue
        s = score_model(job.model_paths[0]); s["clonotype_id"] = job.clonotype_id
        qcs.append(verdict(s, scramble_threshold))

    html = render_report(clons, anns, qcs, metrics=metrics)
    out = Path(run_dir) / "report.html"
    out.write_text(html)
    return str(out)
