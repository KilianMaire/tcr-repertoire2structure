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
from .seqs import build_tcr_seqs, build_mhc_seqs
from .schema import QCResult


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

    # Real V domains (reconstructed from V/J germline + CDR3) and real MHC
    # ectodomains (fetched + cached from IPD/IMGT-HLA), unless the caller
    # injected sequences (offline tests do). build_mhc_seqs may omit an HLA it
    # cannot resolve, so fold only clonotypes whose HLA has a heavy chain.
    seqs = tcr_seqs or build_tcr_seqs([c for c, _ in foldable])
    mhc = mhc_seqs or build_mhc_seqs(sorted({a.hla for _, a in foldable}))
    # build_construct returns None for a class II allele (unrepresentable by the class I
    # construct); drop those so they get no fold job rather than a mis-modelled complex.
    jobs = [j for j in (build_construct(c, a, seqs, mhc)
                        for c, a in foldable if a.hla in mhc) if j is not None]
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
