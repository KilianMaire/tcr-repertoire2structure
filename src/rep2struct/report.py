from __future__ import annotations
from pathlib import Path
from jinja2 import Environment, FileSystemLoader, select_autoescape

_TPL_DIR = Path(__file__).parent / "templates"

# Evidence-type label per QCResult.qc_verdict. presented/not_presented rows are
# binding-score predictions, never a fold or structure claim.
EVIDENCE = {
    "reliable": "structure (reliable)",
    "suspect": "structure (suspect)",
    "qc_failed": "structure (qc failed)",
    "presented": "predicted presentation",
    "not_presented": "predicted presentation",
}

def _msa_note(basis) -> str:
    if basis in ("local", "colab_cpu"):
        return f"MSA {basis}"
    return "MSA-free (reduced confidence)"


def render_report(clonotypes, annotations, qc_results, metrics=None, msa_basis=None) -> str:
    ann = {a.clonotype_id: a for a in annotations}
    qc = {q.clonotype_id: q for q in qc_results}
    msa_basis = msa_basis or {}
    rows = []
    # Keep clonotype input order; never sort by a raw numeric field
    # (e.g. cdr3_pep_atoms) so distances are never presented as a ranking.
    for c in clonotypes:
        a = ann.get(c.id)
        q = qc.get(c.id)
        rows.append({
            "id": c.id, "size": c.size,
            "annotatable": bool(a and a.annotatable),
            "epitope": a.epitope if a else None,
            "hla": a.hla if a else None,
            "tier": a.confidence_tier if a else "n/a",
            "qc": q.qc_verdict if q else "not folded",
            "tool": q.tool if q else None,
            "evidence": EVIDENCE.get(q.qc_verdict, "structure") if q else "n/a",
            "msa_note": _msa_note(msa_basis.get(c.id)) if q else None,
        })
    env = Environment(loader=FileSystemLoader(str(_TPL_DIR)),
                      autoescape=select_autoescape(["html"]))
    return env.get_template("report.html.j2").render(rows=rows, metrics=metrics)
