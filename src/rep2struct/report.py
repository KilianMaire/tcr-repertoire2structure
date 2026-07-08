from __future__ import annotations
from pathlib import Path
from jinja2 import Environment, FileSystemLoader, select_autoescape

_TPL_DIR = Path(__file__).parent / "templates"

def render_report(clonotypes, annotations, qc_results, metrics=None) -> str:
    ann = {a.clonotype_id: a for a in annotations}
    qc = {q.clonotype_id: q for q in qc_results}
    rows = []
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
        })
    env = Environment(loader=FileSystemLoader(str(_TPL_DIR)),
                      autoescape=select_autoescape(["html"]))
    return env.get_template("report.html.j2").render(rows=rows, metrics=metrics)
