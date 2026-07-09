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
    "pose_only": "pose (peptide in groove)",
    "pose_failed": "pose (qc failed)",
}

def msa_basis_from_manifest(manifest: dict) -> str:
    """Per-clonotype MSA manifest {chain_id: {"got_msa": bool}, ...} -> basis token.
    "colab_cpu:k/n" when any chain got an MSA, else "none" (honestly MSA-free)."""
    chains = [v for v in manifest.values() if isinstance(v, dict) and "got_msa" in v]
    n = len(chains)
    k = sum(1 for v in chains if v["got_msa"])
    return f"colab_cpu:{k}/{n}" if k else "none"


def _msa_note(basis) -> str:
    if basis and basis.startswith("colab_cpu"):
        _, _, cnt = basis.partition(":")
        return f"MSA colab_cpu ({cnt} chains)" if cnt else "MSA colab_cpu"
    if basis == "local":
        return "MSA local"
    return "MSA-free (reduced confidence)"


def render_report(clonotypes, annotations, qc_results, metrics=None, msa_basis=None,
                   validity=None) -> str:
    ann = {a.clonotype_id: a for a in annotations}
    qc = {q.clonotype_id: q for q in qc_results}
    msa_basis = msa_basis or {}
    validity = validity or {}
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
            "validity": validity.get(c.id, "n/a"),
        })
    env = Environment(loader=FileSystemLoader(str(_TPL_DIR)),
                      autoescape=select_autoescape(["html"]))
    return env.get_template("report.html.j2").render(rows=rows, metrics=metrics)
