from __future__ import annotations
import warnings
import numpy as np
from .schema import QCResult

def _heavy_by_chain(cif_path):
    from Bio.PDB import MMCIFParser
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = next(MMCIFParser(QUIET=True).get_structure("x", str(cif_path)).get_models())
    out = {}
    for ch in m:
        atoms = [a.coord for r in ch for a in r if a.element != "H"]
        if atoms:
            out[ch.id] = np.array(atoms)
    return out


def load_chains(cif_path):
    return _heavy_by_chain(cif_path)


def common_checks(chains: dict, expected: set) -> dict:
    issues = []
    missing = expected - set(chains)
    if missing:
        issues.append(f"missing chains {sorted(missing)}")
    finite = all(np.isfinite(a).all() for a in chains.values())
    if not finite:
        issues.append("non-finite coords")
    has_peptide = "E" in chains and len(chains["E"]) > 0
    min_inter = None
    if finite and len(chains) >= 2:
        ids = list(chains)
        best = np.inf
        for i in range(len(ids)):
            for j in range(i + 1, len(ids)):
                a, b = chains[ids[i]], chains[ids[j]]
                d = np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(-1))
                best = min(best, float(d.min()))
        min_inter = best
        if min_inter < 0.5:
            issues.append("severe steric clash")
    return {"ok": not issues, "issues": issues, "n_chains": len(chains),
            "has_peptide": has_peptide, "min_interatomic": min_inter}

def score_model(cif_path) -> dict:
    chains = _heavy_by_chain(cif_path)
    if not {"A", "B", "C", "D", "E"}.issubset(chains):
        return {"n_chains": len(chains), "cdr3_pep_atoms": None, "crossing_angle": None, "dockq": None}
    pep = chains["E"]
    beta = chains["B"]
    d = np.sqrt(((beta[:, None, :] - pep[None, :, :]) ** 2).sum(-1))
    cdr3_pep_atoms = float((d < 4.5).sum())
    return {"n_chains": len(chains), "cdr3_pep_atoms": cdr3_pep_atoms,
            "crossing_angle": None, "dockq": None}

def verdict(scores, scramble_threshold: float) -> QCResult:
    cid = scores.get("clonotype_id", "unknown")
    if scores.get("cdr3_pep_atoms") is None:
        return QCResult(cid, "qc_failed", f"model has {scores.get('n_chains')} chains, need 5")
    if scores["cdr3_pep_atoms"] <= scramble_threshold:
        return QCResult(cid, "suspect",
                        "CDR3 to peptide contact not above scramble calibration",
                        cdr3_pep_atoms=scores["cdr3_pep_atoms"])
    return QCResult(cid, "reliable", "CDR3 to peptide contact beats scramble null",
                    cdr3_pep_atoms=scores["cdr3_pep_atoms"])

def verdict_binding(score: float, threshold: float, clonotype_id: str, tool: str) -> QCResult:
    presented = score > threshold
    return QCResult(
        clonotype_id,
        "presented" if presented else "not_presented",
        ("predicted presentation above the score null" if presented
         else "predicted presentation not above the score null"),
        tool=tool,
        calibration_basis="binding_score_null",
    )
