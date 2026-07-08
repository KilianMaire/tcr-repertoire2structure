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
