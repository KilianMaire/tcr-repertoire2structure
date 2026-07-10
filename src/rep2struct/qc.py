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


_THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


def _cdr3b_beta_atoms(cif_path, cdr3b):
    """Heavy atoms of the chain-B residues that spell the CDR3beta loop, located by
    matching the cdr3b one-letter sequence against chain B's residues (same approach
    as model_cdr3b_plddt). Returns an (N,3) array, or None if cdr3b is not found in
    chain B (caller must treat None as unmapped, never silently fall back to the whole
    chain, or the metric would mix two definitions)."""
    from Bio.PDB import MMCIFParser
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = next(MMCIFParser(QUIET=True).get_structure("x", str(cif_path)).get_models())
    chb = next((c for c in m if c.id == "B"), None)
    if chb is None:
        return None
    res = list(chb)
    seq = "".join(_THREE_TO_ONE.get(r.resname, "X") for r in res)
    i = seq.find(cdr3b)
    if i < 0:
        return None
    sel = res[i:i + len(cdr3b)]
    atoms = [a.coord for r in sel for a in r if a.element != "H"]
    return np.array(atoms) if atoms else None


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

def ensemble_contact(paths, cdr3b=None):
    """MEDIAN beta-peptide heavy-atom contact across all Protenix samples of one
    construct. With cdr3b given, contact is restricted to the CDR3beta loop residues
    (the specificity-determining loop, and the only part of chain B that reaches the
    peptide: a live profile of a real fold put 47 of 48 peptide contacts in CDR3beta);
    without it, the whole beta V-domain is used (CDR3beta-dominated but not identical).
    Pass cdr3b whenever the claim is literally "CDR3beta-peptide contact".

    Protenix emits several samples per seed whose docking pose varies wildly (MSA-free
    especially), so a single sample is not representative. The MEDIAN, not the mean, is
    the honest summary: it reflects the TYPICAL pose and is insensitive to a lone
    degenerate sample. Learned from the first live fold, where one scramble sample with
    591 spurious contacts (vs 0 in the other four) dominated the mean and inverted the
    cognate-vs-scramble verdict. Skips any model that does not parse to a full 5-chain
    TCR-pMHC (and, when cdr3b is given, any model whose CDR3beta cannot be located).
    Returns (median|None, n_models, n_valid)."""
    paths = list(paths)
    vals = [v for v in (score_model(p, cdr3b).get("cdr3_pep_atoms") for p in paths) if v is not None]
    if not vals:
        return None, len(paths), 0
    return float(np.median(vals)), len(paths), len(vals)


def score_model(cif_path, cdr3b=None) -> dict:
    chains = _heavy_by_chain(cif_path)
    if not {"A", "B", "C", "D", "E"}.issubset(chains):
        return {"n_chains": len(chains), "cdr3_pep_atoms": None, "crossing_angle": None, "dockq": None}
    pep = chains["E"]
    if cdr3b:
        beta = _cdr3b_beta_atoms(cif_path, cdr3b)
        if beta is None:  # CDR3beta not locatable -> loud None, never silently use whole chain
            return {"n_chains": len(chains), "cdr3_pep_atoms": None, "crossing_angle": None,
                    "dockq": None, "cdr3b_unmapped": True}
    else:
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

def score_pose(chains: dict):
    if "C" not in chains or "E" not in chains:
        return None
    mhc, pep = chains["C"], chains["E"]
    d = np.sqrt(((pep[:, None, :] - mhc[None, :, :]) ** 2).sum(-1))
    return float((d < 4.5).sum())


def mean_confidence(bfactors):
    if not bfactors:
        return None
    return float(sum(bfactors) / len(bfactors))


def verdict_groove(pose_atoms, clonotype_id: str, tool: str,
                   confidence=None) -> QCResult:
    """Honest peptide-groove verdict. Live calibration (2026-07-08) showed groove
    contact does NOT separate a cognate peptide from a scrambled non-binder: MHC-Fine
    seats any peptide in the groove just as deeply, and neither the contact count nor
    plddt discriminates. So there is no reliable/suspect split here. This reports only
    whether an in-groove pose was produced (placement), never a specificity claim."""
    if pose_atoms is None or pose_atoms <= 0:
        return QCResult(clonotype_id, "pose_failed",
                        "no in-groove peptide pose to score", tool=tool,
                        calibration_basis="pose_quality")
    reason = (f"in-groove pose ({pose_atoms:.0f} peptide-MHC contacts, "
              f"confidence {confidence}); this tool seats any peptide in the groove, "
              f"so it is placement only, not binding evidence")
    return QCResult(clonotype_id, "pose_only", reason,
                    tool=tool, calibration_basis="pose_quality")
