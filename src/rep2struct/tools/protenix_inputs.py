from __future__ import annotations
from .construct_io import parse_fasta, scramble_peptide

TOOL = "protenix"
_ORDER = ["A", "B", "C", "D", "E"]


def _to_protenix(name: str, chains: dict) -> list:
    """One Protenix prediction record: chains A-E as protein chains, no covalent
    bonds. Mirrors science/scripts/build_protenix_inputs.py, the recipe that produced the
    validated TABLO folds."""
    seqs = [{"proteinChain": {"sequence": chains[c], "count": 1, "id": [c]}}
            for c in _ORDER if c in chains]
    return [{"name": name, "sequences": seqs, "covalent_bonds": []}]


def build(construct_fasta: str) -> dict:
    """Protenix inputs for the full TCR-pMHC fold (cognate + scramble control).

    Protenix folds the whole A-E construct, so this emits the Protenix JSON for the
    cognate and a scramble control that shuffles only the peptide (chain E),
    everything else identical. The scramble is the per-clonotype null the cdr3-peptide
    QC calibrates against (Honesty Rule 2: a fold does not confirm specificity; a
    cognate must beat its own scramble on CDR3-peptide contact).
    """
    chains = parse_fasta(construct_fasta)
    missing = set("ABCDE") - set(chains)
    if missing:
        raise ValueError(f"protenix construct missing chains {sorted(missing)}")
    sc = dict(chains)
    sc["E"] = scramble_peptide(chains["E"])
    return {"cognate": _to_protenix("cognate", chains),
            "scramble": _to_protenix("scramble", sc)}
