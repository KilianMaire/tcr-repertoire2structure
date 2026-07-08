from __future__ import annotations
from .construct_io import parse_fasta, scramble_peptide

TOOL = "affinetune"


def build(construct_fasta: str) -> dict:
    chains = parse_fasta(construct_fasta)
    def rec(pep):
        return {"mhc": chains["C"], "b2m": chains["D"], "peptide": pep}
    return {"cognate": rec(chains["E"]),
            "scramble": rec(scramble_peptide(chains["E"]))}
