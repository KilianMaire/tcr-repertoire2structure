from __future__ import annotations
from .construct_io import parse_fasta, scramble_peptide, pmhc_only

TOOL = "mhcfine"


def build(construct_fasta: str) -> dict:
    chains = pmhc_only(parse_fasta(construct_fasta))
    sc = dict(chains)
    sc["E"] = scramble_peptide(chains["E"])
    return {"cognate": {"chains": chains}, "scramble": {"chains": sc}}
