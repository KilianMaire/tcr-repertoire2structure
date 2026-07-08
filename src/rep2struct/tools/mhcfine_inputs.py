from __future__ import annotations
from .construct_io import parse_fasta, scramble_peptide

TOOL = "mhcfine"


def build(construct_fasta: str) -> dict:
    """pMHC (class I) input for MHC-Fine (cognate + scramble control).

    MHC-Fine takes a single MHC protein sequence plus the peptide (no TCR, and
    it adds b2m internally), so this emits protein_sequence + peptide_sequence,
    not chain sequences. Source is the shared construct FASTA: chain C is the
    MHC heavy chain, chain E the peptide. The scramble shuffles only the peptide,
    giving mhcfine its own pose-in-groove calibration null.
    """
    chains = parse_fasta(construct_fasta)
    protein, peptide = chains["C"], chains["E"]
    return {
        "cognate": {"protein_sequence": protein, "peptide_sequence": peptide},
        "scramble": {"protein_sequence": protein,
                     "peptide_sequence": scramble_peptide(peptide)},
    }
