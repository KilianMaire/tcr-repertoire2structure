from __future__ import annotations
from .construct_io import scramble_peptide

TOOL = "tcrdock"

# tcrdock reconstructs the TCR:pMHC structure from gene names and CDR3 loops,
# so it consumes a 10-column TSV row, NOT raw chain sequences. These are the
# exact column names its setup_for_alphafold.py expects.
COLUMNS = ["organism", "mhc_class", "mhc", "peptide",
           "va", "ja", "cdr3a", "vb", "jb", "cdr3b"]


def _mhc_allele(hla):
    # tcrdock wants a bare class I allele name, e.g. "A*02:01" (no "HLA-" prefix).
    return hla.replace("HLA-", "").strip() if hla else hla


def _row(clonotype, annotation, peptide) -> dict:
    return {
        "organism": "human",
        "mhc_class": 1,
        "mhc": _mhc_allele(annotation.hla),
        "peptide": peptide,
        # allele-qualified gene when tcr_explorer resolved one, else the bare gene
        "va": clonotype.trav_allele or clonotype.trav,
        "ja": clonotype.traj,
        "cdr3a": clonotype.cdr3a,
        "vb": clonotype.trbv_allele or clonotype.trbv,
        "jb": clonotype.trbj,
        "cdr3b": clonotype.cdr3b,
    }


def build(clonotype, annotation) -> dict:
    """Gene-level TCR:pMHC construct for tcrdock (cognate + scramble control).

    Emits the TSV row tcrdock consumes (gene names + CDR3 + peptide + MHC
    allele), not chain sequences: tcrdock builds the structure from its own
    templates. The scramble row is identical except the peptide is shuffled,
    giving tcrdock its own calibration null (never shared with another tool).
    """
    pep = annotation.epitope
    # A class I structural row with no peptide or no MHC allele is not a valid
    # tcrdock target; fail loud rather than emit a silent mhc=None / crash row.
    if not pep or not annotation.hla:
        raise ValueError(
            f"tcrdock needs both peptide and HLA; got epitope={annotation.epitope!r} "
            f"hla={annotation.hla!r} for {annotation.clonotype_id}")
    row = {
        "cognate": {"row": _row(clonotype, annotation, pep)},
        "scramble": {"row": _row(clonotype, annotation, scramble_peptide(pep))},
    }
    assert list(row["cognate"]["row"].keys()) == COLUMNS  # emitted contract == declared
    return row
