import pytest
from rep2struct.tools import mhcfine_inputs, affinetune_inputs, protenix_inputs

FASTA = ">A\nAAAA\n>B\nBBBB\n>C\nCCCC\n>D\nDDDD\n>E\nSIINFEKL\n"


# tcrdock uses a gene-level TSV contract, not chain sequences; its builder is
# tested in test_tcrdock_inputs.py.


def test_mhcfine_emits_protein_and_peptide_only():
    out = mhcfine_inputs.build(FASTA)
    # mhcfine takes an MHC protein sequence + peptide, not chain sequences
    assert out["cognate"]["protein_sequence"] == "CCCC"   # chain C = MHC heavy
    assert out["cognate"]["peptide_sequence"] == "SIINFEKL"
    assert "chains" not in out["cognate"]
    # TCR (A/B) and b2m (D) are not mhcfine inputs
    assert set(out["cognate"]) == {"protein_sequence", "peptide_sequence"}


def test_mhcfine_scrambles_only_peptide():
    out = mhcfine_inputs.build(FASTA)
    c, s = out["cognate"], out["scramble"]
    assert s["protein_sequence"] == c["protein_sequence"]
    assert s["peptide_sequence"] != c["peptide_sequence"]
    assert sorted(s["peptide_sequence"]) == sorted(c["peptide_sequence"])


def test_affinetune_maps_class_i_fields_and_scrambles():
    out = affinetune_inputs.build(FASTA)
    assert out["cognate"]["mhc"] == "CCCC" and out["cognate"]["b2m"] == "DDDD"
    assert out["cognate"]["peptide"] == "SIINFEKL"
    assert out["scramble"]["peptide"] != "SIINFEKL"


def _protenix_chain(record, chain_id):
    seqs = record[0]["sequences"]
    return next(s["proteinChain"]["sequence"] for s in seqs
               if s["proteinChain"]["id"] == [chain_id])


def test_protenix_emits_all_five_chains():
    out = protenix_inputs.build(FASTA)
    cog = out["cognate"]
    assert [_protenix_chain(cog, c) for c in "ABCDE"] == ["AAAA", "BBBB", "CCCC", "DDDD", "SIINFEKL"]
    assert cog[0]["covalent_bonds"] == []


def test_protenix_scrambles_only_the_peptide():
    out = protenix_inputs.build(FASTA)
    cog, scr = out["cognate"], out["scramble"]
    # A-D identical, only E (peptide) shuffled to the same multiset
    for c in "ABCD":
        assert _protenix_chain(scr, c) == _protenix_chain(cog, c)
    assert _protenix_chain(scr, "E") != "SIINFEKL"
    assert sorted(_protenix_chain(scr, "E")) == sorted("SIINFEKL")


def test_protenix_raises_on_missing_chain():
    with pytest.raises(ValueError):
        protenix_inputs.build(">A\nAAAA\n>E\nSIINFEKL\n")  # missing B, C, D
