from rep2struct.tools import mhcfine_inputs, affinetune_inputs

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
