from rep2struct.tools import mhcfine_inputs, affinetune_inputs

FASTA = ">A\nAAAA\n>B\nBBBB\n>C\nCCCC\n>D\nDDDD\n>E\nSIINFEKL\n"


# tcrdock uses a gene-level TSV contract, not chain sequences; its builder is
# tested in test_tcrdock_inputs.py.


def test_mhcfine_drops_tcr_chains():
    out = mhcfine_inputs.build(FASTA)
    assert set(out["cognate"]["chains"]) == {"C", "D", "E"}
    assert "A" not in out["cognate"]["chains"]


def test_affinetune_maps_class_i_fields_and_scrambles():
    out = affinetune_inputs.build(FASTA)
    assert out["cognate"]["mhc"] == "CCCC" and out["cognate"]["b2m"] == "DDDD"
    assert out["cognate"]["peptide"] == "SIINFEKL"
    assert out["scramble"]["peptide"] != "SIINFEKL"
