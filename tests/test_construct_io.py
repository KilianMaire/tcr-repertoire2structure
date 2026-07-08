from rep2struct.tools.construct_io import parse_fasta, scramble_peptide, pmhc_only


def test_parse_fasta_roundtrip():
    c = parse_fasta(">A\nAAA\n>E\nSII\n")
    assert c == {"A": "AAA", "E": "SII"}


def test_scramble_is_deterministic_and_non_identity():
    assert scramble_peptide("SIINFEKL") == scramble_peptide("SIINFEKL")
    assert scramble_peptide("SIINFEKL") != "SIINFEKL"
    assert sorted(scramble_peptide("SIINFEKL")) == sorted("SIINFEKL")  # same composition


def test_pmhc_only_drops_tcr_chains():
    chains = {"A": "a", "B": "b", "C": "c", "D": "d", "E": "e"}
    assert pmhc_only(chains) == {"C": "c", "D": "d", "E": "e"}

def test_normalize_hla_collapses_multivalue():
    from rep2struct.tools.construct_io import normalize_hla
    assert normalize_hla("HLA-A*02,HLA-A*02:01") == "HLA-A*02:01"
    assert normalize_hla("A*03:01") == "HLA-A*03:01"        # prefix added
    assert normalize_hla("HLA-A*02:01") == "HLA-A*02:01"    # idempotent
    assert normalize_hla(None) is None
