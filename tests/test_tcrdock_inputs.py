from rep2struct.schema import Clonotype, Annotation
from rep2struct.tools import tcrdock_inputs


def _clon():
    return Clonotype(id="c1", trav="TRAV1-2", cdr3a="CAVMDSSYKLIF",
                     trbv="TRBV19", cdr3b="CASSIRSSYEQYF", size=5,
                     trav_allele="TRAV1-2*01", trbv_allele="TRBV19*01",
                     traj="TRAJ33", trbj="TRBJ2-1")


def _ann():
    return Annotation(clonotype_id="c1", annotatable=True, confidence_tier="high",
                      epitope="GILGFVFTL", hla="HLA-A*02:01")


def test_tcrdock_emits_gene_level_tsv_row():
    out = tcrdock_inputs.build(_clon(), _ann())
    r = out["cognate"]["row"]
    assert r["organism"] == "human" and r["mhc_class"] == 1
    assert r["mhc"] == "A*02:01"                     # "HLA-" prefix stripped
    assert r["va"] == "TRAV1-2*01" and r["ja"] == "TRAJ33"
    assert r["vb"] == "TRBV19*01" and r["jb"] == "TRBJ2-1"
    assert r["cdr3a"] == "CAVMDSSYKLIF" and r["cdr3b"] == "CASSIRSSYEQYF"
    assert r["peptide"] == "GILGFVFTL"


def test_tcrdock_scramble_shuffles_only_peptide():
    out = tcrdock_inputs.build(_clon(), _ann())
    c, s = out["cognate"]["row"], out["scramble"]["row"]
    assert s["peptide"] != c["peptide"]
    assert sorted(s["peptide"]) == sorted(c["peptide"])
    # every non-peptide field is identical between cognate and scramble
    assert {k: v for k, v in s.items() if k != "peptide"} == \
           {k: v for k, v in c.items() if k != "peptide"}


def test_tcrdock_falls_back_to_bare_gene_without_allele():
    c = _clon()
    c.trav_allele = None
    c.trbv_allele = None
    out = tcrdock_inputs.build(c, _ann())
    r = out["cognate"]["row"]
    assert r["va"] == "TRAV1-2" and r["vb"] == "TRBV19"


def test_tcrdock_bare_hla_passes_through():
    a = _ann()
    a.hla = "B*07:02"
    out = tcrdock_inputs.build(_clon(), a)
    assert out["cognate"]["row"]["mhc"] == "B*07:02"


def test_tcrdock_normalizes_multivalue_hla():
    # real reference DB returns comma-joined multi-value HLA; keep the specific one
    a = _ann()
    a.hla = "HLA-A*02,HLA-A*02:01"
    out = tcrdock_inputs.build(_clon(), a)
    assert out["cognate"]["row"]["mhc"] == "A*02:01"


def test_tcrdock_bare_single_field_hla_survives():
    a = _ann()
    a.hla = "HLA-A*02"          # no 2-field token available: keep what we have
    out = tcrdock_inputs.build(_clon(), a)
    assert out["cognate"]["row"]["mhc"] == "A*02"


def test_tcrdock_row_keys_match_declared_columns():
    out = tcrdock_inputs.build(_clon(), _ann())
    assert list(out["cognate"]["row"].keys()) == tcrdock_inputs.COLUMNS
    assert list(out["scramble"]["row"].keys()) == tcrdock_inputs.COLUMNS


def test_tcrdock_rejects_missing_peptide_or_hla():
    import pytest
    a = _ann()
    a.epitope = None
    with pytest.raises(ValueError):
        tcrdock_inputs.build(_clon(), a)
    a2 = _ann()
    a2.hla = None
    with pytest.raises(ValueError):
        tcrdock_inputs.build(_clon(), a2)
