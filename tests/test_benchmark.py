from pathlib import Path
from rep2struct import benchmark as bm
from rep2struct.schema import Clonotype, Annotation

FIX = Path(__file__).parent / "fixtures"

def test_is_novel():
    assert bm.is_novel(None) is True
    assert bm.is_novel(0.0) is False
    assert bm.is_novel(1.0) is False
    assert bm.is_novel(1.5) is True

def test_panel_epitopes_sorted_unique():
    truth = {"c1": ("GILGFVFTL", "HLA-A*02:01"),
             "c2": ("GILGFVFTL", "HLA-A*02:01"),
             "c3": ("NLVPMVATV", "HLA-A*02:01")}
    assert bm.panel_epitopes(truth) == [
        ("GILGFVFTL", "HLA-A*02:01"), ("NLVPMVATV", "HLA-A*02:01")]

def test_per_hla_novel_counts():
    clonos = [Clonotype("c1","TRAV1","CAA","TRBV1","CBB",5),
              Clonotype("c2","TRAV1","CAC","TRBV1","CBD",3)]
    truth = {"c1": ("GILGFVFTL","HLA-A*02:01"),
             "c2": ("GILGFVFTL","HLA-A*02:01")}
    anns = [Annotation("c1", False, "unannotatable", tcrdist=None),
            Annotation("c2", True, "high", tcrdist=0.0)]
    out = bm.per_hla_novel_counts(clonos, truth, anns)
    assert out["HLA-A*02:01"]["n_total"] == 2
    assert out["HLA-A*02:01"]["n_novel"] == 1
    assert out["HLA-A*02:01"]["epitopes"]["GILGFVFTL"]["n_novel"] == 1

def test_decoys_same_hla_first():
    panel = [("GILGFVFTL","HLA-A*02:01"), ("NLVPMVATV","HLA-A*02:01"),
             ("GLCTLVAML","HLA-A*02:01"), ("KLGGALQAK","HLA-A*03:01")]
    d = bm.decoys_for("GILGFVFTL", "HLA-A*02:01", panel, k=2)
    assert d == ["GLCTLVAML", "NLVPMVATV"]
    assert "GILGFVFTL" not in d

def test_decoys_same_hla_only_no_cross_fill():
    panel = [("GILGFVFTL","HLA-A*02:01"), ("KLGGALQAK","HLA-A*03:01")]
    assert bm.decoys_for("GILGFVFTL", "HLA-A*02:01", panel, k=3) == []

def test_scramble_preserves_composition_and_differs():
    s = bm.scramble_peptide("GILGFVFTL", seed=1)
    assert sorted(s) == sorted("GILGFVFTL") and s != "GILGFVFTL"

def test_build_panel_constructs_keys_and_peptides():
    clono = Clonotype("c1","TRAV1","CAA","TRBV1","CBB",5)
    tcr_seqs = {"c1": {"A": "AAAA", "B": "BBBB"}}
    mhc_seqs = {"HLA-A*02:01": {"heavy": "HHHH", "b2m": "MMMM"}}
    jobs = bm.build_panel_constructs(clono, "GILGFVFTL", "HLA-A*02:01",
                                     ["NLVPMVATV"], tcr_seqs, mhc_seqs)
    assert set(jobs) == {"GILGFVFTL", "NLVPMVATV", "__scramble__"}
    assert ">E\nGILGFVFTL" in jobs["GILGFVFTL"].construct_fasta
    scr_pep = jobs["__scramble__"].construct_fasta.split(">E\n")[1].split("\n")[0]
    assert sorted(scr_pep) == sorted("GILGFVFTL")
    assert jobs["GILGFVFTL"].clonotype_id == "c1"

def test_contact_and_retrieval_with_fixtures():
    paths = {"COGNATE": [str(FIX/"cognate_min.cif")],
             "DECOY":   [str(FIX/"scramble_min.cif")]}
    contacts = bm.contact_by_epitope(paths)
    assert contacts["COGNATE"] >= contacts["DECOY"]
    res = bm.retrieval_result(contacts, "COGNATE")
    assert res["ranked"][0] == "COGNATE"
    assert res["top1"] is True

def test_auroc_pairs():
    assert bm.auroc([(10.0, [1.0, 2.0])]) == 1.0
    assert bm.auroc([(1.0, [10.0])]) == 0.0
    assert bm.auroc([(5.0, [5.0])]) == 0.5
    assert bm.auroc([]) is None

def test_retrieval_none_cognate_is_not_a_win():
    assert bm.retrieval_result({"COG": None, "D": 5.0}, "COG")["top1"] is False

def test_retrieval_all_none_is_not_a_win():
    assert bm.retrieval_result({"COG": None, "D": None}, "COG")["top1"] is False

def test_retrieval_tie_is_not_a_win():
    assert bm.retrieval_result({"COG": 5.0, "D": 5.0}, "COG")["top1"] is False

def test_retrieval_excludes_scramble_key():
    r = bm.retrieval_result({"COG": 5.0, "D": 1.0, "__scramble__": 99.0}, "COG")
    assert r["top1"] is True and "__scramble__" not in r["ranked"]

def test_cdr3b_plddt_none_when_substring_absent():
    assert bm.model_cdr3b_plddt(str(FIX/"cognate_min.cif"), "AAAA", "ZZZZ") is None

def test_cdr3b_plddt_returns_float_or_none_when_located():
    v = bm.model_cdr3b_plddt(str(FIX/"cognate_min.cif"), "BBBB", "BB")
    assert v is None or isinstance(v, float)

def test_model_cdr3b_plddt_none_when_no_chain_b_seq():
    assert bm.model_cdr3b_plddt(str(FIX/"cognate_min.cif"), "", "CASS") is None

def test_sequence_baseline_top1():
    assert bm.sequence_baseline_top1("GILGFVFTL", "GILGFVFTL") is True
    assert bm.sequence_baseline_top1("NLVPMVATV", "GILGFVFTL") is False
    assert bm.sequence_baseline_top1(None, "GILGFVFTL") is False

def test_bootstrap_ci_all_hits():
    pt, lo, hi = bm.bootstrap_ci([True]*20, n_boot=500, seed=1)
    assert pt == 1.0 and lo == 1.0 and hi == 1.0

def test_bootstrap_ci_bounds_order():
    pt, lo, hi = bm.bootstrap_ci([True, False, True, False], n_boot=500, seed=1)
    assert 0.0 <= lo <= pt <= hi <= 1.0

def test_permutation_p_strong_signal():
    assert bm.permutation_p([True]*10, chance=0.25, n_perm=5000, seed=1) < 0.01
    assert bm.permutation_p([False]*10, chance=0.25, n_perm=5000, seed=1) > 0.5

def test_tcr_blind_accuracy_constant_predictor():
    contacts = [{"A": 9.0, "B": 1.0}, {"A": 8.0, "B": 2.0}]
    assert bm.tcr_blind_prediction(contacts) == "A"
    assert bm.tcr_blind_accuracy(contacts, ["A", "B"]) == 0.5

def test_label_permutation_p_runs():
    contacts = [{"A": 9.0, "B": 1.0}, {"B": 9.0, "A": 1.0}]
    p = bm.label_permutation_p(1.0, contacts, ["A", "B"], n_perm=1000, seed=1)
    assert 0.0 <= p <= 1.0

def test_paired_contrast_cognate_higher():
    r = bm.paired_contrast([(10.0, 2.0), (8.0, 3.0)])
    assert r["frac_cognate_higher"] == 1.0 and r["mean_delta"] > 0

def test_paired_contrast_empty():
    assert bm.paired_contrast([(None, 1.0)])["n"] == 0
