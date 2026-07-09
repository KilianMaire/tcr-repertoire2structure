from rep2struct import benchmark as bm
from rep2struct.schema import Clonotype, Annotation

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
