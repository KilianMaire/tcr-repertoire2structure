from rep2struct.report import render_report
from rep2struct.schema import Clonotype, Annotation, QCResult

def test_report_contains_rows_and_verdicts():
    clons = [Clonotype("c1", "TRAV1-2", "CAVA", "TRBV19", "CASSB", 5)]
    anns = [Annotation("c1", True, "high", 3.0, epitope="GILGFVFTL", hla="HLA-A*02:01")]
    qcs = [QCResult("c1", "reliable", "ok", cdr3_pep_atoms=12.0)]
    html = render_report(clons, anns, qcs, metrics={"precision": 0.8, "recall": 0.6,
                                                    "unannotatable_rate": 0.3, "n": 10})
    assert "GILGFVFTL" in html
    assert "reliable" in html
    assert "0.8" in html  # validation block rendered
    assert html.lstrip().lower().startswith("<!doctype html")

def test_unannotatable_is_shown():
    clons = [Clonotype("c2", "TRAV1", "CAVA", "TRBV2", "CASSX", 2)]
    anns = [Annotation("c2", False, "unannotatable")]
    qcs = []
    html = render_report(clons, anns, qcs)
    assert "unannotatable" in html
