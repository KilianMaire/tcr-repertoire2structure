import re

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


def _fixtures():
    c = Clonotype(id="c1", trav="TRAV1", cdr3a="CAA", trbv="TRBV2", cdr3b="CAB", size=5)
    a = Annotation(clonotype_id="c1", annotatable=True, confidence_tier="high",
                   epitope="SIINFEKL", hla="A*02:01")
    return [c], [a]


def test_binding_row_is_labelled_predicted_presentation():
    c, a = _fixtures()
    q = QCResult("c1", "presented", "predicted presentation above the score null",
                 tool="affinetune", calibration_basis="binding_score_null")
    html = render_report(c, a, [q])
    assert "predicted presentation" in html.lower()
    assert "affinetune" in html


def test_binding_row_not_called_a_structure_or_fold():
    c, a = _fixtures()
    q = QCResult("c1", "presented", "predicted presentation above the score null",
                 tool="affinetune")
    html = render_report(c, a, [q])
    row = [ln for ln in html.splitlines() if "affinetune" in ln]
    assert row and all("fold" not in ln.lower() and "structure" not in ln.lower() for ln in row)


def test_msa_free_row_flagged_reduced_confidence():
    c, a = _fixtures()   # reuse the existing helper in this file
    q = QCResult("c1", "reliable", "ok", tool="protenix", calibration_basis="scramble_null")
    html = render_report(c, a, [q], msa_basis={"c1": "none"})
    assert "reduced confidence" in html.lower()


def test_render_report_msa_basis_is_optional():
    c, a = _fixtures()
    q = QCResult("c1", "reliable", "ok", tool="protenix")
    render_report(c, a, [q])  # must not raise without msa_basis


def test_pose_row_labelled_pose_not_fold_or_structure():
    c, a = _fixtures()
    q = QCResult("c1", "pose_reliable",
                 "peptide in groove contact 20 beats scramble null; model confidence 88.0",
                 tool="mhcfine", calibration_basis="groove_scramble_null")
    html = render_report(c, a, [q])
    # match the whole table row for this clonotype, independent of template line layout
    row = re.search(r"<tr[^>]*>.*?mhcfine.*?</tr>", html, re.S)
    assert row, "no table row rendered for the mhcfine result"
    block = row.group(0).lower()
    assert "pose" in block
    assert "fold" not in block and "structure" not in block


def test_validity_column_renders_when_supplied():
    c, a = _fixtures()
    q = QCResult("c1", "reliable", "ok", tool="protenix")
    html = render_report(c, a, [q], validity={"c1": "valid"})
    assert "valid" in html.lower()


def test_render_report_validity_is_optional():
    c, a = _fixtures()
    q = QCResult("c1", "reliable", "ok", tool="protenix")
    render_report(c, a, [q])  # must not raise without validity


def test_pose_failed_row_labelled_pose_not_structure():
    c, a = _fixtures()
    q = QCResult("c1", "pose_failed", "no MHC-peptide pose to score", tool="mhcfine")
    html = render_report(c, a, [q])
    row = re.search(r"<tr[^>]*>.*?mhcfine.*?</tr>", html, re.S)
    assert row, "no table row rendered for the mhcfine result"
    block = row.group(0).lower()
    assert "pose" in block
    assert "structure" not in block and "fold" not in block


def test_qc_failed_row_for_structure_tools_still_says_structure():
    c, a = _fixtures()
    q = QCResult("c1", "qc_failed", "model has 3 chains, need 5", tool="protenix")
    html = render_report(c, a, [q])
    row = re.search(r"<tr[^>]*>.*?protenix.*?</tr>", html, re.S)
    assert row, "no table row rendered for the protenix result"
    assert "structure (qc failed)" in row.group(0).lower()
