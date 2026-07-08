from pathlib import Path
from rep2struct.qc import score_model, verdict
FIX = Path(__file__).parent / "fixtures"

def test_cognate_is_reliable():
    s = score_model(FIX / "cognate_min.cif")
    r = verdict(s, scramble_threshold=1.0)
    assert s["cdr3_pep_atoms"] > 1.0
    assert r.qc_verdict == "reliable"

def test_scramble_is_suspect():
    s = score_model(FIX / "scramble_min.cif")
    r = verdict(s, scramble_threshold=1.0)
    assert r.qc_verdict == "suspect"

def test_threshold_boundary_is_suspect():
    r = verdict({"cdr3_pep_atoms": 5.0, "clonotype_id": "b"}, scramble_threshold=5.0)
    assert r.qc_verdict == "suspect"

def test_three_chain_model_is_qc_failed():
    s = score_model(FIX / "threechain_min.cif")
    assert s["cdr3_pep_atoms"] is None
    r = verdict(s, scramble_threshold=1.0)
    assert r.qc_verdict == "qc_failed"
