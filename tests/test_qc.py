from pathlib import Path
from rep2struct.qc import score_model, verdict, verdict_binding
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


def test_binding_verdict_presented_and_not_presented():
    hi = verdict_binding(0.9, 0.5, "c1", tool="affinetune")
    lo = verdict_binding(0.3, 0.5, "c2", tool="affinetune")
    assert hi.qc_verdict == "presented" and lo.qc_verdict == "not_presented"


def test_binding_verdict_is_honest_not_a_fold():
    hi = verdict_binding(0.9, 0.5, "c1", tool="affinetune")
    lo = verdict_binding(0.3, 0.5, "c2", tool="affinetune")
    for r in (hi, lo):
        assert "presentation" in r.reason.lower()
        assert "fold" not in r.reason.lower() and "structure" not in r.reason.lower()
    assert hi.tool == "affinetune" and hi.calibration_basis == "binding_score_null"
    assert hi.cdr3_pep_atoms is None and hi.dockq is None


def test_binding_verdict_boundary_equality_is_not_presented():
    assert verdict_binding(0.5, 0.5, "c1", tool="affinetune").qc_verdict == "not_presented"


import numpy as np
from rep2struct.qc import common_checks


def _chain(*xyz):
    return np.array(xyz, dtype=float)


def test_common_checks_passes_a_sane_two_chain_model():
    chains = {"C": _chain([0, 0, 0], [10, 0, 0]), "E": _chain([5, 0, 0], [6, 0, 0])}
    r = common_checks(chains, expected={"C", "E"})
    assert r["ok"] and r["issues"] == [] and r["has_peptide"] and r["n_chains"] == 2


def test_common_checks_flags_missing_chain():
    chains = {"C": _chain([0, 0, 0])}
    r = common_checks(chains, expected={"C", "E"})
    assert not r["ok"] and any("missing" in i for i in r["issues"]) and not r["has_peptide"]


def test_common_checks_flags_nonfinite_coords():
    chains = {"C": _chain([0, 0, 0]), "E": _chain([np.nan, 0, 0])}
    r = common_checks(chains, expected={"C", "E"})
    assert not r["ok"] and any("non-finite" in i for i in r["issues"])


def test_common_checks_flags_severe_clash():
    # two different chains with atoms 0.1A apart
    chains = {"C": _chain([0, 0, 0]), "E": _chain([0.1, 0, 0])}
    r = common_checks(chains, expected={"C", "E"})
    assert not r["ok"] and any("clash" in i for i in r["issues"])


from rep2struct.qc import score_pose, mean_confidence, verdict_groove


def test_score_pose_counts_peptide_mhc_contacts():
    chains = {"C": _chain([0, 0, 0], [100, 0, 0]), "E": _chain([1, 0, 0], [50, 0, 0])}
    # only the E atom at (1,0,0) is within 4.5A of a C atom
    assert score_pose(chains) == 1.0
    assert score_pose({"E": _chain([0, 0, 0])}) is None  # no MHC heavy chain


def test_mean_confidence():
    assert mean_confidence([80.0, 90.0]) == 85.0
    assert mean_confidence(None) is None
    assert mean_confidence([]) is None


def test_verdict_groove_is_a_pose_not_a_fold():
    hi = verdict_groove(20.0, 10.0, "c1", tool="mhcfine", confidence=88.0)
    lo = verdict_groove(5.0, 10.0, "c2", tool="mhcfine")
    assert hi.qc_verdict == "pose_reliable" and lo.qc_verdict == "pose_suspect"
    for r in (hi, lo):
        assert "pose" in r.reason.lower()
        for bad in ("fold", "structure", "recognition"):
            assert bad not in r.reason.lower()
    assert hi.tool == "mhcfine" and hi.calibration_basis == "groove_scramble_null"
    none_verdict = verdict_groove(None, 10.0, "c3", tool="mhcfine")
    assert none_verdict.qc_verdict == "pose_failed"
    assert "pose" in none_verdict.reason.lower()
