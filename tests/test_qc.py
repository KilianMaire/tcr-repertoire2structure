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


def test_ensemble_contact_summarizes_valid_models_and_skips_bad():
    from rep2struct.qc import ensemble_contact
    single = score_model(FIX / "cognate_min.cif")["cdr3_pep_atoms"]
    # two copies of the same cognate -> median equals the single-model contact
    med, n_models, n_valid = ensemble_contact([str(FIX / "cognate_min.cif")] * 2)
    assert med == single and n_models == 2 and n_valid == 2
    # a 3-chain model is not a valid TCR-pMHC and is skipped, not counted
    med2, n_models2, n_valid2 = ensemble_contact(
        [str(FIX / "cognate_min.cif"), str(FIX / "threechain_min.cif")])
    assert med2 == single and n_models2 == 2 and n_valid2 == 1


def test_ensemble_contact_median_ignores_a_degenerate_outlier(monkeypatch):
    # Regression from the first live fold: scramble samples were [0, 591, 0, 0, 0]; the mean
    # (118) was dominated by the lone 591 pose and flipped the verdict. The MEDIAN reflects
    # the typical pose (0) and ignores the outlier.
    from rep2struct import qc
    scores = {"s0": 0.0, "s1": 591.0, "s2": 0.0, "s3": 0.0, "s4": 0.0}
    monkeypatch.setattr(qc, "score_model", lambda p: {"cdr3_pep_atoms": scores[p], "n_chains": 5})
    med, n_models, n_valid = qc.ensemble_contact(list(scores))
    assert n_models == 5 and n_valid == 5
    assert med == 0.0                                  # median ignores the lone 591
    assert sum(scores.values()) / 5 > 100              # the mean would have been dominated by it


def test_ensemble_contact_all_invalid_is_none():
    from rep2struct.qc import ensemble_contact
    mean, n_models, n_valid = ensemble_contact([str(FIX / "threechain_min.cif")])
    assert mean is None and n_models == 1 and n_valid == 0


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


def test_verdict_groove_is_pose_only_never_specificity():
    # Honest semantics: any in-groove pose is "pose_only" (placement), regardless of
    # contact count. There is NO reliable/suspect split, because live calibration
    # showed groove contact does not separate a binder from a scrambled non-binder.
    v = verdict_groove(20.0, "c1", tool="mhcfine", confidence=88.0)
    assert v.qc_verdict == "pose_only"
    assert "pose" in v.reason.lower() and "placement" in v.reason.lower()
    for bad in ("fold", "structure"):
        assert bad not in v.reason.lower()
    assert v.tool == "mhcfine" and v.calibration_basis == "pose_quality"
    # a high-contact scramble-like pose gets the SAME verdict as a low-contact one
    assert verdict_groove(5.0, "c1b", tool="mhcfine").qc_verdict == "pose_only"
    # no pose / peptide not in the groove -> honest failure
    for empty in (None, 0.0):
        f = verdict_groove(empty, "c2", tool="mhcfine")
        assert f.qc_verdict == "pose_failed"
        assert "pose" in f.reason.lower()
        assert f.calibration_basis == "pose_quality"
