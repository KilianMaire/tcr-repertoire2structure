import asyncio
from pathlib import Path
from rep2struct import agent_tools as at

FIX = Path(__file__).parent / "fixtures"
_ROOT = Path(__file__).resolve().parents[1]


def _run(coro):
    return asyncio.run(coro)


# The @tool decorator (claude_agent_sdk 0.2.113) wraps the async function in an
# SdkMcpTool dataclass (name, description, input_schema, handler). The plain
# async function passed to build_server() lives at `.handler`; call that
# directly with the args dict to run the real tool logic offline.
def _call(sdk_tool, args):
    return _run(sdk_tool.handler(args))


def test_ingest_then_annotate_tools(tmp_path):
    def sim(cdr3_a, v_a, cdr3_b, v_b, species="human", top_k=5):
        if cdr3_b == "CASSIRSSYEQYF":
            return ([{"epitope": "GILGFVFTL", "mhc": "HLA-A*02:01",
                      "antigen": "Flu M1", "distance": 3.0}], "tcrdist", 1, [])
        return ([], "tcrdist", 0, [])
    at.configure(sim_fn=sim, assign_fn=lambda g, species="human", chain=None: g + "*01")
    rd = str(tmp_path / "run")
    ing = _call(at.ingest_repertoire, {"run_dir": rd, "csv_path": str(FIX / "tenx_tiny.csv")})
    assert ing["structuredContent"]["clonotypes"] >= 1
    ann = _call(at.annotate_specificity, {"run_dir": rd})
    tiers = ann["structuredContent"]["tiers"]
    assert tiers.get("high", 0) >= 1


def test_prep_and_select_skips_unfoldable(tmp_path):
    # One clonotype is annotatable with an hla, the other is unannotatable
    # (hla=None). build_construct would KeyError on mhc_seqs[None] for the
    # unannotatable one if prep_and_select didn't filter to foldable entries
    # first (the same guard as Task 10's pipeline.run_pipeline).
    def sim(cdr3_a, v_a, cdr3_b, v_b, species="human", top_k=5):
        if cdr3_b == "CASSIRSSYEQYF":
            return ([{"epitope": "GILGFVFTL", "mhc": "HLA-A*02:01",
                      "antigen": "Flu M1", "distance": 3.0}], "tcrdist", 1, [])
        return ([], "tcrdist", 0, [])
    at.configure(sim_fn=sim, assign_fn=lambda g, species="human", chain=None: g + "*01")
    rd = str(tmp_path / "run3")
    _call(at.ingest_repertoire, {"run_dir": rd, "csv_path": str(FIX / "tenx_tiny.csv")})
    _call(at.annotate_specificity, {"run_dir": rd})
    out = _call(at.prep_and_select, {"run_dir": rd, "top_n": 10})
    jobs = out["structuredContent"]["jobs"]
    listed = _call(at.list_fold_jobs, {"run_dir": rd})
    assert len(jobs) == len(listed["structuredContent"]["jobs"])
    # every prepared job must have come from a foldable (annotatable + hla) entry
    assert len(jobs) >= 1


def test_qc_tool_flags_scramble(tmp_path):
    rd = str(tmp_path / "run2")
    at.configure()
    # seed a fold job + model path pointing at the scramble fixture
    _call(at.ingest_repertoire, {"run_dir": rd, "csv_path": str(FIX / "tenx_tiny.csv")})
    _call(at.record_fold_result, {"run_dir": rd, "clonotype_id": "z",
                                  "model_paths": [str(FIX / "scramble_min.cif")]})
    out = _call(at.qc_structure, {"run_dir": rd, "clonotype_id": "z", "scramble_threshold": 1.0})
    assert out["structuredContent"]["qc_verdict"] == "suspect"


def _clean_five_chains():
    # a well-separated synthetic 5-chain TCR-pMHC so the common-gate passes; the ensemble
    # contact itself is read from the real CIF fixtures, this only satisfies the validity gate.
    import numpy as np
    return {c: np.array([[i * 40.0, 0, 0], [i * 40.0 + 10, 0, 0]], float)
            for i, c in enumerate("ABCDE")}


def test_qc_protenix_ensembles_and_uses_own_scramble_null(tmp_path, monkeypatch):
    # Protenix (cdr3_peptide) records a cognate + scramble pair per clonotype, named
    # {cid}_cognate / {cid}_scramble. QC must ensemble the cognate samples and calibrate
    # against this clonotype's OWN scramble ensemble, ignoring the caller's scalar.
    import shutil
    from rep2struct import qc
    from rep2struct.runstate import RunState
    monkeypatch.setattr(qc, "load_chains", lambda p: _clean_five_chains())  # bypass the toy-fixture clash gate
    rd = str(tmp_path / "run")
    cog = tmp_path / "z_cognate.cif"; shutil.copy(FIX / "cognate_min.cif", cog)   # high contact
    scr = tmp_path / "z_scramble.cif"; shutil.copy(FIX / "scramble_min.cif", scr)  # low contact
    _run(at.record_fold_result.handler(
        {"run_dir": rd, "clonotype_id": "z", "model_paths": [str(cog), str(scr)],
         "tool": "protenix"}))
    # a deliberately huge scalar threshold (999) would force suspect if it were used;
    # the derived scramble null (lower than the cognate contact) must win -> reliable.
    res = _run(at.qc_structure.handler(
        {"run_dir": rd, "clonotype_id": "z", "scramble_threshold": 999.0,
         "output_type": "structure", "tool": "protenix"}))
    assert res["structuredContent"]["qc_verdict"] == "reliable"
    stored = RunState(rd).read_stage("qc")
    assert stored[0]["tool"] == "protenix" and stored[0]["calibration_basis"] == "scramble_null"


def test_qc_protenix_cognate_losing_to_its_scramble_is_suspect(tmp_path, monkeypatch):
    # Beat-the-null gate: swap the roles so the "cognate" is the low-contact scramble
    # fixture and the "scramble" is the high-contact cognate fixture -> suspect.
    import shutil
    from rep2struct import qc
    monkeypatch.setattr(qc, "load_chains", lambda p: _clean_five_chains())
    rd = str(tmp_path / "run")
    cog = tmp_path / "z_cognate.cif"; shutil.copy(FIX / "scramble_min.cif", cog)   # low contact
    scr = tmp_path / "z_scramble.cif"; shutil.copy(FIX / "cognate_min.cif", scr)   # high contact
    _run(at.record_fold_result.handler(
        {"run_dir": rd, "clonotype_id": "z", "model_paths": [str(cog), str(scr)],
         "tool": "protenix"}))
    res = _run(at.qc_structure.handler(
        {"run_dir": rd, "clonotype_id": "z", "scramble_threshold": 0.0,
         "output_type": "structure", "tool": "protenix"}))
    assert res["structuredContent"]["qc_verdict"] == "suspect"


def test_build_fold_notebook_writes_wired_protenix_ipynb(tmp_path):
    import json
    from rep2struct.runstate import RunState
    from rep2struct.schema import FoldJob
    rd = str(tmp_path / "run")
    fasta = ">A\nAAAA\n>B\nBBBB\n>C\nCCCCMHC\n>D\nDDDD\n>E\nSIINFEKL\n"
    RunState(rd).write_stage("foldjobs", [FoldJob(clonotype_id="c1", construct_fasta=fasta)])
    res = _call(at.build_fold_notebook, {"run_dir": rd, "clonotype_id": "c1", "tool": "protenix"})
    nb = json.loads(Path(res["structuredContent"]["notebook_path"]).read_text())
    src = "".join(s for cell in nb["cells"] for s in cell["source"])
    assert "NotImplementedError" not in src              # protenix is wired
    assert "c1_cognate" in src and "c1_scramble" in src   # cognate + scramble pair, keyed by cid
    assert "protenix pred" in src and "--use_msa false" in src
    assert "SIINFEKL" in src                              # peptide embedded


def test_render_final_report(tmp_path):
    at.configure()
    rd = str(tmp_path / "run4")
    _call(at.ingest_repertoire, {"run_dir": rd, "csv_path": str(FIX / "tenx_tiny.csv")})
    _call(at.annotate_specificity, {"run_dir": rd})
    out = _call(at.render_final_report, {"run_dir": rd})
    report_path = Path(out["structuredContent"]["report_path"])
    assert report_path.exists()


def test_record_fold_result_keeps_tool_and_is_back_compatible(tmp_path):
    from rep2struct.runstate import RunState
    rd = str(tmp_path / "run")
    _run(at.record_fold_result.handler(
        {"run_dir": rd, "clonotype_id": "c1", "model_paths": ["c1.cif"], "tool": "tcrdock"}))
    done = RunState(rd).read_stage("folds")
    assert done["c1"]["tool"] == "tcrdock" and done["c1"]["paths"] == ["c1.cif"]


def test_qc_structure_binding_path_uses_binding_verdict(tmp_path):
    from rep2struct.runstate import RunState
    rd = str(tmp_path / "run")
    score_file = tmp_path / "c1.score"
    score_file.write_text("0.9")
    _run(at.record_fold_result.handler(
        {"run_dir": rd, "clonotype_id": "c1", "model_paths": [str(score_file)],
         "tool": "affinetune"}))
    res = _run(at.qc_structure.handler(
        {"run_dir": rd, "clonotype_id": "c1", "scramble_threshold": 0.5,
         "output_type": "binding_score", "tool": "affinetune"}))
    assert res["structuredContent"]["qc_verdict"] == "presented"
    stored = RunState(rd).read_stage("qc")
    assert stored[0]["tool"] == "affinetune"


def test_binding_threshold_is_set_from_the_groups_own_scramble_null(tmp_path):
    # The tcrdock/affinetune verdict must use the group's OWN scramble score (the sibling
    # {cid}_scramble.score the fold writes) as the threshold, NOT the caller-supplied number.
    # Validated tcrdock flu M1 null: cognate -11.219 beats scramble -20.574 -> presented.
    from rep2struct.runstate import RunState
    rd = str(tmp_path / "run")
    (tmp_path / "c1_cognate.score").write_text("-11.219")
    (tmp_path / "c1_scramble.score").write_text("-20.574")
    _run(at.record_fold_result.handler(
        {"run_dir": rd, "clonotype_id": "c1", "model_paths": [str(tmp_path / "c1_cognate.score")],
         "tool": "tcrdock"}))
    # a deliberately WRONG explicit threshold (999) would force not_presented if it were used;
    # the derived null (-20.574) must win and give presented.
    res = _run(at.qc_structure.handler(
        {"run_dir": rd, "clonotype_id": "c1", "scramble_threshold": 999.0,
         "output_type": "structure", "tool": "tcrdock"}))
    assert res["structuredContent"]["qc_verdict"] == "presented"


def test_binding_verdict_rejects_a_cognate_that_loses_to_its_own_scramble(tmp_path):
    # Beat-the-null gate: a cognate scoring BELOW its own scramble null is not_presented,
    # even if the caller passes a lax explicit threshold.
    rd = str(tmp_path / "run")
    (tmp_path / "c1_cognate.score").write_text("-25.0")   # worse (lower) than its scramble
    (tmp_path / "c1_scramble.score").write_text("-20.574")
    _run(at.record_fold_result.handler(
        {"run_dir": rd, "clonotype_id": "c1", "model_paths": [str(tmp_path / "c1_cognate.score")],
         "tool": "tcrdock"}))
    res = _run(at.qc_structure.handler(
        {"run_dir": rd, "clonotype_id": "c1", "scramble_threshold": -999.0,
         "output_type": "structure", "tool": "tcrdock"}))
    assert res["structuredContent"]["qc_verdict"] == "not_presented"


def test_binding_threshold_falls_back_to_explicit_when_no_scramble_sibling(tmp_path):
    # Back-compat: a bare {cid}.score with no sibling scramble uses the caller's threshold.
    rd = str(tmp_path / "run")
    (tmp_path / "c1.score").write_text("0.9")
    _run(at.record_fold_result.handler(
        {"run_dir": rd, "clonotype_id": "c1", "model_paths": [str(tmp_path / "c1.score")],
         "tool": "affinetune"}))
    res = _run(at.qc_structure.handler(
        {"run_dir": rd, "clonotype_id": "c1", "scramble_threshold": 0.5,
         "output_type": "binding_score", "tool": "affinetune"}))
    assert res["structuredContent"]["qc_verdict"] == "presented"


def test_output_type_is_derived_from_tool_not_agent_arg(tmp_path):
    # A binding tool recorded, but the caller passes the WRONG output_type="structure".
    # Derivation from the tool must still take the binding path.
    rd = str(tmp_path / "run")
    score_file = tmp_path / "c1.score"
    score_file.write_text("0.9")
    _run(at.record_fold_result.handler(
        {"run_dir": rd, "clonotype_id": "c1", "model_paths": [str(score_file)],
         "tool": "affinetune"}))
    res = _run(at.qc_structure.handler(
        {"run_dir": rd, "clonotype_id": "c1", "scramble_threshold": 0.5,
         "output_type": "structure", "tool": "affinetune"}))  # deliberately wrong output_type
    assert res["structuredContent"]["qc_verdict"] == "presented"


def test_output_type_for_helper_defaults_to_structure():
    from rep2struct import structure_tools as st
    assert st.output_type_for("affinetune") == "binding_score"
    assert st.output_type_for("protenix") == "structure"
    assert st.output_type_for("unknown-tool") == "structure"


def test_build_server():
    server = at.build_server()
    assert server["name"] == "rep2struct"


def test_list_structure_tools_returns_registry():
    res = _run(at.list_structure_tools.handler({"run_dir": "/tmp/whatever"}))
    names = {t["name"] for t in res["structuredContent"]["tools"]}
    assert names == {"protenix", "af3", "mhcfine", "tcrdock", "affinetune"}


def test_qc_structure_common_gate_fails_closed_on_missing_chains(tmp_path):
    # a mhcfine pose recorded but the CIF is a fixture missing chain C/E -> qc_failed
    rd = str(tmp_path / "run")
    bad = str(_ROOT / "tests" / "fixtures" / "threechain_min.cif")  # implementer: pick a fixture that lacks the expected mhcfine chains
    _run(at.record_fold_result.handler(
        {"run_dir": rd, "clonotype_id": "c1", "model_paths": [bad], "tool": "mhcfine"}))
    res = _run(at.qc_structure.handler(
        {"run_dir": rd, "clonotype_id": "c1", "scramble_threshold": 1.0,
         "output_type": "structure", "tool": "mhcfine"}))
    # mhcfine's qc_metric is peptide_groove, so a common-gate failure is an honest
    # pose_failed (never "structure (qc failed)") -- see follow-up B.
    assert res["structuredContent"]["qc_verdict"] == "pose_failed"


def test_qc_structure_dispatches_peptide_groove_for_mhcfine(tmp_path, monkeypatch):
    import numpy as np
    from rep2struct import qc
    from rep2struct.runstate import RunState
    rd = str(tmp_path / "run")
    cif = tmp_path / "c1.cif"; cif.write_text("dummy")
    _run(at.record_fold_result.handler(
        {"run_dir": rd, "clonotype_id": "c1", "model_paths": [str(cif)], "tool": "mhcfine"}))
    # synthetic pMHC: C=MHC heavy, D=b2m, E=peptide near C, all well separated (passes the gate)
    chains = {
        "C": np.array([[0, 0, 0], [20, 0, 0]], float),
        "D": np.array([[50, 0, 0]], float),
        "E": np.array([[1, 0, 0], [2, 0, 0]], float),
    }
    monkeypatch.setattr(qc, "load_chains", lambda p: chains)
    res = _run(at.qc_structure.handler(
        {"run_dir": rd, "clonotype_id": "c1", "scramble_threshold": 0.5,
         "output_type": "structure", "tool": "mhcfine"}))
    # the peptide_groove dispatch yields an honest pose-only verdict (never specificity)
    assert res["structuredContent"]["qc_verdict"] == "pose_only"
    stored = RunState(rd).read_stage("qc")
    assert stored[0]["tool"] == "mhcfine" and stored[0]["calibration_basis"] == "pose_quality"


def test_build_fold_notebook_writes_wired_mhcfine_ipynb(tmp_path):
    import json
    from rep2struct.runstate import RunState
    from rep2struct.schema import FoldJob
    rd = str(tmp_path / "run")
    fasta = ">A\nAAAA\n>B\nBBBB\n>C\nCCCCMHC\n>D\nDDDD\n>E\nSIINFEKL\n"
    RunState(rd).write_stage("foldjobs", [FoldJob(clonotype_id="c1", construct_fasta=fasta)])
    res = _call(at.build_fold_notebook, {"run_dir": rd, "clonotype_id": "c1", "tool": "mhcfine"})
    path = res["structuredContent"]["notebook_path"]
    assert Path(path).exists()
    nb = json.loads(Path(path).read_text())
    assert nb["nbformat"] == 4
    src = "".join(s for cell in nb["cells"] for s in cell["source"])
    assert "NotImplementedError" not in src          # mhcfine is wired
    assert "SIINFEKL" in src and "CCCCMHC" in src      # MHC heavy + peptide embedded
    assert "c1_cognate" in src and "c1_scramble" in src  # keys prefixed by clonotype id
    assert "np.string_" in src and "kalign" in src     # validated shim recipe carried through
    assert "numpy<2" not in src                          # stock numpy 2 kept (no ABI-poisoning downgrade)


def test_build_fold_notebook_unknown_job_is_reported(tmp_path):
    from rep2struct.runstate import RunState
    RunState(str(tmp_path / "run")).write_stage("foldjobs", [])
    res = _call(at.build_fold_notebook, {"run_dir": str(tmp_path / "run"),
                                         "clonotype_id": "nope", "tool": "mhcfine"})
    assert "no fold job" in res["content"][0]["text"].lower()
    assert "structuredContent" not in res


def test_prep_and_select_stamps_group_id(tmp_path):
    from rep2struct.runstate import RunState
    from rep2struct.schema import Clonotype, Annotation
    rd = str(tmp_path / "run")
    clon = Clonotype(id="c1", trav="TRAV1", cdr3a="CAA", trbv="TRBV2", cdr3b="CAB",
                     size=5, traj="TRAJ1", trbj="TRBJ1")
    ann = Annotation(clonotype_id="c1", annotatable=True, confidence_tier="high",
                     epitope="SIINFEKL", hla="A*02:01")
    RunState(rd).write_stage("ingest", [clon])
    RunState(rd).write_stage("annotate", [ann])
    at.configure(assign_fn=lambda c: c)  # no allele network call
    try:
        _run(at.prep_and_select.handler({"run_dir": rd, "top_n": 5}))
    finally:
        at.configure()
    jobs = RunState(rd).read_stage("foldjobs")
    assert jobs and all(j["group_id"] == "c1_tcr_human_structure" for j in jobs)
    assert jobs and all(j["msa_ref"] == "" for j in jobs)  # no runners injected -> MSA-free default; locks the build_msa stamping loop
    assert all(j["msa_basis"] == "none" for j in jobs)


def test_qc_persists_validity_summary_for_structure(tmp_path, monkeypatch):
    import numpy as np
    from rep2struct import qc
    from rep2struct.runstate import RunState
    rd = str(tmp_path / "run")
    cif = tmp_path / "c1.cif"; cif.write_text("x")
    _run(at.record_fold_result.handler({"run_dir": rd, "clonotype_id": "c1",
        "model_paths": [str(cif)], "tool": "mhcfine"}))
    chains = {"C": np.array([[0, 0, 0], [20, 0, 0]], float),
              "D": np.array([[50, 0, 0]], float),
              "E": np.array([[1, 0, 0], [2, 0, 0]], float)}
    monkeypatch.setattr(qc, "load_chains", lambda p: chains)
    _run(at.qc_structure.handler({"run_dir": rd, "clonotype_id": "c1",
        "scramble_threshold": 0.5, "output_type": "structure", "tool": "mhcfine"}))
    assert RunState(rd).read_stage("validity")["c1"] == "valid"


def test_qc_persists_validity_na_for_binding_score(tmp_path):
    from rep2struct.runstate import RunState
    rd = str(tmp_path / "run")
    score_path = tmp_path / "c1.score"; score_path.write_text("0.9")
    _run(at.record_fold_result.handler({"run_dir": rd, "clonotype_id": "c1",
        "model_paths": [str(score_path)], "tool": "affinetune"}))
    _run(at.qc_structure.handler({"run_dir": rd, "clonotype_id": "c1",
        "scramble_threshold": 0.5, "output_type": "binding_score", "tool": "affinetune"}))
    assert RunState(rd).read_stage("validity")["c1"] == "n/a (binding score)"
