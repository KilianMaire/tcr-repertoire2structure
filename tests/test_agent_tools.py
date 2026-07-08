import asyncio
from pathlib import Path
from rep2struct import agent_tools as at

FIX = Path(__file__).parent / "fixtures"


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


def test_build_server():
    server = at.build_server()
    assert server["name"] == "rep2struct"


def test_list_structure_tools_returns_registry():
    res = _run(at.list_structure_tools.handler({"run_dir": "/tmp/whatever"}))
    names = {t["name"] for t in res["structuredContent"]["tools"]}
    assert names == {"protenix", "af3", "mhcfine", "tcrdock", "affinetune"}


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
