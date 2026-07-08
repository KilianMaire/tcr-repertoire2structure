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


def test_build_server():
    server = at.build_server()
    assert server["name"] == "rep2struct"
