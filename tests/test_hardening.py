"""Regression tests for the iteration-1 hardening pass: routing persistence, class II
refusal, the poly-G stub honesty guard, and the corrupt-score-file guard."""
import asyncio

from rep2struct import agent_tools as at
from rep2struct.foldprep import mhc_class_of, build_construct
from rep2struct.runstate import RunState
from rep2struct.schema import FoldJob, Annotation, Clonotype


def _run(coro):
    return asyncio.run(coro)


# --- Bug 1: the strategist's choice is persisted onto the jobs ---------------

def test_assign_group_tool_persists_tool_onto_group_jobs(tmp_path):
    rd = str(tmp_path)
    rs = RunState(rd)
    rs.write_stage("foldjobs", [
        FoldJob(clonotype_id="c1", construct_fasta="x", group_id="G1"),
        FoldJob(clonotype_id="c2", construct_fasta="x", group_id="G1"),
        FoldJob(clonotype_id="c3", construct_fasta="x", group_id="G2"),
    ])
    r = _run(at.assign_group_tool.handler({"run_dir": rd, "group_id": "G1", "tool": "tcrdock"}))
    assert r["structuredContent"]["assigned"] == 2
    jobs = {j["clonotype_id"]: j for j in rs.read_stage("foldjobs")}
    assert jobs["c1"]["tool"] == "tcrdock" and jobs["c2"]["tool"] == "tcrdock"
    assert jobs["c3"]["tool"] is None  # other group untouched
    # and list_fold_jobs now surfaces it in the text the executor filters on
    listed = _run(at.list_fold_jobs.handler({"run_dir": rd}))
    assert "tool=tcrdock" in listed["content"][0]["text"]


# --- scale guards: annotate caps by size, report shows only selected ------------

def test_annotate_caps_by_size_and_reports_skipped(tmp_path, monkeypatch):
    monkeypatch.setenv("R2S_ANNOTATE_CAP", "2")
    at.configure(sim_fn=lambda *a, **k: ([], "tcrdist", 0, []))
    rd = str(tmp_path)
    # size-sorted, as parse_10x would leave them
    RunState(rd).write_stage("ingest", [
        Clonotype(id=f"c{i}", trav="TRAV1", cdr3a="CAAA", trbv="TRBV1", cdr3b="CBBB", size=s)
        for i, s in enumerate([50, 40, 30, 20])])
    r = _run(at.annotate_specificity.handler({"run_dir": rd}))
    sc = r["structuredContent"]
    assert sc["annotated"] == 2 and sc["skipped"] == 2   # only the two largest annotated
    assert "not annotated" in r["content"][0]["text"]
    at.configure()


def test_report_shows_only_selected_clonotypes(tmp_path):
    rd = str(tmp_path)
    at.configure()
    RunState(rd).write_stage("ingest", [
        Clonotype(id=f"c{i}", trav="TRAV1", cdr3a="CAAA", trbv="TRBV1", cdr3b="CBBB", size=9)
        for i in range(3)])
    RunState(rd).write_stage("annotate", [
        Annotation(clonotype_id=f"c{i}", annotatable=False, confidence_tier="unannotatable")
        for i in range(3)])
    RunState(rd).write_stage("foldjobs", [
        FoldJob(clonotype_id="c0", construct_fasta="x", group_id="G")])
    r = _run(at.render_final_report.handler({"run_dir": rd}))
    html = (tmp_path / "report.html").read_text()
    assert "c0" in html and "c1" not in html and "c2" not in html  # only the folded one


# --- prep_and_select is an immutable checkpoint: resume must not clobber it ------

def test_prep_and_select_does_not_clobber_existing_foldjobs(tmp_path):
    rd = str(tmp_path)
    rs = RunState(rd)
    # a resume-state foldjobs stage carrying the strategist's persisted tool + MSA
    rs.write_stage("foldjobs", [
        FoldJob(clonotype_id="c1", construct_fasta="x", group_id="G1", tool="tcrdock",
                msa_basis="colab_cpu:2/5")])
    r = _run(at.prep_and_select.handler({"run_dir": rd, "top_n": 50}))
    assert r["structuredContent"].get("reused") is True
    jobs = rs.read_stage("foldjobs")
    assert len(jobs) == 1 and jobs[0]["tool"] == "tcrdock"        # tag preserved
    assert jobs[0]["msa_basis"] == "colab_cpu:2/5"                # MSA preserved


# --- Class II is derived from the HLA and refused, not silently mis-modelled ---

def test_mhc_class_of_hla():
    assert mhc_class_of("HLA-A*02:01") == 1
    assert mhc_class_of("HLA-B*07:02") == 1
    assert mhc_class_of("HLA-DRB1*15:01") == 2
    assert mhc_class_of("HLA-DQB1*06:02") == 2
    assert mhc_class_of(None) == 1


def test_build_construct_refuses_class_ii():
    c = Clonotype(id="c1", trav="TRAV1", cdr3a="CAAA", trbv="TRBV1", cdr3b="CBBB", size=3)
    a2 = Annotation(clonotype_id="c1", annotatable=True, confidence_tier="high",
                    epitope="PKYVKQNTLKLAT", hla="HLA-DRB1*15:01")
    seqs = {"c1": {"A": "AAAA", "B": "BBBB", "reconstructed": True}}
    mhc = {"HLA-DRB1*15:01": {"heavy": "H", "b2m": "M"}}
    assert build_construct(c, a2, seqs, mhc) is None  # class II -> no class I construct
    a1 = Annotation(clonotype_id="c1", annotatable=True, confidence_tier="high",
                    epitope="GILGFVFTL", hla="HLA-A*02:01")
    job = build_construct(c, a1, {"c1": {"A": "AAAA", "B": "BBBB", "reconstructed": True}},
                          {"HLA-A*02:01": {"heavy": "H", "b2m": "M"}})
    assert job is not None and job.mhc_class == 1 and job.tcr_reconstructed is True


# --- Bug 2: a poly-G stub V-domain can never be graded reliable/presented -----

def _seed_score_fold(tmp_path, cid, reconstructed, cognate="-11.0", scramble="-20.0"):
    rd = str(tmp_path)
    RunState(rd).write_stage("foldjobs", [
        FoldJob(clonotype_id=cid, construct_fasta="x", group_id="G", tool="tcrdock",
                tcr_reconstructed=reconstructed)])
    out = tmp_path / "out"; out.mkdir(exist_ok=True)
    (out / f"{cid}_cognate.score").write_text(cognate + "\n")
    (out / f"{cid}_scramble.score").write_text(scramble + "\n")
    _run(at.record_local_folds.handler({"run_dir": rd, "tool": "tcrdock"}))
    return rd


def test_stub_v_domain_is_downgraded_from_presented(tmp_path):
    rd = _seed_score_fold(tmp_path, "cS", reconstructed=False)
    res = _run(at.qc_structure.handler(
        {"run_dir": rd, "clonotype_id": "cS", "scramble_threshold": 0.0,
         "output_type": "structure", "tool": "tcrdock"}))
    assert res["structuredContent"]["qc_verdict"] == "suspect"
    assert "stub" in res["structuredContent"]["reason"]


def test_real_v_domain_still_presented(tmp_path):
    rd = _seed_score_fold(tmp_path, "cR", reconstructed=True)
    res = _run(at.qc_structure.handler(
        {"run_dir": rd, "clonotype_id": "cR", "scramble_threshold": 0.0,
         "output_type": "structure", "tool": "tcrdock"}))
    assert res["structuredContent"]["qc_verdict"] == "presented"


# --- corrupt .score fails one clonotype, does not crash the QC tool -----------

def test_corrupt_score_file_yields_qc_failed(tmp_path):
    rd = str(tmp_path)
    RunState(rd).write_stage("foldjobs", [
        FoldJob(clonotype_id="cB", construct_fasta="x", group_id="G", tool="tcrdock")])
    out = tmp_path / "out"; out.mkdir()
    (out / "cB_cognate.score").write_text("not_a_number\n")
    (out / "cB_scramble.score").write_text("-20.0\n")
    _run(at.record_local_folds.handler({"run_dir": rd, "tool": "tcrdock"}))
    res = _run(at.qc_structure.handler(
        {"run_dir": rd, "clonotype_id": "cB", "scramble_threshold": 0.0,
         "output_type": "structure", "tool": "tcrdock"}))
    assert res["structuredContent"]["qc_verdict"] == "qc_failed"
    assert "unreadable" in res["structuredContent"]["reason"]
