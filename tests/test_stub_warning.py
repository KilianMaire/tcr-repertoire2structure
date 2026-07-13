"""The poly-G stub V-domain must be flagged to the user BEFORE folding, not only in QC."""
import asyncio

from rep2struct import agent_tools as at
from rep2struct.agent_tools import _stub_warning
from rep2struct.runstate import RunState
from rep2struct.schema import FoldJob


def _run(coro):
    return asyncio.run(coro)


def test_stub_warning_empty_when_all_real():
    assert _stub_warning([]) == ""
    assert _stub_warning([("c1", True), ("c2", True)]) == ""


def test_stub_warning_names_the_stub_clonotypes():
    w = _stub_warning([("c1", False), ("c2", True), ("c3", False)])
    assert "STUB" in w
    assert "2 of 3" in w
    assert "c1" in w and "c3" in w
    assert "c2" not in w


def _seed(tmp_path, jobs):
    RunState(str(tmp_path)).write_stage("foldjobs", jobs)
    return str(tmp_path)


def test_list_fold_jobs_surfaces_stub_before_folding(tmp_path):
    rd = _seed(tmp_path, [
        FoldJob(clonotype_id="cReal", construct_fasta="x", group_id="G", tool="protenix",
                tcr_reconstructed=True),
        FoldJob(clonotype_id="cStub", construct_fasta="x", group_id="G", tool="protenix",
                tcr_reconstructed=False),
    ])
    res = _run(at.list_fold_jobs.handler({"run_dir": rd}))
    text = res["content"][0]["text"]
    # the per-job line flags the stub, and the summary warns up front
    assert "cStub" in text and "stub_V=True" in text
    assert "cReal" in text and "stub_V=False" in text
    assert "STUB V-domain" in text or "STUB" in text


def test_list_fold_jobs_no_warning_when_all_real(tmp_path):
    rd = _seed(tmp_path, [
        FoldJob(clonotype_id="cReal", construct_fasta="x", group_id="G", tool="protenix",
                tcr_reconstructed=True),
    ])
    res = _run(at.list_fold_jobs.handler({"run_dir": rd}))
    assert "STUB" not in res["content"][0]["text"]
