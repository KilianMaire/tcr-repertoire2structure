import asyncio
from pathlib import Path

from rep2struct import agent_tools
from rep2struct.runstate import RunState


def _seed_job(rd):
    fasta = ">A\nAAAA\n>B\nBBBB\n>C\nCCCC\n>D\nDDDD\n>E\nGILGFVFTL\n"
    RunState(rd).write_stage("foldjobs", [{"clonotype_id": "c0",
                                           "construct_fasta": fasta,
                                           "group_id": "g0"}])


def _call(rd, route):
    return asyncio.run(agent_tools.build_fold_artifact.handler(
        {"run_dir": rd, "clonotype_id": "c0", "tool": "protenix",
         "compute_route": route}))


def test_colab_route_writes_a_notebook(tmp_path):
    rd = str(tmp_path); _seed_job(rd)
    r = _call(rd, "colab")
    sc = r["structuredContent"]
    assert sc["artifact_kind"] == "colab_notebook"
    assert sc["route_wired"] is True
    assert sc["artifact_path"].endswith("c0_protenix.ipynb")
    assert Path(sc["artifact_path"]).exists()


def test_local_gpu_route_writes_a_bash_script(tmp_path):
    rd = str(tmp_path); _seed_job(rd)
    r = _call(rd, "local_gpu")
    sc = r["structuredContent"]
    assert sc["artifact_kind"] == "bash_script"
    assert sc["route_wired"] is True
    assert sc["artifact_path"].endswith("c0_protenix.sh")
    body = Path(sc["artifact_path"]).read_text()
    assert body.startswith("#!/usr/bin/env bash")
    assert "protenix pred" in body


def test_unwired_ssh_route_still_scripts_but_flags_not_wired(tmp_path):
    rd = str(tmp_path); _seed_job(rd)
    r = _call(rd, "ssh")
    sc = r["structuredContent"]
    assert sc["artifact_kind"] == "bash_script"
    assert sc["route_wired"] is False
    assert Path(sc["artifact_path"]).exists()
