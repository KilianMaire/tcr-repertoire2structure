import asyncio
from pathlib import Path

from rep2struct import agent_tools, intake
from rep2struct.intake import IntakeSpec
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


def test_local_gpu_route_uses_working_path_from_intake(tmp_path):
    rd = str(tmp_path); _seed_job(rd)
    intake.save_intake(rd, IntakeSpec("10x_vdj", "/data/c.csv", "which epitope?",
                                      "local_gpu", {"working_path": "/scratch/h100run"}))
    r = _call(rd, "local_gpu")
    sc = r["structuredContent"]
    body = Path(sc["artifact_path"]).read_text()
    assert 'cd "/scratch/h100run"' in body


def test_unwired_ssh_route_still_scripts_but_flags_not_wired(tmp_path):
    rd = str(tmp_path); _seed_job(rd)
    r = _call(rd, "ssh")
    sc = r["structuredContent"]
    assert sc["artifact_kind"] == "bash_script"
    assert sc["route_wired"] is False
    assert Path(sc["artifact_path"]).exists()


def _fasta(pep):
    return f">A\nAAAA\n>B\nBBBB\n>C\nCCCC\n>D\nDDDD\n>E\n{pep}\n"


def test_group_artifact_batches_all_clonotypes_into_one_notebook(tmp_path):
    rd = str(tmp_path)
    RunState(rd).write_stage("foldjobs", [
        {"clonotype_id": "c0", "construct_fasta": _fasta("GILGFVFTL"), "group_id": "g0"},
        {"clonotype_id": "c1", "construct_fasta": _fasta("NLVPMVATV"), "group_id": "g0"},
        {"clonotype_id": "c2", "construct_fasta": _fasta("ELAGIGILTV"), "group_id": "g1"},
    ])
    r = asyncio.run(agent_tools.build_group_artifact.handler(
        {"run_dir": rd, "group_id": "g0", "tool": "protenix", "compute_route": "colab"}))
    sc = r["structuredContent"]
    assert sc["artifact_kind"] == "colab_notebook"
    assert sc["n_artifacts"] == 1 and sc["n_clonotypes"] == 2
    art = sc["artifacts"][0]
    assert sorted(art["clonotypes"]) == ["c0", "c1"]          # both group-0 clonotypes, not c2
    assert art["artifact_path"].endswith("g0_protenix.ipynb")  # named by group, one file
    body = Path(art["artifact_path"]).read_text()
    assert "c0_cognate" in body and "c1_cognate" in body      # merged into one notebook
    assert "c2_cognate" not in body                            # other group excluded


def test_group_artifact_skips_already_folded_clonotypes(tmp_path):
    rd = str(tmp_path)
    RunState(rd).write_stage("foldjobs", [
        {"clonotype_id": "c0", "construct_fasta": _fasta("GILGFVFTL"), "group_id": "g0"},
        {"clonotype_id": "c1", "construct_fasta": _fasta("NLVPMVATV"), "group_id": "g0"},
    ])
    RunState(rd).write_stage("folds", {"c0": {"paths": ["x.cif"], "tool": "protenix"}})
    r = asyncio.run(agent_tools.build_group_artifact.handler(
        {"run_dir": rd, "group_id": "g0", "tool": "protenix", "compute_route": "colab"}))
    sc = r["structuredContent"]
    assert sc["n_clonotypes"] == 1
    assert sc["artifacts"][0]["clonotypes"] == ["c1"]  # c0 already done, only c1 pending


def test_group_artifact_shards_beyond_max_batch(tmp_path):
    rd = str(tmp_path)
    RunState(rd).write_stage("foldjobs", [
        {"clonotype_id": f"c{i}", "construct_fasta": _fasta("GILGFVFTL"), "group_id": "g0"}
        for i in range(40)])
    r = asyncio.run(agent_tools.build_group_artifact.handler(
        {"run_dir": rd, "group_id": "g0", "tool": "protenix", "compute_route": "colab"}))
    sc = r["structuredContent"]
    assert sc["n_clonotypes"] == 40 and sc["n_artifacts"] == 3   # 16 + 16 + 8
    assert all(len(a["clonotypes"]) <= 16 for a in sc["artifacts"])
    assert sc["artifacts"][0]["artifact_path"].endswith("g0_protenix_part1.ipynb")


def test_group_bash_route_fails_loud_for_non_protenix_tool(tmp_path):
    rd = str(tmp_path)
    RunState(rd).write_stage("foldjobs", [
        {"clonotype_id": "c0", "construct_fasta": _fasta("GILGFVFTL"), "group_id": "g0",
         "tool": "tcrdock"}])
    r = asyncio.run(agent_tools.build_group_artifact.handler(
        {"run_dir": rd, "group_id": "g0", "tool": "tcrdock", "compute_route": "local_gpu"}))
    sc = r["structuredContent"]
    # local_gpu is a "wired" route, but only Protenix has a bash runner: must NOT claim wired
    assert sc["route_wired"] is False
    body = Path(sc["artifacts"][0]["artifact_path"]).read_text()
    assert "no bash runner wired for tcrdock" in body and "exit 1" in body


def test_group_artifact_rejects_tool_mismatch_with_persisted(tmp_path):
    rd = str(tmp_path)
    RunState(rd).write_stage("foldjobs", [
        {"clonotype_id": "c0", "construct_fasta": _fasta("GILGFVFTL"), "group_id": "g0",
         "tool": "tcrdock"}])  # strategist persisted tcrdock
    r = asyncio.run(agent_tools.build_group_artifact.handler(
        {"run_dir": rd, "group_id": "g0", "tool": "protenix", "compute_route": "colab"}))
    assert "mismatch" in r["content"][0]["text"]  # executor passed protenix, refused
