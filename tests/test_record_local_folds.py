# tests/test_record_local_folds.py
import asyncio
from pathlib import Path

from rep2struct import agent_tools
from rep2struct.runstate import RunState


def _write_cif(p):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("data_x\n_atom_site.group_PDB\n")


def test_scan_finds_cognate_and_scramble_by_directory(tmp_path):
    out = tmp_path / "out"
    _write_cif(out / "c0_cognate" / "preds" / "cognate_sample_0.cif")
    _write_cif(out / "c0_scramble" / "preds" / "scramble_sample_0.cif")
    _write_cif(out / "c1_cognate" / "preds" / "cognate_sample_0.cif")
    found = agent_tools.scan_recorded_folds(str(tmp_path))
    assert set(found) == {"c0", "c1"}
    c0 = [p for p in found["c0"]["paths"]]
    assert any("c0_cognate" in p for p in c0)
    assert any("c0_scramble" in p for p in c0)
    assert found["c0"]["tool"] == "protenix"


def test_record_local_folds_writes_the_folds_stage(tmp_path):
    out = tmp_path / "out"
    _write_cif(out / "c0_cognate" / "preds" / "s0.cif")
    r = asyncio.run(agent_tools.record_local_folds.handler(
        {"run_dir": str(tmp_path), "tool": "protenix"}))
    assert r["structuredContent"]["recorded"] == 1
    done = RunState(str(tmp_path)).read_stage("folds")
    assert "c0" in done and done["c0"]["tool"] == "protenix"


def test_scan_finds_binding_score_folds_by_filename(tmp_path):
    # tcrdock/affinetune write sibling {cid}_cognate.score / {cid}_scramble.score files in
    # one folder, not _cognate/_scramble dirs; the scan must key on the filename.
    out = tmp_path / "out"
    out.mkdir(parents=True)
    (out / "c0_cognate.score").write_text("-11.219\n")
    (out / "c0_scramble.score").write_text("-20.574\n")
    (out / "c1_cognate.score").write_text("-9.5\n")
    found = agent_tools.scan_recorded_folds(str(tmp_path), tool="tcrdock")
    assert set(found) == {"c0", "c1"}
    assert found["c0"]["tool"] == "tcrdock"
    # only the cognate path is recorded; QC's _scramble_null finds the scramble sibling
    assert found["c0"]["paths"] == [str(out / "c0_cognate.score")]


def test_record_local_folds_records_score_tool_from_disk(tmp_path):
    out = tmp_path / "out"
    out.mkdir(parents=True)
    (out / "c0_cognate.score").write_text("-11.219\n")
    (out / "c0_scramble.score").write_text("-20.574\n")
    r = asyncio.run(agent_tools.record_local_folds.handler(
        {"run_dir": str(tmp_path), "tool": "tcrdock"}))
    assert r["structuredContent"]["recorded"] == 1
    done = RunState(str(tmp_path)).read_stage("folds")
    assert done["c0"]["tool"] == "tcrdock"
    assert done["c0"]["paths"][0].endswith("c0_cognate.score")


def test_record_local_folds_does_not_clobber_existing(tmp_path):
    rs = RunState(str(tmp_path))
    rs.write_stage("folds", {"c0": {"paths": ["cloud/complete.cif"], "tool": "protenix"}})
    out = tmp_path / "out"
    _write_cif(out / "c0_cognate" / "preds" / "partial.cif")  # disk scan would find only this
    r = asyncio.run(agent_tools.record_local_folds.handler(
        {"run_dir": str(tmp_path), "tool": "protenix"}))
    assert r["structuredContent"]["recorded"] == 0
    done = rs.read_stage("folds")
    assert done["c0"]["paths"] == ["cloud/complete.cif"]  # prior complete record preserved
