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
