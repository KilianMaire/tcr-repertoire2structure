import json
from rep2struct.runstate import RunState

def test_write_read_roundtrip(tmp_path):
    rs = RunState(tmp_path / "run1")
    assert not rs.stage_done("ingest")
    rs.write_stage("ingest", {"clonotypes": [{"id": "c1", "size": 3}]})
    assert rs.stage_done("ingest")
    got = rs.read_stage("ingest")
    assert got["clonotypes"][0]["id"] == "c1"

def test_path_for_is_under_run_dir(tmp_path):
    rs = RunState(tmp_path / "run2")
    p = rs.path_for("fold")
    assert str(p).startswith(str(tmp_path / "run2"))
