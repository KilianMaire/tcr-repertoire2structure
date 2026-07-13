import json
from pathlib import Path

from rep2struct.transcript import Transcript, MAX_TEXT


def _lines(run_dir):
    p = Path(run_dir) / "transcript.jsonl"
    return [json.loads(l) for l in p.read_text().splitlines() if l.strip()]


def test_record_writes_wellformed_jsonl_with_incrementing_seq(tmp_path):
    t = Transcript(str(tmp_path))
    t.record("user", text="which epitope?")
    t.record("agent_text", agent="orchestrator", text="delegating")
    t.record("tool_use", agent="orchestrator", tool="prep_and_select", text="top_n=1")
    rows = _lines(tmp_path)
    assert [r["seq"] for r in rows] == [0, 1, 2]
    assert rows[0] == {"seq": 0, "kind": "user", "agent": None, "tool": None,
                       "text": "which epitope?"}
    assert rows[2]["tool"] == "prep_and_select"


def test_text_is_truncated_at_bound(tmp_path):
    t = Transcript(str(tmp_path))
    t.record("tool_result", text="x" * (MAX_TEXT + 500))
    row = _lines(tmp_path)[0]
    assert len(row["text"]) <= MAX_TEXT
    assert row["text"].endswith("…") or len(row["text"]) == MAX_TEXT


def test_seq_continues_across_reopen(tmp_path):
    Transcript(str(tmp_path)).record("run_start")
    Transcript(str(tmp_path)).record("agent_text", text="second session")
    rows = _lines(tmp_path)
    assert [r["seq"] for r in rows] == [0, 1]


def test_ordering_preserved_across_many_records(tmp_path):
    t = Transcript(str(tmp_path))
    for i in range(50):
        t.record("agent_text", text=str(i))
    rows = _lines(tmp_path)
    assert [r["text"] for r in rows] == [str(i) for i in range(50)]
    assert [r["seq"] for r in rows] == list(range(50))


def test_io_error_is_swallowed_never_crashes_the_run(tmp_path):
    # Best effort: losing a transcript line must never propagate into the pipeline.
    t = Transcript(str(tmp_path))
    # Make the target unwritable by turning the file path into a directory.
    (Path(tmp_path) / "transcript.jsonl").mkdir()
    t.record("agent_text", text="should not raise")  # must not raise
