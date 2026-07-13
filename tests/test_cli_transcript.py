import json
from pathlib import Path

from claude_agent_sdk import TextBlock, ToolUseBlock, ToolResultBlock

from rep2struct import cli
from rep2struct.transcript import Transcript


def _lines(run_dir):
    p = Path(run_dir) / "transcript.jsonl"
    return [json.loads(l) for l in p.read_text().splitlines() if l.strip()]


def test_record_blocks_captures_text_tooluse_and_toolresult(tmp_path, capsys):
    t = Transcript(str(tmp_path))
    blocks = [
        TextBlock(text="routing this group to tcrdock"),
        ToolUseBlock(id="1", name="build_fold_notebook",
                     input={"tool": "tcrdock", "clonotype_id": "7a5e14a1017e"}),
        ToolResultBlock(tool_use_id="1",
                        content="wrote notebooks/7a5e14a1017e_tcrdock.ipynb",
                        is_error=False),
    ]
    text = cli.record_blocks(blocks, t, agent="structure-strategist", is_assistant=True)

    rows = _lines(tmp_path)
    assert [r["kind"] for r in rows] == ["agent_text", "tool_use", "tool_result"]
    assert rows[0]["agent"] == "structure-strategist"
    assert rows[1]["tool"] == "build_fold_notebook"
    assert "tcrdock" in rows[1]["text"]
    assert "7a5e14a1017e_tcrdock.ipynb" in rows[2]["text"]
    # assistant text is still returned (and printed) for the terminal path
    assert text == ["routing this group to tcrdock"]
    assert "routing this group" in capsys.readouterr().out


def test_record_blocks_non_assistant_records_but_returns_no_text(tmp_path):
    t = Transcript(str(tmp_path))
    blocks = [ToolResultBlock(tool_use_id="9", content="done", is_error=False)]
    text = cli.record_blocks(blocks, t, agent="orchestrator", is_assistant=False)
    assert text == []
    assert _lines(tmp_path)[0]["kind"] == "tool_result"


def test_record_blocks_tolerates_no_transcript(tmp_path):
    # A run with no transcript sink must still stream text and not crash.
    blocks = [TextBlock(text="hi")]
    assert cli.record_blocks(blocks, None, agent="orchestrator", is_assistant=True) == ["hi"]
