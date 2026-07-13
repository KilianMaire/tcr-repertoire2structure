"""Best effort transcript of the agent pipeline.

Every meaningful message the run produces (user turns, agent text, tool calls and
results, sub agent handoffs) is appended as one JSON line to
`<run_dir>/transcript.jsonl`. The live viewer tails that file. Writing is best
effort: a failure to record a line must never propagate into the pipeline, which
is the actual product.
"""
from __future__ import annotations

import json
from pathlib import Path

MAX_TEXT = 2000


def _truncate(text: str) -> str:
    text = "" if text is None else str(text)
    if len(text) <= MAX_TEXT:
        return text
    return text[: MAX_TEXT - 1] + "…"


class Transcript:
    def __init__(self, run_dir: str) -> None:
        self.path = Path(run_dir) / "transcript.jsonl"
        try:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            self._seq = self._existing_line_count()
        except Exception:
            self._seq = 0

    def _existing_line_count(self) -> int:
        if not self.path.exists():
            return 0
        with open(self.path) as fh:
            return sum(1 for line in fh if line.strip())

    def record(self, kind: str, agent: str | None = None,
               tool: str | None = None, text: str = "") -> None:
        try:
            line = {"seq": self._seq, "kind": kind, "agent": agent,
                    "tool": tool, "text": _truncate(text)}
            with open(self.path, "a") as fh:
                fh.write(json.dumps(line) + "\n")
            self._seq += 1
        except Exception:
            # Best effort: never let observability break a run.
            pass
