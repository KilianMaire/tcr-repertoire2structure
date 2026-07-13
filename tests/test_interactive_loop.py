import asyncio
import json
from pathlib import Path

from rep2struct import cli
from rep2struct.transcript import Transcript


class _FakeClient:
    """Minimal stand in: receive_response yields nothing, query records the text."""
    def __init__(self):
        self.queries = []

    async def receive_response(self):
        return
        yield  # make it an async generator

    async def query(self, text):
        self.queries.append(text)


def _answers_from(seq):
    it = iter(seq)

    async def source():
        try:
            return next(it)
        except StopIteration:
            return None
    return source


def _user_lines(run_dir):
    p = Path(run_dir) / "transcript.jsonl"
    rows = [json.loads(l) for l in p.read_text().splitlines() if l.strip()]
    return [r["text"] for r in rows if r["kind"] == "user"]


def test_interactive_loop_feeds_answers_until_source_ends(tmp_path):
    client = _FakeClient()
    t = Transcript(str(tmp_path))
    asyncio.run(cli._interactive_loop(client, t, _answers_from(["human", "A*02:01", "colab"])))
    assert client.queries == ["human", "A*02:01", "colab"]
    assert _user_lines(tmp_path) == ["human", "A*02:01", "colab"]


def test_interactive_loop_stops_on_blank_answer(tmp_path):
    client = _FakeClient()
    t = Transcript(str(tmp_path))
    asyncio.run(cli._interactive_loop(client, t, _answers_from(["human", "   ", "never"])))
    assert client.queries == ["human"]


def test_interactive_loop_stops_on_none(tmp_path):
    client = _FakeClient()
    t = Transcript(str(tmp_path))
    asyncio.run(cli._interactive_loop(client, t, _answers_from([])))
    assert client.queries == []
