import json
import threading
import urllib.request
from pathlib import Path

from rep2struct import viewer
from rep2struct.transcript import Transcript


def test_read_transcript_empty_when_absent(tmp_path):
    assert viewer.read_transcript(str(tmp_path)) == []


def test_read_transcript_parses_rows_and_skips_malformed(tmp_path):
    t = Transcript(str(tmp_path))
    t.record("agent_text", agent="orchestrator", text="hello")
    # append a malformed line by hand; the reader must skip it, not crash
    with open(Path(tmp_path) / "transcript.jsonl", "a") as fh:
        fh.write("{not json}\n")
    rows = viewer.read_transcript(str(tmp_path))
    assert len(rows) == 1
    assert rows[0]["text"] == "hello"


def test_page_html_is_self_contained(tmp_path):
    html = viewer.page_html()
    assert "<html" in html.lower()
    # no external requests allowed (offline, CSP-friendly demo surface)
    assert "http://" not in html and "https://" not in html
    assert "/transcript.jsonl" in html  # it polls the transcript endpoint


def test_server_serves_page_and_transcript(tmp_path):
    t = Transcript(str(tmp_path))
    t.record("user", text="which epitope?")
    httpd = viewer.build_server(str(tmp_path), port=0)
    threading.Thread(target=httpd.handle_request, daemon=True).start()
    port = httpd.server_address[1]
    body = urllib.request.urlopen(f"http://127.0.0.1:{port}/transcript.jsonl", timeout=5).read()
    rows = json.loads(body)
    assert rows[0]["text"] == "which epitope?"
    httpd.server_close()
