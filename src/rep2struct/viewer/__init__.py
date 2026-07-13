"""Live viewer for the agent pipeline transcript.

`python -m rep2struct.viewer <run_dir> [--port N]` starts a tiny stdlib server
that serves a self contained timeline page and the run's transcript. The page
polls the transcript and appends new events live. Using a local server avoids the
file:// fetch restriction; no third party dependency is involved.
"""
from __future__ import annotations

import json
from functools import partial
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path

_HTML = Path(__file__).with_name("index.html")


def page_html() -> str:
    return _HTML.read_text()


def read_transcript(run_dir: str) -> list:
    """Parse the run's transcript into a list of event dicts, skipping any
    malformed line. Returns [] when the transcript does not exist yet."""
    path = Path(run_dir) / "transcript.jsonl"
    if not path.exists():
        return []
    rows = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            rows.append(json.loads(line))
        except json.JSONDecodeError:
            continue
    return rows


class _Handler(BaseHTTPRequestHandler):
    def __init__(self, *args, run_dir="", **kwargs):
        self._run_dir = run_dir
        super().__init__(*args, **kwargs)

    def log_message(self, *args):  # keep the console quiet during a demo
        pass

    def _send(self, body: bytes, content_type: str):
        self.send_response(200)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        if self.path.startswith("/transcript.jsonl"):
            body = json.dumps(read_transcript(self._run_dir)).encode()
            self._send(body, "application/json")
        else:
            self._send(page_html().encode(), "text/html; charset=utf-8")


def build_server(run_dir: str, port: int = 8000) -> HTTPServer:
    handler = partial(_Handler, run_dir=run_dir)
    return HTTPServer(("127.0.0.1", port), handler)


def serve(run_dir: str, port: int = 8000) -> None:
    httpd = build_server(run_dir, port)
    url = f"http://127.0.0.1:{httpd.server_address[1]}"
    print(f"R2S viewer on {url}  (transcript: {run_dir}/transcript.jsonl)")
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        httpd.server_close()
