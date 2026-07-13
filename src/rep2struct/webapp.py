"""Interactive web front end for the R2S pipeline.

`python -m rep2struct.webapp [BASE_DIR] [--port N]` serves a chat landing page:
drop a 10x CSV to start a run, answer the intake agent in the chat, and watch the
multi agent orchestration stream into the timeline. Past runs stay in the sidebar,
each backed by its own run_dir on disk.

Kept deliberately in the standard library (threaded http.server plus one asyncio
loop per run in a worker thread, bridged by a queue). One run is active at a time,
which is all a demo needs; browsing a past run is read only.
"""
from __future__ import annotations

import asyncio
import json
import threading
from functools import partial
from http.server import ThreadingHTTPServer, BaseHTTPRequestHandler
from pathlib import Path
from urllib.parse import urlparse, parse_qs

from .viewer import read_transcript

_HTML = Path(__file__).parent / "viewer" / "app.html"


def app_html() -> str:
    return _HTML.read_text()


def next_run_dir(base_dir: str) -> str:
    """Sequential run_N under base_dir, so runs sort and names are deterministic."""
    base = Path(base_dir)
    base.mkdir(parents=True, exist_ok=True)
    n = sum(1 for p in base.glob("run_*") if p.is_dir())
    return str(base / f"run_{n}")


def save_upload(run_dir: str, filename: str, data: bytes) -> str:
    Path(run_dir).mkdir(parents=True, exist_ok=True)
    name = Path(filename).name or "input.csv"
    if not name.lower().endswith(".csv"):
        name = name + ".csv"
    out = Path(run_dir) / name
    out.write_bytes(data)
    return str(out)


def list_runs(base_dir: str) -> list:
    base = Path(base_dir)
    if not base.exists():
        return []
    runs = []
    for p in sorted(base.glob("run_*")):
        if p.is_dir():
            runs.append({"name": p.name, "events": len(read_transcript(str(p)))})
    return runs


class _Session:
    def __init__(self, run_dir, loop, queue, thread):
        self.run_dir = run_dir
        self.loop = loop
        self.queue = queue
        self.thread = thread
        self.awaiting = False  # True while the run is blocked waiting for a user answer


class AppState:
    """Owns the base dir and the single active run. The runner is injected so the
    HTTP bridge can be tested without the real agent client."""

    def __init__(self, base_dir: str, runner=None):
        self.base_dir = base_dir
        self._runner = runner or _default_runner()
        self._session: _Session | None = None

    def start(self, csv_bytes: bytes, filename: str) -> str:
        run_dir = next_run_dir(self.base_dir)
        csv_path = save_upload(run_dir, filename, csv_bytes)
        loop = asyncio.new_event_loop()
        queue: asyncio.Queue = asyncio.Queue()

        async def source():
            # Mark the run as waiting for the user while blocked on the queue, so
            # the UI can tell "agent working" from "your turn" honestly.
            if self._session:
                self._session.awaiting = True
            try:
                return await queue.get()
            finally:
                if self._session:
                    self._session.awaiting = False

        def worker():
            asyncio.set_event_loop(loop)
            try:
                loop.run_until_complete(
                    self._runner(run_dir, answer_source=source, data_path=csv_path))
            finally:
                loop.close()

        thread = threading.Thread(target=worker, daemon=True)
        self._session = _Session(run_dir, loop, queue, thread)
        thread.start()
        return run_dir

    def answer(self, text: str) -> bool:
        """Enqueue a user answer onto the active run's loop. Returns False if no
        run is active."""
        s = self._session
        if not s:
            return False
        s.loop.call_soon_threadsafe(s.queue.put_nowait, text)
        return True

    def is_running(self) -> bool:
        return bool(self._session and self._session.thread.is_alive())

    def is_awaiting(self) -> bool:
        """True when the active run is blocked waiting for the user's next answer."""
        return bool(self.is_running() and self._session.awaiting)

    def current_run_name(self) -> str | None:
        return Path(self._session.run_dir).name if self._session else None


def _default_runner():
    # Imported lazily so the module loads without the SDK during tests.
    from .cli import run_session
    return run_session


class _Handler(BaseHTTPRequestHandler):
    def __init__(self, *args, state=None, **kwargs):
        self._state = state
        super().__init__(*args, **kwargs)

    def log_message(self, *args):
        pass

    def _send(self, body: bytes, ctype: str, code: int = 200):
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _json(self, obj, code: int = 200):
        self._send(json.dumps(obj).encode(), "application/json", code)

    def do_GET(self):
        parsed = urlparse(self.path)
        if parsed.path == "/runs":
            self._json(list_runs(self._state.base_dir))
        elif parsed.path == "/transcript.jsonl":
            run = parse_qs(parsed.query).get("run", [None])[0]
            run_dir = (str(Path(self._state.base_dir) / run) if run
                       else (self._state._session.run_dir if self._state._session else ""))
            self._json(read_transcript(run_dir) if run_dir else [])
        elif parsed.path == "/status":
            self._json({"running": self._state.is_running(),
                        "awaiting": self._state.is_awaiting(),
                        "current": self._state.current_run_name()})
        else:
            self._send(app_html().encode(), "text/html; charset=utf-8")

    def do_POST(self):
        length = int(self.headers.get("Content-Length", 0))
        body = self.rfile.read(length) if length else b""
        parsed = urlparse(self.path)
        if parsed.path == "/start":
            filename = self.headers.get("X-Filename", "input.csv")
            run_dir = self._state.start(body, filename)
            self._json({"run": Path(run_dir).name})
        elif parsed.path == "/answer":
            ok = self._state.answer(body.decode("utf-8"))
            self._json({"ok": ok}, 200 if ok else 409)
        else:
            self._json({"error": "not found"}, 404)


def build_app_server(base_dir: str, port: int = 8000, runner=None) -> ThreadingHTTPServer:
    state = AppState(base_dir, runner=runner)
    handler = partial(_Handler, state=state)
    return ThreadingHTTPServer(("127.0.0.1", port), handler)


def serve(base_dir: str = "runs", port: int = 8000, open_browser: bool = False) -> None:
    import webbrowser

    httpd = build_app_server(base_dir, port)
    url = f"http://127.0.0.1:{httpd.server_address[1]}"
    print(f"R2S web app on {url}  (runs under {base_dir}/)")
    if open_browser:
        # Fire once the server is accepting; a short timer beats racing serve_forever.
        threading.Timer(0.5, lambda: webbrowser.open(url)).start()
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        httpd.server_close()


def parse_args(args):
    """Return (base_dir, port, open_browser) from the argv tail. Pure and tested."""
    args = list(args)
    open_browser = True
    if "--no-browser" in args:
        open_browser = False
        args.remove("--no-browser")
    port = 8000
    if "--port" in args:
        i = args.index("--port")
        port = int(args[i + 1])
        del args[i:i + 2]
    base_dir = args[0] if args else "runs"
    return base_dir, port, open_browser


def main():
    import sys
    args = sys.argv[1:]
    if args and args[0] in ("-h", "--help"):
        print("usage: r2s [BASE_DIR] [--port N] [--no-browser]\n\n"
              "Starts the R2S web app and opens it in your browser. Drop a 10x CSV\n"
              "to start a run, answer the intake agent in the chat, watch the\n"
              "orchestration stream. Runs are stored under BASE_DIR (default runs/).\n"
              "  --port N       port to bind (default 8000)\n"
              "  --no-browser   do not open a browser (print the URL only)")
        return
    base_dir, port, open_browser = parse_args(args)
    serve(base_dir, port=port, open_browser=open_browser)


if __name__ == "__main__":
    main()
