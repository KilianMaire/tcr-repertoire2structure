import json
import time
from pathlib import Path

from rep2struct import webapp
from rep2struct.transcript import Transcript


def _read_users(run_dir):
    p = Path(run_dir) / "transcript.jsonl"
    if not p.exists():
        return []
    rows = [json.loads(l) for l in p.read_text().splitlines() if l.strip()]
    return [r["text"] for r in rows if r["kind"] == "user"]


async def _fake_runner(run_dir, answer_source=None, data_path=None):
    """Stands in for cli.run_session: echoes each answer into the transcript until
    a blank/None answer ends it. Records data_path so the test can assert wiring."""
    t = Transcript(run_dir)
    t.record("run_start", text=f"fake run, data_path={data_path}")
    while True:
        a = await answer_source()
        if not a or not a.strip():
            break
        t.record("user", text=a)
    t.record("run_end", text="fake done")


def test_next_run_dir_is_sequential(tmp_path):
    a = webapp.next_run_dir(str(tmp_path))
    Path(a).mkdir(parents=True)
    b = webapp.next_run_dir(str(tmp_path))
    assert Path(a).name == "run_0"
    assert Path(b).name == "run_1"


def test_save_upload_writes_the_csv(tmp_path):
    rd = str(tmp_path / "run_0")
    p = webapp.save_upload(rd, "contigs.csv", b"barcode,cdr3\n")
    assert Path(p).read_bytes() == b"barcode,cdr3\n"
    assert p.endswith(".csv")


def test_list_runs_reports_event_counts(tmp_path):
    rd = str(tmp_path / "run_0")
    t = Transcript(rd)
    t.record("run_start", text="x")
    t.record("user", text="y")
    runs = webapp.list_runs(str(tmp_path))
    assert runs == [{"name": "run_0", "events": 2}]


def test_start_and_answer_bridge_drives_the_runner(tmp_path):
    state = webapp.AppState(str(tmp_path), runner=_fake_runner)
    run_dir = state.start(b"barcode,cdr3\n", "contigs.csv")
    # feed two answers then a blank to end the fake run
    state.answer("human")
    state.answer("A*02:01")
    state.answer("")
    # wait for the runner thread to drain the queue and finish
    deadline = time.time() + 5
    while state.is_running() and time.time() < deadline:
        time.sleep(0.02)
    assert _read_users(run_dir) == ["human", "A*02:01"]
    # data_path was threaded into the runner
    rows = Path(run_dir, "transcript.jsonl").read_text()
    assert "contigs.csv" in rows


def test_awaiting_flips_true_while_blocked_then_false(tmp_path):
    state = webapp.AppState(str(tmp_path), runner=_fake_runner)
    state.start(b"barcode,cdr3\n", "contigs.csv")
    # the fake runner reaches `await answer_source()` and blocks: awaiting is True
    deadline = time.time() + 5
    while not state.is_awaiting() and time.time() < deadline:
        time.sleep(0.02)
    assert state.is_awaiting() is True
    # ending the run clears awaiting
    state.answer("")
    deadline = time.time() + 5
    while state.is_running() and time.time() < deadline:
        time.sleep(0.02)
    assert state.is_awaiting() is False


def test_server_serves_landing_page(tmp_path):
    import threading, urllib.request
    httpd = webapp.build_app_server(str(tmp_path), port=0, runner=_fake_runner)
    threading.Thread(target=httpd.handle_request, daemon=True).start()
    port = httpd.server_address[1]
    body = urllib.request.urlopen(f"http://127.0.0.1:{port}/", timeout=5).read().decode()
    assert "<html" in body.lower()
    httpd.server_close()
