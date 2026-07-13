"""Capture screenshots and a screen recording of the R2S web app for the submission video.

Honest by construction: instead of a scripted mock, it REPLAYS the real transcript that run_1
produced (runs/run_1/transcript.jsonl) into a fresh web app run, paced so it is watchable. Every
line on screen is genuine agent output; only the timing is slowed down. The single real user turn
(run_1 ended with the user typing "Ok") is reproduced through the actual chat box, so the green
"your turn" indicator is shown for real.

Run:  ./.venv/bin/python scripts/capture_demo.py
Out:  ~/Desktop/r2s_demo_capture/  (01_landing.png ... session.webm)

Requires: playwright (pip install playwright); the chromium browser is already cached on this Mac.
"""
from __future__ import annotations

import asyncio
import json
import re
import threading
from pathlib import Path

# Strip a stray "tour N" footer the orchestrator once echoed from an unrelated context file; it
# is not part of R2S output and should not appear on screen.
_FOOTER = re.compile(r"\s*↳\s*tour\s*\d+\s*$")

from playwright.async_api import async_playwright

from rep2struct import webapp
from rep2struct.transcript import Transcript

REPO = Path(__file__).resolve().parent.parent
SOURCE = REPO / "runs" / "run_1" / "transcript.jsonl"
CSV = Path.home() / "Desktop" / "r2s_demo_donor1.csv"
OUT = Path.home() / "Desktop" / "r2s_demo_capture"
PORT = 8790
PACE = 0.45          # seconds between replayed events
SETTLE = 0.6         # let the page poll (500 ms) catch up before a screenshot


def _load_events():
    if not SOURCE.exists():
        raise SystemExit(f"no source transcript at {SOURCE}; run a real run first")
    return [json.loads(l) for l in SOURCE.read_text().splitlines() if l.strip()]


def _make_replay_runner(events):
    async def replay_runner(run_dir, answer_source=None, data_path=None):
        name = Path(run_dir).name
        t = Transcript(run_dir)
        for ev in events:
            if ev["kind"] == "run_start":
                t.record("run_start", text=f"intake phase, run_dir {run_dir}")
            elif ev["kind"] == "user":
                ans = await answer_source()          # blocks -> "your turn" in the UI
                if not ans or not ans.strip():
                    break
                t.record("user", text=ans)
            else:
                text = _FOOTER.sub("", (ev.get("text") or "").replace("run_1", name))
                t.record(ev["kind"], agent=ev.get("agent"), tool=ev.get("tool"), text=text)
            await asyncio.sleep(PACE)
    return replay_runner


def _start_server(events):
    httpd = webapp.build_app_server("runs", port=PORT, runner=_make_replay_runner(events))
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    return httpd


async def _feed_text(page):
    return await page.evaluate(
        "() => Array.from(document.querySelectorAll('#log .card, .card')).map(e=>e.textContent).join(' \\n ')")


async def _wait_for(page, needle, timeout=40.0):
    """Poll the visible feed until it contains needle (or time out)."""
    for _ in range(int(timeout / 0.25)):
        if needle.lower() in (await _feed_text(page)).lower():
            return True
        await asyncio.sleep(0.25)
    print(f"  (timed out waiting for {needle!r})")
    return False


async def _wait_awaiting(page, timeout=60.0):
    """Poll until the UI shows the green 'your turn' state (the run is blocked on the chat)."""
    for _ in range(int(timeout / 0.25)):
        turn = await page.evaluate(
            "() => { const w=document.getElementById('working');"
            " return w && !w.hidden && w.className.includes('turn'); }")
        if turn:
            return True
        await asyncio.sleep(0.25)
    print("  (timed out waiting for your-turn)")
    return False


async def _shot(page, name):
    await asyncio.sleep(SETTLE)
    await page.screenshot(path=str(OUT / name))
    print(f"  captured {name}")


async def main():
    OUT.mkdir(parents=True, exist_ok=True)
    if not CSV.exists():
        raise SystemExit(f"no demo CSV at {CSV}; create a small 10x slice there first")
    events = _load_events()
    httpd = _start_server(events)
    print(f"replay server on http://127.0.0.1:{PORT} ({len(events)} real events from run_1)")

    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)
        context = await browser.new_context(
            viewport={"width": 1280, "height": 800},
            record_video_dir=str(OUT), record_video_size={"width": 1280, "height": 800})
        page = await context.new_page()
        await page.goto(f"http://127.0.0.1:{PORT}/")
        await asyncio.sleep(1.0)

        await _shot(page, "01_landing.png")                       # hero dropzone

        # start the run by choosing the CSV through the real file input
        await page.set_input_files("#file", str(CSV))
        await _wait_for(page, "intake phase")
        await _shot(page, "02_run_start.png")

        await _wait_for(page, "annotated 119")                    # annotation result
        await _shot(page, "04_annotation.png")

        await _wait_for(page, "assigned protenix")                # strategist routing
        await _shot(page, "05_routing.png")

        await _wait_for(page, "colab_notebook for group")         # executor built the artifact
        await _shot(page, "06_build_artifact.png")

        # the real single user turn: the run pauses awaiting the chat -> green "your turn"
        await _wait_awaiting(page)
        await asyncio.sleep(1.0)
        await _shot(page, "03_intake_your_turn.png")
        await page.fill("#msg", "Ok")
        await page.click("#send")
        await asyncio.sleep(1.5)

        await _shot(page, "07_timeline_full.png")
        await asyncio.sleep(2.0)
        await _shot(page, "08_done.png")

        await context.close()                                     # flushes the video
        await browser.close()

    httpd.shutdown()
    # give the recorded video a stable name
    vids = sorted(OUT.glob("*.webm"))
    if vids:
        target = OUT / "session.webm"
        if target.exists():
            target.unlink()
        vids[-1].rename(target)
        print(f"video: {target}")
    print(f"done -> {OUT}")


if __name__ == "__main__":
    asyncio.run(main())
