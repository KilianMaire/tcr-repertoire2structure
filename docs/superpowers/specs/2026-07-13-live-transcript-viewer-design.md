# Live transcript viewer design (2026-07-13)

## Purpose

Give the R2S agent pipeline a live, honest, watchable window. The multi agent
orchestration (orchestrator delegating to intake, structure strategist, fold
executors, QC, report) is R2S's differentiator, but today it only streams plain
assistant text to the terminal and hides every tool call and sub agent handoff.
This feature records the real message stream to a file and renders it as a live
timeline in the browser, so a run can be shown as it happens without staging
anything.

Primary near term use is the hackathon demo video, but the writer is a genuine
observability feature of the tool, shipped in the repo and tested.

## Honesty contract

The viewer renders only what the agents actually emitted. It tails a JSONL file
that the live run appends to. Nothing is hand fed. If the writer is disabled or
the file is absent, the viewer shows a waiting state rather than inventing events.

## Components

### 1. `src/rep2struct/transcript.py` (the writer)

`Transcript(run_dir)` with `record(kind, agent=None, tool=None, text="")`. Each
call appends one JSON line to `run_dir/transcript.jsonl` and flushes. Line shape:

```
{"seq": <int>, "kind": <str>, "agent": <str|null>, "tool": <str|null>, "text": <str>}
```

`kind` is one of: `user`, `agent_text`, `tool_use`, `tool_result`,
`subagent_start`, `subagent_stop`, `run_start`, `run_end`.

Rules:
- `seq` is a monotonically increasing counter derived from the current line
  count, so ordering is stable and the viewer can dedupe.
- `text` is truncated to a bounded length (tool args and results can be large).
- Best effort: every `record` is wrapped so a write failure never propagates.
  The pipeline is the product; losing a transcript line must not fail a run.
- No wall clock is required; `seq` carries order. A timestamp field is optional
  and omitted in v1 to keep the writer trivial and deterministic in tests.

### 2. Hook into `src/rep2struct/cli.py`

`_drain` already iterates the message stream. Extend it to also record, for each
block: `TextBlock` to `agent_text`, `ToolUseBlock` to `tool_use` (tool name plus
a one line arg summary), `ToolResultBlock` to `tool_result` (truncated). Record
each user query (the initial prompt and every stdin answer) as `user`. Register
`SubagentStart` and `SubagentStop` hooks in `build_options` to record
`subagent_start` and `subagent_stop` with the sub agent name when the SDK
surfaces them; the `tool_use` capture already shows delegations even if the hook
path is unavailable, so sub agent visibility does not depend on the hooks.

The writer is passed in (dependency injected) so `_drain` stays testable with a
fake stream and an in memory transcript.

### 3. `src/rep2struct/viewer/` (the page and server)

`python -m rep2struct.viewer <run_dir> [--port N]` starts a stdlib
`http.server` that serves `index.html` at `/` and the run's `transcript.jsonl`
at `/transcript.jsonl`. Using a local server avoids the `file://` fetch and CORS
restrictions. `index.html` is self contained (inline CSS and JS, no external
requests) and polls `/transcript.jsonl` every 400 ms, appending only lines whose
`seq` it has not yet rendered.

Rendering: a vertical live timeline. One card per event, colored by lane:
`you`, `orchestrator`, `intake`, `structure-strategist`, `executor`, `qc`,
`report`. Tool calls render as a compact pill (`build_fold_notebook -> path`),
sub agent handoffs as a divider (`orchestrator -> structure-strategist`). New
cards fade in, the page auto scrolls, a header shows the run dir and a live event
count. Dark theme, monospace for tool args, legible at 1080p screen capture.

## Failure handling

- Writer swallows all IO errors.
- Viewer server: if `transcript.jsonl` is missing, serve an empty array so the
  page shows the waiting state; return 200 with `[]` rather than 404.
- Malformed line in the JSONL: the page skips it and continues.

## Testing

- `transcript.py`: record writes well formed JSONL, seq increments, text is
  truncated at the bound, a forced IO error is swallowed (best effort), ordering
  is preserved across many records.
- `_drain` integration: feed a fake message stream containing a `TextBlock`, a
  `ToolUseBlock`, and a `ToolResultBlock`, assert the corresponding
  `agent_text`, `tool_use`, `tool_result` lines land in the transcript in order,
  and that user queries are recorded.
- Viewer server: smoke test that it serves `index.html` at `/` with 200 and
  returns the transcript (or `[]`) at `/transcript.jsonl`. The page itself is
  verified manually.

## Scope and non goals

- v1 does not persist timestamps or token usage; `seq` ordering only.
- v1 does not replay or seek; it appends live.
- No websocket; polling is sufficient for a single local viewer.
- The lane mapping is a small static dictionary from agent or tool name to lane;
  unknown agents fall back to a neutral lane.
