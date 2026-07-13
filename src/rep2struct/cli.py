from __future__ import annotations

import asyncio
import sys

from claude_agent_sdk import (
    ClaudeSDKClient, AssistantMessage, TextBlock, ToolUseBlock, ToolResultBlock,
    ResultMessage,
)

from . import intake
from .agents import build_options, intake_orchestrator_prompt, orchestrator_prompt
from .transcript import Transcript


def plan_from_run_dir(run_dir: str) -> str:
    """intake when the interview has not run, else run. The seam the tests exercise."""
    return intake.next_phase(run_dir)


def _tool_summary(block) -> str:
    """One line summary of a tool call's arguments for the timeline."""
    inp = getattr(block, "input", None) or {}
    return ", ".join(f"{k}={str(v)[:80]}" for k, v in inp.items())


def _result_summary(block) -> str:
    content = getattr(block, "content", "")
    if isinstance(content, list):
        content = " ".join(str(getattr(c, "text", c)) for c in content)
    return str(content)


def record_blocks(blocks, transcript, agent, is_assistant):
    """Dispatch one message's content blocks: print+return assistant text for the
    terminal, and record every block (text, tool use, tool result) to the
    transcript when one is provided. Returns the list of assistant text chunks."""
    text = []
    for block in blocks:
        if isinstance(block, TextBlock):
            if is_assistant:
                print(block.text, end="", flush=True)
                text.append(block.text)
            if transcript:
                transcript.record("agent_text", agent=agent, text=block.text)
        elif isinstance(block, ToolUseBlock):
            if transcript:
                transcript.record("tool_use", agent=agent, tool=block.name,
                                  text=_tool_summary(block))
        elif isinstance(block, ToolResultBlock):
            if transcript:
                transcript.record("tool_result", agent=agent,
                                  text=_result_summary(block))
    return text


async def _drain(client, transcript=None, agent="orchestrator"):
    """Stream one agent turn to the terminal and return the concatenated assistant
    text. When a transcript is given, every block is also recorded to it."""
    text = []
    async for msg in client.receive_response():
        blocks = getattr(msg, "content", None) or []
        text.extend(record_blocks(blocks, transcript, agent,
                                  isinstance(msg, AssistantMessage)))
        if isinstance(msg, ResultMessage):
            break
    print(flush=True)
    return "".join(text)


async def _stdin_answer():
    """Default answer source for the CLI: block on stdin off the event loop."""
    try:
        return await asyncio.to_thread(input, "> ")
    except EOFError:
        return None


async def _interactive_loop(client, transcript, answer_source):
    """Drive the intake conversation: drain one agent turn, ask the answer source
    for the next user answer, stop on a blank or absent answer. The answer source
    is injected (stdin for the CLI, a queue for the web app) so the loop itself is
    transport agnostic and unit testable."""
    while True:
        await _drain(client, transcript)
        answer = await answer_source()
        if not answer or not answer.strip():
            break
        transcript.record("user", text=answer)
        await client.query(answer)


async def run_session(run_dir: str, top_n: int = 8, answer_source=None,
                      data_path: str | None = None) -> None:
    """Phase 0/A/B driver. In the intake phase the intake-agent interviews the user turn by
    turn, then the orchestrator runs in handoff mode and stops at the artifacts. A rerun on the
    same run_dir enters the run phase and resumes through the checkpoint. top_n is the selection
    depth. answer_source injects where user answers come from (stdin by default, a queue for the
    web app); data_path, when set, tells the intake agent where the dropped file already is."""
    answer_source = answer_source or _stdin_answer
    transcript = Transcript(run_dir)
    if plan_from_run_dir(run_dir) == "intake":
        transcript.record("run_start", text=f"intake phase, run_dir {run_dir}")
        opts = build_options(run_dir, mode="handoff")
        async with ClaudeSDKClient(options=opts) as client:
            data_hint = (f" The data file has already been provided at {data_path}; do not ask "
                         f"for the path, confirm it and ask the remaining intake questions."
                         if data_path else "")
            kickoff = (
                f"Run the intake interview for run_dir {run_dir} using the intake-agent, then "
                f"call record_intake. After the brief is recorded, proceed to the handoff "
                f"orchestration: ingest_repertoire, annotate_specificity, prep_and_select with "
                f"top_n {top_n}, then the structure-strategist and each tool's executor to build "
                f"the fold artifacts, and stop." + data_hint)
            await client.query(kickoff)
            await _interactive_loop(client, transcript, answer_source)
        transcript.record("run_end", text="intake phase complete")
        return
    # run phase: resume from the persisted intake + checkpoint
    transcript.record("run_start", text=f"run phase, run_dir {run_dir}")
    spec = intake.load_intake(run_dir)
    opts = build_options(run_dir, mode="handoff")
    prompt = (intake_orchestrator_prompt(run_dir, spec, top_n=top_n) if spec
              else orchestrator_prompt("", run_dir, top_n))
    async with ClaudeSDKClient(options=opts) as client:
        await client.query(
            prompt + "\nIf fold artifacts have already been run, call record_local_folds "
                     "first, then proceed to QC and the report.")
        await _drain(client, transcript)
    transcript.record("run_end", text="run phase complete")


def main():
    import asyncio
    import os
    args = [a for a in sys.argv[1:]]
    if args and args[0] in ("-h", "--help"):
        print("usage: python -m rep2struct [RUN_DIR] [--top-n N]\n\n"
              "Runs the repertoire-to-structure session on RUN_DIR (default runs/session).\n"
              "First run: the intake agent interviews you (data, question, compute route),\n"
              "then builds fold artifacts and stops. Fold them, then rerun the SAME RUN_DIR\n"
              "to resume through QC and the report.\n\n"
              "  --top-n N   how many top-ranked clonotypes to fold (default 8; or env R2S_TOP_N)")
        return
    top_n = int(os.environ.get("R2S_TOP_N", "8"))
    if "--top-n" in args:
        i = args.index("--top-n")
        top_n = int(args[i + 1])
        del args[i:i + 2]
    run_dir = args[0] if args else "runs/session"
    asyncio.run(run_session(run_dir, top_n=top_n))


if __name__ == "__main__":
    main()
