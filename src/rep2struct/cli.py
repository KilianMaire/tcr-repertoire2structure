from __future__ import annotations

import sys

from claude_agent_sdk import ClaudeSDKClient, AssistantMessage, TextBlock, ResultMessage

from . import intake
from .agents import build_options, intake_orchestrator_prompt, orchestrator_prompt


def plan_from_run_dir(run_dir: str) -> str:
    """intake when the interview has not run, else run. The seam the tests exercise."""
    return intake.next_phase(run_dir)


async def _drain(client):
    """Stream one agent turn to the terminal and return the concatenated assistant text."""
    text = []
    async for msg in client.receive_response():
        if isinstance(msg, AssistantMessage):
            for block in msg.content:
                if isinstance(block, TextBlock):
                    print(block.text, end="", flush=True)
                    text.append(block.text)
        if isinstance(msg, ResultMessage):
            break
    print(flush=True)
    return "".join(text)


async def run_session(run_dir: str, top_n: int = 8) -> None:
    """Phase 0/A/B driver. In the intake phase the intake-agent interviews the user turn by
    turn (stdin), then the orchestrator runs in handoff mode and stops at the artifacts. A
    rerun on the same run_dir enters the run phase and resumes through the checkpoint. top_n is
    the selection depth (how many top-ranked clonotypes to fold)."""
    if plan_from_run_dir(run_dir) == "intake":
        opts = build_options(run_dir, mode="handoff")
        async with ClaudeSDKClient(options=opts) as client:
            await client.query(
                f"Run the intake interview for run_dir {run_dir} using the intake-agent, "
                f"then call record_intake and proceed to the handoff orchestration.")
            # Interactive loop: the agent asks, the user answers on stdin, until the run ends.
            while True:
                await _drain(client)
                try:
                    answer = input("> ")
                except EOFError:
                    break
                if not answer.strip():
                    break
                await client.query(answer)
        return
    # run phase: resume from the persisted intake + checkpoint
    spec = intake.load_intake(run_dir)
    opts = build_options(run_dir, mode="handoff")
    prompt = (intake_orchestrator_prompt(run_dir, spec, top_n=top_n) if spec
              else orchestrator_prompt("", run_dir, top_n))
    async with ClaudeSDKClient(options=opts) as client:
        await client.query(
            prompt + "\nIf fold artifacts have already been run, call record_local_folds "
                     "first, then proceed to QC and the report.")
        await _drain(client)


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
