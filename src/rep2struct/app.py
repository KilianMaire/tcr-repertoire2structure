from __future__ import annotations
import asyncio
from claude_agent_sdk import query, ResultMessage
from .agents import build_options, orchestrator_prompt


async def run(csv_path, run_dir, top_n=8):
    opts = build_options(run_dir)
    result = None
    async for message in query(prompt=orchestrator_prompt(csv_path, run_dir, top_n), options=opts):
        if isinstance(message, ResultMessage) and getattr(message, "subtype", None) == "success":
            result = message.result
    return result


def main():
    import sys
    csv_path, run_dir = sys.argv[1], sys.argv[2]
    top_n = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    print(asyncio.run(run(csv_path, run_dir, top_n)))


if __name__ == "__main__":
    main()
