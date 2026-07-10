from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

from . import compute_routes


@dataclass
class IntakeSpec:
    data_type: str
    input_path: str
    question: str
    compute_route: str
    route_params: dict = field(default_factory=dict)


def _strip_secrets(route: str, params: dict) -> dict:
    """Drop any secret field the route declares (defence in depth: secrets should
    never have been put here, but never let one reach disk)."""
    secret = set(compute_routes.by_name(route).secret_fields)
    return {k: v for k, v in params.items() if k not in secret}


def save_intake(run_dir: str, spec: IntakeSpec) -> str:
    Path(run_dir).mkdir(parents=True, exist_ok=True)
    out = Path(run_dir) / "intake.json"
    payload = {
        "data_type": spec.data_type,
        "input_path": spec.input_path,
        "question": spec.question,
        "compute_route": spec.compute_route,
        "route_params": _strip_secrets(spec.compute_route, spec.route_params),
    }
    out.write_text(json.dumps(payload, indent=2))
    return str(out)


def load_intake(run_dir: str) -> IntakeSpec | None:
    p = Path(run_dir) / "intake.json"
    if not p.exists():
        return None
    d = json.loads(p.read_text())
    return IntakeSpec(d["data_type"], d["input_path"], d["question"],
                      d["compute_route"], d.get("route_params", {}))


def next_phase(run_dir: str) -> str:
    """intake when the interview has not run yet, else run (the orchestrator's own
    checkpoint decides fresh fold vs resume)."""
    return "run" if (Path(run_dir) / "intake.json").exists() else "intake"
