# Intake agent and plastic compute routes Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Put a plastic conversational intake agent in front of the R2S orchestrator that asks data, question, and compute environment, then builds the right fold artifact for the chosen route and stops for the user to fold.

**Architecture:** A new open `compute_routes` registry (mirroring `structure_tools`) and an `IntakeSpec` carry the interview result. A route aware `build_fold_artifact` tool emits a Colab notebook or a local bash script. A handoff mode on the orchestrator skips Playwright: executors build the artifact and stop. A `ClaudeSDKClient` CLI runs the multi turn interview, and a rerun on the same run dir resumes through checkpoint.

**Tech Stack:** Python 3.11, `claude-agent-sdk` 0.2.113 (`tool`, `create_sdk_mcp_server`, `ClaudeSDKClient`, `ClaudeAgentOptions`, `AgentDefinition`), pytest. Reuses `structure_tools`, `tools/protenix_inputs`, `tools/notebook`, `runstate`, `agent_tools`.

## Global Constraints

- Docs and any generated prose obey the no dash punctuation rule (use periods, commas, parentheses).
- This repo KEEPS the `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>` trailer on commits (hackathon repo).
- Tests run with `./.venv/bin/python -m pytest -q` from the repo root.
- Branch: `feat/intake-agent-compute-routes` (already created off `main`, spec already committed).
- Secrets rule: an SSH password collected in the interview lives only in session memory, never written to `run_dir/intake.json` or any log or artifact.
- Never Boltz. Protenix stays the default workhorse. The existing autonomous Playwright path stays intact (additive only).
- Registry honesty: an unwired compute route must yield the fallback bash script plus a plain not wired statement, never a fabricated remote run.

---

### Task 1: `compute_routes` open registry

**Files:**
- Create: `src/rep2struct/compute_routes.py`
- Test: `tests/test_compute_routes.py`

**Interfaces:**
- Produces: `ComputeRoute` dataclass; `REGISTRY: list[ComputeRoute]`; `get_default() -> ComputeRoute`; `by_name(name: str) -> ComputeRoute`; `recommend(context: str = "") -> ComputeRoute`; `artifact_kind_for(name: str) -> str`; `is_wired(name: str) -> bool`; `required_fields_for(name: str) -> tuple[str, ...]`; `as_dicts() -> list[dict]`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_compute_routes.py
from rep2struct import compute_routes as cr


def test_default_is_colab_the_simplest():
    d = cr.get_default()
    assert d.name == "colab" and d.is_default


def test_registry_has_the_four_v1_routes():
    assert {r.name for r in cr.REGISTRY} == {"colab", "local_gpu", "ssh", "server"}


def test_colab_and_local_gpu_are_wired_ssh_and_server_are_not():
    assert cr.is_wired("colab") and cr.is_wired("local_gpu")
    assert not cr.is_wired("ssh") and not cr.is_wired("server")


def test_artifact_kind_per_route():
    assert cr.artifact_kind_for("colab") == "colab_notebook"
    assert cr.artifact_kind_for("local_gpu") == "bash_script"
    assert cr.artifact_kind_for("ssh") == "bash_script"


def test_ssh_requires_connection_fields_and_marks_password_secret():
    ssh = cr.by_name("ssh")
    assert ssh.required_fields == ("host", "user", "remote_path")
    assert ssh.secret_fields == ("password",)
    # colab needs nothing extra
    assert cr.by_name("colab").required_fields == ()


def test_recommend_falls_back_to_default_when_user_does_not_know():
    assert cr.recommend("I don't know").name == "colab"
    assert cr.recommend("").name == "colab"


def test_as_dicts_exposes_fields_but_never_a_secret_value():
    d = {x["name"]: x for x in cr.as_dicts()}
    assert d["ssh"]["required_fields"] == ["host", "user", "remote_path"]
    assert d["ssh"]["secret_fields"] == ["password"]
    assert d["ssh"]["wired"] is False and d["colab"]["is_default"] is True
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./.venv/bin/python -m pytest tests/test_compute_routes.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'rep2struct.compute_routes'`

- [ ] **Step 3: Write minimal implementation**

```python
# src/rep2struct/compute_routes.py
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ComputeRoute:
    name: str
    description: str
    required_fields: tuple[str, ...]   # non secret connection fields the agent must collect
    secret_fields: tuple[str, ...]     # collected in session memory only, never persisted
    artifact_kind: str                 # "colab_notebook" | "bash_script"
    wired: bool                        # True: R2S produces a runnable artifact the user runs
    is_default: bool = False


REGISTRY: list[ComputeRoute] = [
    ComputeRoute(
        name="colab",
        description="Google Colab. The simplest option: R2S writes a self contained "
                    "notebook the user uploads and runs on a Colab GPU.",
        required_fields=(),
        secret_fields=(),
        artifact_kind="colab_notebook",
        wired=True,
        is_default=True,
    ),
    ComputeRoute(
        name="local_gpu",
        description="A GPU machine the user has a shell on (for example a local H100). "
                    "R2S writes a bash script that folds with Protenix and leaves the "
                    "CIFs on that machine.",
        required_fields=("working_path",),
        secret_fields=(),
        artifact_kind="bash_script",
        wired=True,
    ),
    ComputeRoute(
        name="ssh",
        description="A remote host reached over SSH. R2S collects the connection details "
                    "and hands the user a job script plus scp/sbatch instructions. The SSH "
                    "runner itself is not wired yet.",
        required_fields=("host", "user", "remote_path"),
        secret_fields=("password",),
        artifact_kind="bash_script",
        wired=False,
    ),
    ComputeRoute(
        name="server",
        description="A shared server the user names by address. Same handoff as ssh; the "
                    "server runner is not wired yet.",
        required_fields=("address", "path"),
        secret_fields=(),
        artifact_kind="bash_script",
        wired=False,
    ),
]


def get_default() -> ComputeRoute:
    return next(r for r in REGISTRY if r.is_default)


def by_name(name: str) -> ComputeRoute:
    return next(r for r in REGISTRY if r.name == name)


def recommend(context: str = "") -> ComputeRoute:
    """The simplest route for a user who does not know: the default (Colab)."""
    return get_default()


def artifact_kind_for(name: str) -> str:
    return by_name(name).artifact_kind


def is_wired(name: str) -> bool:
    return by_name(name).wired


def required_fields_for(name: str) -> tuple[str, ...]:
    return by_name(name).required_fields


def as_dicts() -> list[dict]:
    return [
        {
            "name": r.name,
            "description": r.description,
            "required_fields": list(r.required_fields),
            "secret_fields": list(r.secret_fields),
            "artifact_kind": r.artifact_kind,
            "wired": r.wired,
            "is_default": r.is_default,
        }
        for r in REGISTRY
    ]
```

- [ ] **Step 4: Run test to verify it passes**

Run: `./.venv/bin/python -m pytest tests/test_compute_routes.py -q`
Expected: PASS (7 passed)

- [ ] **Step 5: Commit**

```bash
git add src/rep2struct/compute_routes.py tests/test_compute_routes.py
git commit -m "feat(intake): open compute-route registry (colab/local_gpu/ssh/server)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: `IntakeSpec` persistence (no secrets) and phase detection

**Files:**
- Create: `src/rep2struct/intake.py`
- Test: `tests/test_intake.py`

**Interfaces:**
- Consumes: `compute_routes.by_name` (to know a route's `secret_fields`).
- Produces: `IntakeSpec` dataclass with fields `data_type: str, input_path: str, question: str, compute_route: str, route_params: dict`; `save_intake(run_dir: str, spec: IntakeSpec) -> str` (returns the json path, strips any secret field from `route_params` before writing); `load_intake(run_dir: str) -> IntakeSpec | None`; `next_phase(run_dir: str) -> str` returning `"intake"` when no `intake.json` exists, else `"run"`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_intake.py
from rep2struct import intake
from rep2struct.intake import IntakeSpec


def test_round_trip(tmp_path):
    spec = IntakeSpec("10x_vdj", "/data/contigs.csv", "which epitope?",
                      "local_gpu", {"working_path": "/scratch/run"})
    p = intake.save_intake(str(tmp_path), spec)
    assert p.endswith("intake.json")
    got = intake.load_intake(str(tmp_path))
    assert got == spec


def test_secret_never_written_to_disk(tmp_path):
    # A password slipped into route_params must NOT reach intake.json.
    spec = IntakeSpec("10x_vdj", "/data/c.csv", "q", "ssh",
                      {"host": "hpc.uni.dk", "user": "kilian",
                       "remote_path": "/work", "password": "hunter2"})
    p = intake.save_intake(str(tmp_path), spec)
    with open(p) as fh:
        raw = fh.read()
    assert "hunter2" not in raw
    assert "password" not in raw
    got = intake.load_intake(str(tmp_path))
    assert "password" not in got.route_params
    assert got.route_params["host"] == "hpc.uni.dk"


def test_next_phase_is_intake_when_no_file_then_run(tmp_path):
    assert intake.next_phase(str(tmp_path)) == "intake"
    intake.save_intake(str(tmp_path),
                       IntakeSpec("d", "/i.csv", "q", "colab", {}))
    assert intake.next_phase(str(tmp_path)) == "run"


def test_load_returns_none_when_absent(tmp_path):
    assert intake.load_intake(str(tmp_path)) is None
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./.venv/bin/python -m pytest tests/test_intake.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'rep2struct.intake'`

- [ ] **Step 3: Write minimal implementation**

```python
# src/rep2struct/intake.py
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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `./.venv/bin/python -m pytest tests/test_intake.py -q`
Expected: PASS (4 passed)

- [ ] **Step 5: Commit**

```bash
git add src/rep2struct/intake.py tests/test_intake.py
git commit -m "feat(intake): IntakeSpec persistence (secret-stripping) + phase detection

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: Protenix bash script builder (local_gpu / handoff artifact)

**Files:**
- Create: `src/rep2struct/tools/protenix_script.py`
- Test: `tests/test_protenix_script.py`

**Interfaces:**
- Consumes: nothing new (it is handed the already shaped `{key: <protenix JSON>}` inputs by `_fold_inputs`, exactly like `_protenix_notebook`).
- Produces: `build(inputs: dict, working_path: str = ".") -> str` returning a self contained bash script. The script writes each embedded record to `inputs/{key}.json` and folds it with the SAME proven command as `_protenix_notebook` (`protenix pred -i inputs/{key}.json -o out/{key} -s 101 -n protenix_base_default_v1.0.0 --use_default_params true`), leaving CIFs under `out/{key}` (the `{cid}_cognate` / `{cid}_scramble` layout QC needs). No repatriation (CIFs are already local).

- [ ] **Step 1: Write the failing test**

```python
# tests/test_protenix_script.py
from rep2struct.tools import protenix_inputs, protenix_script

FASTA = (">A\nAAAA\n>B\nBBBB\n>C\nCCCC\n>D\nDDDD\n>E\nGILGFVFTL\n")


def _inputs():
    built = protenix_inputs.build(FASTA)
    return {f"c0_{k}": v for k, v in built.items()}


def test_script_is_bash_and_writes_inputs_and_folds():
    s = protenix_script.build(_inputs(), working_path="/scratch/run")
    assert s.startswith("#!/usr/bin/env bash")
    assert "set -euo pipefail" in s
    # both constructs are referenced
    assert "c0_cognate" in s and "c0_scramble" in s
    # the proven fold command and the local out/ layout QC reads
    assert "protenix pred -i inputs/" in s
    assert "-o out/" in s
    assert "--use_default_params true" in s
    # runs in the working path the user gave, no browser repatriation
    assert "/scratch/run" in s
    assert "files.download" not in s


def test_embeds_input_json_so_no_external_files_needed():
    s = protenix_script.build(_inputs())
    assert "GILGFVFTL" in s  # the peptide from the embedded record
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./.venv/bin/python -m pytest tests/test_protenix_script.py -q`
Expected: FAIL with `ImportError: cannot import name 'protenix_script'`

- [ ] **Step 3: Write minimal implementation**

```python
# src/rep2struct/tools/protenix_script.py
from __future__ import annotations

import json

TOOL = "protenix"

# The proven Protenix invocation, identical to tools/notebook._protenix_notebook cell 3.
_FOLD_CMD = ("protenix pred -i inputs/{key}.json -o out/{key} -s 101 "
             "-n protenix_base_default_v1.0.0 --use_default_params true")


def build(inputs: dict, working_path: str = ".") -> str:
    """A self contained bash script that folds each embedded Protenix record on a machine
    the user has a shell on. INPUTS is {key: <protenix prediction JSON>} (cognate + scramble,
    keys prefixed by clonotype id). CIFs land under out/{key} in working_path, the same
    {cid}_cognate / {cid}_scramble layout the CDR3-peptide QC calibrates on. MSA free by
    design (mirrors the documented reliable Protenix path); no browser repatriation because
    the outputs are already local."""
    embedded = json.dumps(inputs)
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        f"cd {working_path}",
        "pip install -q protenix",
        "mkdir -p inputs out",
        "python - <<'PY'",
        "import json, os",
        f"INPUTS = json.loads({embedded!r})",
        "os.makedirs('inputs', exist_ok=True)",
        "for key, obj in INPUTS.items():",
        "    json.dump(obj, open(f'inputs/{key}.json', 'w'))",
        "print('wrote', len(INPUTS), 'inputs:', sorted(INPUTS))",
        "PY",
        "for f in inputs/*.json; do",
        '  key=$(basename "$f" .json)',
        f"  {_FOLD_CMD.format(key='$key')}",
        "done",
        "echo DONE",
    ]
    return "\n".join(lines) + "\n"
```

- [ ] **Step 4: Run test to verify it passes**

Run: `./.venv/bin/python -m pytest tests/test_protenix_script.py -q`
Expected: PASS (2 passed)

- [ ] **Step 5: Commit**

```bash
git add src/rep2struct/tools/protenix_script.py tests/test_protenix_script.py
git commit -m "feat(intake): Protenix bash-script builder for the local_gpu route

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: route aware `build_fold_artifact` MCP tool

**Files:**
- Modify: `src/rep2struct/agent_tools.py` (add the tool near `build_fold_notebook:161-189`, register it in `build_server:311-316`)
- Test: `tests/test_build_fold_artifact.py`

**Interfaces:**
- Consumes: `compute_routes.artifact_kind_for`, `compute_routes.is_wired`; the existing `_fold_inputs`, `tools.notebook.build_notebook`, `tools.protenix_script.build`; `RunState`; the persisted `foldjobs` stage.
- Produces: MCP tool `build_fold_artifact` with schema `{run_dir: str, clonotype_id: str, tool: str, compute_route: str}`. It dispatches on the route's `artifact_kind`: `colab_notebook` writes `<run_dir>/notebooks/{cid}_{tool}.ipynb` (JSON notebook) and `bash_script` writes `<run_dir>/scripts/{cid}_{tool}.sh`. For an unwired route it still writes the bash script and sets `structuredContent.route_wired=False` with a plain not wired note. Returns `structuredContent` with `artifact_path`, `artifact_kind`, `route_wired`, `clonotype_id`, `tool`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_build_fold_artifact.py
import asyncio
from pathlib import Path

from rep2struct import agent_tools
from rep2struct.runstate import RunState


def _seed_job(rd):
    fasta = ">A\nAAAA\n>B\nBBBB\n>C\nCCCC\n>D\nDDDD\n>E\nGILGFVFTL\n"
    RunState(rd).write_stage("foldjobs", [{"clonotype_id": "c0",
                                           "construct_fasta": fasta,
                                           "group_id": "g0"}])


def _call(rd, route):
    return asyncio.run(agent_tools.build_fold_artifact.handler(
        {"run_dir": rd, "clonotype_id": "c0", "tool": "protenix",
         "compute_route": route}))


def test_colab_route_writes_a_notebook(tmp_path):
    rd = str(tmp_path); _seed_job(rd)
    r = _call(rd, "colab")
    sc = r["structuredContent"]
    assert sc["artifact_kind"] == "colab_notebook"
    assert sc["route_wired"] is True
    assert sc["artifact_path"].endswith("c0_protenix.ipynb")
    assert Path(sc["artifact_path"]).exists()


def test_local_gpu_route_writes_a_bash_script(tmp_path):
    rd = str(tmp_path); _seed_job(rd)
    r = _call(rd, "local_gpu")
    sc = r["structuredContent"]
    assert sc["artifact_kind"] == "bash_script"
    assert sc["route_wired"] is True
    assert sc["artifact_path"].endswith("c0_protenix.sh")
    body = Path(sc["artifact_path"]).read_text()
    assert body.startswith("#!/usr/bin/env bash")
    assert "protenix pred" in body


def test_unwired_ssh_route_still_scripts_but_flags_not_wired(tmp_path):
    rd = str(tmp_path); _seed_job(rd)
    r = _call(rd, "ssh")
    sc = r["structuredContent"]
    assert sc["artifact_kind"] == "bash_script"
    assert sc["route_wired"] is False
    assert Path(sc["artifact_path"]).exists()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./.venv/bin/python -m pytest tests/test_build_fold_artifact.py -q`
Expected: FAIL with `AttributeError: module 'rep2struct.agent_tools' has no attribute 'build_fold_artifact'`

- [ ] **Step 3: Write minimal implementation**

Add this import at the top of `agent_tools.py` with the other imports:

```python
from . import compute_routes
```

Add the tool after `build_fold_notebook` (after line 189):

```python
@tool("build_fold_artifact",
      "Build the fold artifact (Colab notebook or local bash script) for one clonotype, "
      "chosen by the compute route, write it under the run dir, and return its path.",
      {"run_dir": str, "clonotype_id": str, "tool": str, "compute_route": str})
async def build_fold_artifact(args):
    import json as _json
    from .tools.notebook import build_notebook
    from .tools.protenix_script import build as build_script
    rs = RunState(args["run_dir"])
    jobs = rs.read_stage("foldjobs") if rs.stage_done("foldjobs") else []
    cid = args["clonotype_id"]
    job = next((j for j in jobs if j["clonotype_id"] == cid), None)
    if job is None:
        return _txt(f"no fold job for {cid}")
    tool = args["tool"]
    route = args["compute_route"]
    kind = compute_routes.artifact_kind_for(route)
    wired = compute_routes.is_wired(route)
    # tcrdock needs the gene-level Clonotype+Annotation (see build_fold_notebook); read them
    # back only when needed so the artifact inputs match the notebook path exactly.
    clon = ann = None
    if tool == "tcrdock":
        clon = next((c for c in _load(args["run_dir"], "ingest", Clonotype) if c.id == cid), None)
        ann = next((a for a in _load(args["run_dir"], "annotate", Annotation)
                    if a.clonotype_id == cid), None)
    inputs = _fold_inputs(tool, job, cid, clon, ann)
    if kind == "colab_notebook":
        nb = build_notebook(tool, inputs)
        nb_dir = Path(args["run_dir"]) / "notebooks"
        nb_dir.mkdir(parents=True, exist_ok=True)
        out = nb_dir / f"{cid}_{tool}.ipynb"
        out.write_text(_json.dumps(nb, indent=1))
    else:  # bash_script (local_gpu, and the honest ssh/server handoff)
        working = job.get("working_path") or "."
        script = build_script(inputs, working_path=working)
        sc_dir = Path(args["run_dir"]) / "scripts"
        sc_dir.mkdir(parents=True, exist_ok=True)
        out = sc_dir / f"{cid}_{tool}.sh"
        out.write_text(script)
    note = "" if wired else f" (route '{route}' runner not wired; run the script yourself)"
    r = _txt(f"{kind} for {cid} ({tool}) via {route} written to {out}{note}")
    r["structuredContent"] = {"artifact_path": str(out), "artifact_kind": kind,
                              "route_wired": wired, "clonotype_id": cid, "tool": tool}
    return r
```

Register it in `build_server` (add to the tools list at line 312):

```python
    return create_sdk_mcp_server(name="rep2struct", version="0.1.0", tools=[
        ingest_repertoire, annotate_specificity, prep_and_select, list_fold_jobs,
        list_structure_tools, build_fold_notebook, build_fold_artifact, record_fold_result,
        qc_structure, render_final_report,
    ])
```

- [ ] **Step 4: Run test to verify it passes**

Run: `./.venv/bin/python -m pytest tests/test_build_fold_artifact.py -q`
Expected: PASS (3 passed)

- [ ] **Step 5: Commit**

```bash
git add src/rep2struct/agent_tools.py tests/test_build_fold_artifact.py
git commit -m "feat(intake): route-aware build_fold_artifact (notebook vs bash script)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 5: `record_local_folds` MCP tool for checkpoint resume

**Files:**
- Modify: `src/rep2struct/agent_tools.py` (add helper `scan_recorded_folds` and the tool; register in `build_server`)
- Test: `tests/test_record_local_folds.py`

**Interfaces:**
- Consumes: `RunState`; the `folds` stage shape written by `record_fold_result` (`{cid: {"paths": [...], "tool": ...}}`).
- Produces: module function `scan_recorded_folds(run_dir: str, tool: str = "protenix") -> dict` that globs `<run_dir>/out/*_cognate/**/*.cif` and `*_scramble` siblings and returns `{cid: {"paths": [...], "tool": tool}}` keyed by the clonotype id parsed from the `{cid}_cognate` directory name; and MCP tool `record_local_folds` with schema `{run_dir: str, tool: str}` that merges the scan into the `folds` stage and returns how many clonotypes it recorded. This is the resume path for the local_gpu route (CIFs already on disk) and for a Colab download the user unzipped into `<run_dir>/out`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_record_local_folds.py
import asyncio
from pathlib import Path

from rep2struct import agent_tools
from rep2struct.runstate import RunState


def _write_cif(p):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("data_x\n_atom_site.group_PDB\n")


def test_scan_finds_cognate_and_scramble_by_directory(tmp_path):
    out = tmp_path / "out"
    _write_cif(out / "c0_cognate" / "preds" / "cognate_sample_0.cif")
    _write_cif(out / "c0_scramble" / "preds" / "scramble_sample_0.cif")
    _write_cif(out / "c1_cognate" / "preds" / "cognate_sample_0.cif")
    found = agent_tools.scan_recorded_folds(str(tmp_path))
    assert set(found) == {"c0", "c1"}
    c0 = [p for p in found["c0"]["paths"]]
    assert any("c0_cognate" in p for p in c0)
    assert any("c0_scramble" in p for p in c0)
    assert found["c0"]["tool"] == "protenix"


def test_record_local_folds_writes_the_folds_stage(tmp_path):
    out = tmp_path / "out"
    _write_cif(out / "c0_cognate" / "preds" / "s0.cif")
    r = asyncio.run(agent_tools.record_local_folds.handler(
        {"run_dir": str(tmp_path), "tool": "protenix"}))
    assert r["structuredContent"]["recorded"] == 1
    done = RunState(str(tmp_path)).read_stage("folds")
    assert "c0" in done and done["c0"]["tool"] == "protenix"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./.venv/bin/python -m pytest tests/test_record_local_folds.py -q`
Expected: FAIL with `AttributeError: module 'rep2struct.agent_tools' has no attribute 'scan_recorded_folds'`

- [ ] **Step 3: Write minimal implementation**

Add near the other module helpers in `agent_tools.py` (after `_scramble_null`):

```python
def scan_recorded_folds(run_dir, tool="protenix"):
    """Find fold outputs already on disk under <run_dir>/out and group them per clonotype.
    Protenix marks the construct in the DIRECTORY ({cid}_cognate / {cid}_scramble), not the
    filename, so parse the cid from the top-level out/ dir. Resume path: the local_gpu bash
    route wrote CIFs here directly, and a Colab download unzipped here has the same layout."""
    out = Path(run_dir) / "out"
    found = {}
    if not out.exists():
        return found
    for d in sorted(out.iterdir()):
        if not d.is_dir():
            continue
        for suffix in ("_cognate", "_scramble"):
            if d.name.endswith(suffix):
                cid = d.name[: -len(suffix)]
                cifs = sorted(str(p) for p in d.rglob("*.cif"))
                if cifs:
                    found.setdefault(cid, {"paths": [], "tool": tool})
                    found[cid]["paths"].extend(cifs)
    return found
```

Add the tool after `record_fold_result` (after line 200):

```python
@tool("record_local_folds",
      "Scan <run_dir>/out for fold CIFs already on disk (local_gpu run, or a Colab download "
      "unzipped there) and record them per clonotype so QC can proceed.",
      {"run_dir": str, "tool": str})
async def record_local_folds(args):
    rs = RunState(args["run_dir"])
    done = rs.read_stage("folds") if rs.stage_done("folds") else {}
    found = scan_recorded_folds(args["run_dir"], args.get("tool", "protenix"))
    done.update(found)
    rs.write_stage("folds", done)
    r = _txt(f"recorded {len(found)} clonotypes from disk: {sorted(found)}")
    r["structuredContent"] = {"recorded": len(found), "clonotypes": sorted(found)}
    return r
```

Register it in `build_server` (append `record_local_folds` to the tools list).

- [ ] **Step 4: Run test to verify it passes**

Run: `./.venv/bin/python -m pytest tests/test_record_local_folds.py -q`
Expected: PASS (2 passed)

- [ ] **Step 5: Commit**

```bash
git add src/rep2struct/agent_tools.py tests/test_record_local_folds.py
git commit -m "feat(intake): record_local_folds resume path (scan out/ for CIFs)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 6: `list_compute_routes` intake tool

**Files:**
- Modify: `src/rep2struct/agent_tools.py` (add the tool; register in `build_server`)
- Test: `tests/test_list_compute_routes.py`

**Interfaces:**
- Consumes: `compute_routes.as_dicts`.
- Produces: MCP tool `list_compute_routes` with schema `{run_dir: str}` (run_dir kept for tool call uniformity, unused) returning `structuredContent={"routes": compute_routes.as_dicts()}`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_list_compute_routes.py
import asyncio
from rep2struct import agent_tools


def test_lists_the_routes_with_required_and_wired_flags(tmp_path):
    r = asyncio.run(agent_tools.list_compute_routes.handler({"run_dir": str(tmp_path)}))
    routes = {x["name"]: x for x in r["structuredContent"]["routes"]}
    assert set(routes) == {"colab", "local_gpu", "ssh", "server"}
    assert routes["colab"]["is_default"] is True
    assert routes["ssh"]["required_fields"] == ["host", "user", "remote_path"]
    assert routes["ssh"]["wired"] is False
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./.venv/bin/python -m pytest tests/test_list_compute_routes.py -q`
Expected: FAIL with `AttributeError: module 'rep2struct.agent_tools' has no attribute 'list_compute_routes'`

- [ ] **Step 3: Write minimal implementation**

Add after `list_structure_tools` (after line 128) in `agent_tools.py`:

```python
@tool("list_compute_routes",
      "List the compute routes (Colab, local GPU, SSH, server), the fields each needs, and "
      "whether its runner is wired, so the intake agent asks the right questions.",
      {"run_dir": str})
async def list_compute_routes(args):
    r = _txt("compute route registry")
    r["structuredContent"] = {"routes": compute_routes.as_dicts()}
    return r
```

Register it in `build_server` (append `list_compute_routes`).

- [ ] **Step 4: Run test to verify it passes**

Run: `./.venv/bin/python -m pytest tests/test_list_compute_routes.py -q`
Expected: PASS (1 passed)

- [ ] **Step 5: Commit**

```bash
git add src/rep2struct/agent_tools.py tests/test_list_compute_routes.py
git commit -m "feat(intake): list_compute_routes tool for the intake agent

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 7: intake agent, handoff executors, and `build_options(mode)`

**Files:**
- Modify: `src/rep2struct/agents.py` (add `intake-agent`, a handoff executor prompt, a `mode` param on `build_options`, and a handoff orchestrator prompt)
- Test: `tests/test_agents_intake.py`

**Interfaces:**
- Consumes: `build_agents`, `_executor`, `build_options`, `orchestrator_prompt` (existing).
- Produces: `build_agents()` gains an `"intake-agent"` (opus) whose tools are `["mcp__rep2struct__list_compute_routes", "mcp__rep2struct__ingest_repertoire", "Agent"]` (no Playwright); `_executor(name, tool, mode="auto")` where `mode="handoff"` swaps the Playwright drive for a `build_fold_artifact` then stop prompt; `build_options(run_dir, mode="auto")` where `mode="handoff"` builds handoff executors and DROPS every `mcp__playwright__*` entry from `allowed_tools` while adding `mcp__rep2struct__build_fold_artifact` and `mcp__rep2struct__list_compute_routes` and `mcp__rep2struct__record_local_folds`; `intake_orchestrator_prompt(run_dir, spec)` that threads the IntakeSpec into the existing `orchestrator_prompt`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_agents_intake.py
from rep2struct import agents


def test_intake_agent_exists_without_playwright():
    a = agents.build_agents()
    assert "intake-agent" in a
    assert not any("playwright" in t for t in a["intake-agent"].tools)


def test_handoff_options_drop_playwright_and_add_artifact_tool():
    opts = agents.build_options("/tmp/run", mode="handoff")
    assert not any("playwright" in t for t in opts.allowed_tools)
    assert "mcp__rep2struct__build_fold_artifact" in opts.allowed_tools
    assert "mcp__rep2struct__record_local_folds" in opts.allowed_tools


def test_auto_options_keep_playwright_unchanged():
    opts = agents.build_options("/tmp/run")  # default mode="auto"
    assert any("playwright" in t for t in opts.allowed_tools)


def test_handoff_executor_prompt_builds_artifact_and_does_not_drive_colab():
    ex = agents._executor("protenix-agent", "protenix", mode="handoff")
    assert "build_fold_artifact" in ex.prompt
    assert "playwright" not in ex.prompt.lower()
    assert "Ctrl+Enter" not in ex.prompt


def test_intake_orchestrator_prompt_threads_the_spec():
    from rep2struct.intake import IntakeSpec
    spec = IntakeSpec("10x", "/data/c.csv", "which epitope?", "local_gpu",
                      {"working_path": "/scratch"})
    p = agents.intake_orchestrator_prompt("/tmp/run", spec)
    assert "/data/c.csv" in p and "which epitope?" in p and "local_gpu" in p
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./.venv/bin/python -m pytest tests/test_agents_intake.py -q`
Expected: FAIL with `AttributeError` on `intake-agent` / `mode` / `intake_orchestrator_prompt`

- [ ] **Step 3: Write minimal implementation**

In `agents.py`, add `from .intake import IntakeSpec` and `from . import compute_routes` at the top. Replace `_executor` and `build_options`, and add the intake agent plus the handoff prompt. Change `_executor` signature to accept `mode`:

```python
_HANDOFF_EXEC_TOOLS = ["mcp__rep2struct__list_fold_jobs",
                       "mcp__rep2struct__build_fold_artifact"]


def _handoff_executor(name, tool):
    return AgentDefinition(
        description=f"{name}: builds the {tool} fold artifact for its group and stops "
                    f"for the user to run it.",
        prompt=(
            f"You prepare the {tool} folds for the jobs assigned to your group; you do NOT "
            f"run them. Call list_fold_jobs, and for each job whose tool is '{tool}' AND "
            f"whose done is false: call build_fold_artifact with tool='{tool}' and the run's "
            f"compute_route to write that job's self-contained artifact (a Colab notebook or "
            f"a local bash script, chosen by the route) and get its path. Report each "
            f"artifact path to the user with one line on how to run it: for a Colab notebook, "
            f"upload it and run the cells; for a bash script, run it on the target machine. "
            f"If build_fold_artifact reports route_wired false, say plainly the route runner "
            f"is not wired and hand over the script for the user to run. Never fabricate a "
            f"model or a score, and never open a browser."),
        tools=list(_HANDOFF_EXEC_TOOLS),
        model="sonnet",
    )


def _executor(name, tool, mode="auto"):
    if mode == "handoff":
        return _handoff_executor(name, tool)
    return AgentDefinition(
        description=f"{name}: folds the {tool} group by driving its Colab notebook through the browser.",
        prompt=(
            # ... UNCHANGED existing auto prompt (keep the current body verbatim) ...
        ),
        tools=list(_EXEC_TOOLS),
        model="sonnet",
    )
```

Keep the existing auto prompt body exactly as it is today (lines 13 to 43); only the signature gains `mode` and the handoff branch is new.

Add the intake agent inside `build_agents`, and thread `mode` through the executor loop:

```python
def build_agents(mode="auto"):
    agents = {
        "intake-agent": AgentDefinition(
            description="Conversational intake: asks data, question, and compute environment, "
                        "then hands a structured brief to the orchestrator.",
            prompt=(
                "You run the R2S intake interview. Ask the user, one question at a time and "
                "branching on their answers: (1) what kind of data they have, (2) where the "
                "input file is (a path, or a file they dropped into the run folder), (3) what "
                "question or task they want answered, (4) which compute environment they can "
                "run predictions on. For (4) call list_compute_routes and ask ONLY the "
                "required_fields of the route they name (for ssh: host, user, an SSH key "
                "first and a password only as a last resort, remote_path). If the user does "
                "not know, propose the default route (Colab), the simplest option. Never "
                "write a password anywhere. When you have all four, summarize the brief back "
                "and confirm before the run proceeds."),
            tools=["mcp__rep2struct__list_compute_routes",
                   "mcp__rep2struct__ingest_repertoire", "Agent"],
            model="opus",
        ),
        # ... the existing structure-strategist, qc-agent, report-agent UNCHANGED ...
    }
    for t in structure_tools.REGISTRY:
        agents[f"{t.name}-agent"] = _executor(f"{t.name}-agent", t.name, mode=mode)
    return agents
```

Update `build_options` to take `mode` and adjust `allowed_tools`:

```python
def build_options(run_dir, mode="auto"):
    base = [
        "Agent",
        "mcp__rep2struct__ingest_repertoire", "mcp__rep2struct__annotate_specificity",
        "mcp__rep2struct__prep_and_select", "mcp__rep2struct__list_structure_tools",
        "mcp__rep2struct__list_fold_jobs", "mcp__rep2struct__list_compute_routes",
        "mcp__rep2struct__record_fold_result", "mcp__rep2struct__record_local_folds",
        "mcp__rep2struct__qc_structure", "mcp__rep2struct__render_final_report",
    ]
    if mode == "handoff":
        allowed = base + ["mcp__rep2struct__build_fold_artifact"]
        servers = {"rep2struct": build_server()}
    else:
        allowed = base + ["mcp__rep2struct__build_fold_notebook", "mcp__playwright__*"]
        servers = {"rep2struct": build_server(),
                   "playwright": {"command": "npx", "args": ["@playwright/mcp@latest"]}}
    return ClaudeAgentOptions(
        mcp_servers=servers,
        agents=build_agents(mode=mode),
        allowed_tools=allowed,
        permission_mode="bypassPermissions",
    )
```

Add the handoff orchestrator prompt:

```python
def intake_orchestrator_prompt(run_dir, spec):
    return (
        f"Run the repertoire to structure pipeline on {spec.input_path} with run_dir "
        f"{run_dir}.\nUser question steering the routing: {spec.question}\n"
        f"Compute route: {spec.compute_route} (route params {spec.route_params}).\n"
        f"1. Call ingest_repertoire, then annotate_specificity (honest annotation, keep "
        f"unannotatable as is).\n"
        f"2. Call prep_and_select with top_n 8.\n"
        f"3. Delegate to the structure-strategist to route each group to a tool; each tool's "
        f"executor builds the fold artifact for compute_route '{spec.compute_route}' and "
        f"STOPS. Present the plan (which group to which tool, which artifact files) and let "
        f"the user confirm before generating artifacts.\n"
        f"4. Tell the user to run the artifacts and rerun R2S on this run_dir to resume.")
```

- [ ] **Step 4: Run test to verify it passes**

Run: `./.venv/bin/python -m pytest tests/test_agents_intake.py tests/test_agents_config.py -q`
Expected: PASS (existing `test_agents_config.py` still green, new file passes)

- [ ] **Step 5: Commit**

```bash
git add src/rep2struct/agents.py tests/test_agents_intake.py
git commit -m "feat(intake): intake-agent + handoff executors + build_options(mode)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 8: interactive CLI entry point (`python -m rep2struct`)

**Files:**
- Create: `src/rep2struct/cli.py`
- Create: `src/rep2struct/__main__.py`
- Test: `tests/test_cli.py`

**Interfaces:**
- Consumes: `intake.next_phase`, `intake.load_intake`, `agents.build_options`, `agents.intake_orchestrator_prompt`, `claude_agent_sdk.ClaudeSDKClient`.
- Produces: `plan_from_run_dir(run_dir: str) -> str` (pure: returns `"intake"` or `"run"` via `next_phase`, the seam the test exercises); `async def run_session(run_dir: str)` (the interactive loop, validated manually, not unit tested); `main()` reading `sys.argv[1]` as run_dir.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_cli.py
from rep2struct import cli, intake
from rep2struct.intake import IntakeSpec


def test_plan_selects_intake_then_run(tmp_path):
    assert cli.plan_from_run_dir(str(tmp_path)) == "intake"
    intake.save_intake(str(tmp_path), IntakeSpec("d", "/i.csv", "q", "colab", {}))
    assert cli.plan_from_run_dir(str(tmp_path)) == "run"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./.venv/bin/python -m pytest tests/test_cli.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'rep2struct.cli'`

- [ ] **Step 3: Write minimal implementation**

```python
# src/rep2struct/cli.py
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


async def run_session(run_dir: str) -> None:
    """Phase 0/A/B driver. In the intake phase the intake-agent interviews the user turn by
    turn (stdin), then the orchestrator runs in handoff mode and stops at the artifacts. A
    rerun on the same run_dir enters the run phase and resumes through the checkpoint."""
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
    prompt = (intake_orchestrator_prompt(run_dir, spec) if spec
              else orchestrator_prompt("", run_dir, 8))
    async with ClaudeSDKClient(options=opts) as client:
        await client.query(
            prompt + "\nIf fold artifacts have already been run, call record_local_folds "
                     "first, then proceed to QC and the report.")
        await _drain(client)


def main():
    import asyncio
    run_dir = sys.argv[1] if len(sys.argv) > 1 else "runs/session"
    asyncio.run(run_session(run_dir))


if __name__ == "__main__":
    main()
```

```python
# src/rep2struct/__main__.py
from .cli import main

if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run test to verify it passes**

Run: `./.venv/bin/python -m pytest tests/test_cli.py -q`
Expected: PASS (1 passed)

- [ ] **Step 5: Run the full suite to confirm no regressions**

Run: `./.venv/bin/python -m pytest -q`
Expected: PASS (all prior tests plus the new ones green)

- [ ] **Step 6: Commit**

```bash
git add src/rep2struct/cli.py src/rep2struct/__main__.py tests/test_cli.py
git commit -m "feat(intake): interactive ClaudeSDKClient CLI (python -m rep2struct)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Manual validation (after the suite is green)

The interview and the handoff stop cannot be unit tested (they need a live SDK session). Validate manually once, no live fold required:

1. `./.venv/bin/python -m rep2struct runs/intake_smoke` and answer the four questions (data type, a real 10x CSV path from `Data/dataset`, a question, choose `local_gpu` with a working path). Confirm the agent asks branched questions, never asks for a Colab detail on the local_gpu route, and stops after writing `runs/intake_smoke/scripts/*.sh` and `runs/intake_smoke/intake.json` with no secret in the json.
2. Rerun `./.venv/bin/python -m rep2struct runs/intake_smoke`. Confirm it enters the run phase (does not re-interview) and, with no CIFs present, reports the fold artifacts as pending rather than fabricating results.

## Self-review notes (author)

- Spec section 1 (architecture): Tasks 7 and 8 build the phases and the handoff mode. Covered.
- Spec section 2 (compute_routes registry): Task 1. Covered.
- Spec section 3 (IntakeSpec + intake tools): Tasks 2 and 6. `recommend_route` resolved to a tested `compute_routes.recommend` helper plus agent-prompt reasoning (open decision 2), not a separate MCP tool, to keep the tool surface minimal.
- Spec section 4 (handoff artifact dispatch + propose then confirm): Tasks 3, 4, 7. Covered.
- Spec section 5 (resume by checkpoint): Task 5 (`record_local_folds`) plus the run-phase prompt in Task 8. Covered.
- Testing section: every deterministic layer has an offline test; the interview is left to manual validation, as the spec states.
- Open decision 1 (entry point) resolved to `python -m rep2struct` (no console-script packaging change). Open decision 3 (resume interview) resolved to the same session via `next_phase`, no separate agent.
