# Structure Strategist Routing Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the hardwired Protenix `fold-agent` with a `structure-strategist` agent that reasons over a data registry of structure tools, partitions clonotypes into homogeneous groups, and delegates each group to a per-tool executor, with per-group QC calibration, honest reporting, and MSA restored as a pre-fold artifact.

**Architecture:** A pure-data `StructureTool` registry is the single source of truth for tool validity domains. Deterministic helpers expose facts (matching tools, default, coverage) to a reasoning strategist agent; the strategist picks one tool per homogeneous group and delegates to a small per-tool executor agent. QC calibrates per group and branches on the tool's output type; the report encodes honesty guards in its schema. MSA is computed outside the fold runtime and embedded, killing the remote-MSA-throttle failure.

**Tech Stack:** Python 3.11, dataclasses, claude-agent-sdk 0.2.113, Playwright MCP, pytest. No new heavy deps in the deterministic layer; mmseqs2 is an external binary invoked only in the MSA prep path.

## Global Constraints

- Python 3.11; run tests with `./.venv/bin/python -m pytest -q` from the repo root.
- Protenix is the default tool; NEVER Boltz.
- Docs and report copy obey the no-dash-as-punctuation rule (use commas, periods, parentheses).
- This repo KEEPS the `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>` commit trailer.
- A separate session currently drives Protenix; do NOT modify the Protenix Colab fold code path or `hla_ectodomains.json`.
- Honesty rules are enforced in SCHEMAS, not just prose: a `binding_score` result is never labelled "fold"/"structure"; groups are never cross-compared on raw distances; a clean fold never confirms specificity.
- Full pytest suite must stay green (currently 33/33).

---

### Task 1: StructureTool registry and helpers

**Files:**
- Create: `src/rep2struct/structure_tools.py`
- Test: `tests/test_structure_tools.py`

**Interfaces:**
- Consumes: nothing.
- Produces:
  - `StructureTool` dataclass: `name: str`, `validity: dict`, `output_type: str`, `strengths: str`, `limits: str`, `colab_adapter: str`, `is_default: bool = False`. `validity` has keys `mhc_class: set[int]`, `needs_tcr: Optional[bool]`, `species: str`.
  - `REGISTRY: list[StructureTool]`.
  - `get_default() -> StructureTool`.
  - `tools_for(mhc_class: int, has_tcr: bool, species: str, output_needed: str) -> list[StructureTool]`.
  - `is_covered(mhc_class: int, has_tcr: bool, species: str, output_needed: str) -> bool`.
  - `as_dicts() -> list[dict]` (JSON-safe, sets become sorted lists) for the MCP tool.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_structure_tools.py
from rep2struct import structure_tools as st


def test_default_is_protenix():
    d = st.get_default()
    assert d.name == "protenix" and d.is_default


def test_registry_has_the_five_v1_tools():
    names = {t.name for t in st.REGISTRY}
    assert names == {"protenix", "af3", "mhcfine", "tcrdock", "affinetune"}


def test_affinetune_is_binding_score_the_rest_structure():
    by = {t.name: t for t in st.REGISTRY}
    assert by["affinetune"].output_type == "binding_score"
    assert by["protenix"].output_type == "structure"


def test_tools_for_class_ii_structure_returns_protenix_not_mhcfine():
    got = {t.name for t in st.tools_for(2, has_tcr=True, species="human", output_needed="structure")}
    assert "protenix" in got and "mhcfine" not in got  # mhcfine is class I only


def test_tools_for_binding_score_class_ii_returns_affinetune():
    got = {t.name for t in st.tools_for(2, has_tcr=False, species="mouse", output_needed="binding_score")}
    assert got == {"affinetune"}


def test_is_covered_true_for_structure_false_when_no_match():
    assert st.is_covered(1, True, "human", "structure") is True
    # binding_score for a species affinetune does not list -> not covered
    assert st.is_covered(1, False, "reptile", "binding_score") is False


def test_as_dicts_is_json_safe():
    import json
    json.dumps(st.as_dicts())  # must not raise
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./.venv/bin/python -m pytest tests/test_structure_tools.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'rep2struct.structure_tools'`

- [ ] **Step 3: Write minimal implementation**

```python
# src/rep2struct/structure_tools.py
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional


@dataclass(frozen=True)
class StructureTool:
    name: str
    validity: dict            # {"mhc_class": set[int], "needs_tcr": Optional[bool], "species": str}
    output_type: str          # "structure" | "binding_score"
    strengths: str
    limits: str
    colab_adapter: str
    is_default: bool = False


REGISTRY: list[StructureTool] = [
    StructureTool(
        name="protenix",
        validity={"mhc_class": {1, 2}, "needs_tcr": None, "species": "any"},
        output_type="structure",
        strengths="full 3-chain TCR-pMHC fold, default workhorse",
        limits="imposes canonical geometry even on non-binders (basis of skeptical QC)",
        colab_adapter="protenix_colab",
        is_default=True,
    ),
    StructureTool(
        name="af3",
        validity={"mhc_class": {1, 2}, "needs_tcr": None, "species": "any"},
        output_type="structure",
        strengths="AlphaFold3-class accuracy when weights are available",
        limits="gated model weights; only if the user has them",
        colab_adapter="af3_colab",
    ),
    StructureTool(
        name="mhcfine",
        validity={"mhc_class": {1}, "needs_tcr": False, "species": "any"},
        output_type="structure",
        strengths="most precise class I peptide pose (RMSD 0.66A)",
        limits="class I only, no TCR",
        colab_adapter="mhcfine_colab",
    ),
    StructureTool(
        name="tcrdock",
        validity={"mhc_class": {1, 2}, "needs_tcr": True, "species": "any"},
        output_type="structure",
        strengths="TCR:pMHC interface and V-domain anchoring",
        limits="template-coverage limited; class II not systematically benchmarked",
        colab_adapter="tcrdock_colab",
    ),
    StructureTool(
        name="affinetune",
        validity={"mhc_class": {1, 2}, "needs_tcr": False, "species": "any"},
        output_type="binding_score",
        strengths="is-this-peptide-presented classifier, class I and II",
        limits="returns a presentation score, not a structure",
        colab_adapter="affinetune_colab",
    ),
]


def get_default() -> StructureTool:
    return next(t for t in REGISTRY if t.is_default)


def _species_ok(tool: StructureTool, species: str) -> bool:
    return tool.validity["species"] == "any" or tool.validity["species"] == species


def _tcr_ok(tool: StructureTool, has_tcr: bool) -> bool:
    need = tool.validity["needs_tcr"]
    return need is None or need == has_tcr


def tools_for(mhc_class: int, has_tcr: bool, species: str, output_needed: str) -> list[StructureTool]:
    return [
        t for t in REGISTRY
        if t.output_type == output_needed
        and mhc_class in t.validity["mhc_class"]
        and _tcr_ok(t, has_tcr)
        and _species_ok(t, species)
    ]


def is_covered(mhc_class: int, has_tcr: bool, species: str, output_needed: str) -> bool:
    return len(tools_for(mhc_class, has_tcr, species, output_needed)) > 0


def as_dicts() -> list[dict]:
    out = []
    for t in REGISTRY:
        v = dict(t.validity)
        v["mhc_class"] = sorted(v["mhc_class"])
        out.append({
            "name": t.name, "validity": v, "output_type": t.output_type,
            "strengths": t.strengths, "limits": t.limits, "is_default": t.is_default,
        })
    return out
```

- [ ] **Step 4: Run test to verify it passes**

Run: `./.venv/bin/python -m pytest tests/test_structure_tools.py -q`
Expected: PASS (7 passed)

- [ ] **Step 5: Commit**

```bash
git add src/rep2struct/structure_tools.py tests/test_structure_tools.py
git commit -m "feat: StructureTool registry with validity-domain helpers

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: FoldJob tags and homogeneous grouping

**Files:**
- Modify: `src/rep2struct/schema.py:29-35` (FoldJob)
- Create: `src/rep2struct/grouping.py`
- Test: `tests/test_grouping.py`

**Interfaces:**
- Consumes: `FoldJob` from Task's schema change.
- Produces:
  - `FoldJob` gains fields: `mhc_class: int = 1`, `has_tcr: bool = True`, `species: str = "human"`, `output_needed: str = "structure"`, `tool: Optional[str] = None`, `group_id: Optional[str] = None`. (Existing fields, `msa_ref`, keep their positions; add the new ones after `model_paths`.)
  - `grouping.group_key(job: FoldJob) -> str` returns `f"c{mhc_class}_{'tcr' if has_tcr else 'notcr'}_{species}_{output_needed}"`.
  - `grouping.partition(jobs: list[FoldJob]) -> dict[str, list[FoldJob]]` groups by `group_key`, stamping each job's `group_id`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_grouping.py
from rep2struct.schema import FoldJob
from rep2struct import grouping


def _job(cid, mhc_class=1, has_tcr=True, species="human", output_needed="structure"):
    return FoldJob(clonotype_id=cid, construct_fasta=">A\nAAAA",
                   mhc_class=mhc_class, has_tcr=has_tcr, species=species,
                   output_needed=output_needed)


def test_group_key_is_stable_and_descriptive():
    assert grouping.group_key(_job("x")) == "c1_tcr_human_structure"
    assert grouping.group_key(_job("y", mhc_class=2, has_tcr=False,
                                    species="mouse", output_needed="binding_score")) \
        == "c2_notcr_mouse_binding_score"


def test_partition_splits_and_stamps_group_id():
    jobs = [_job("a"), _job("b"), _job("c", mhc_class=2)]
    groups = grouping.partition(jobs)
    assert set(groups) == {"c1_tcr_human_structure", "c2_tcr_human_structure"}
    assert len(groups["c1_tcr_human_structure"]) == 2
    assert all(j.group_id == "c1_tcr_human_structure" for j in groups["c1_tcr_human_structure"])
    assert groups["c2_tcr_human_structure"][0].group_id == "c2_tcr_human_structure"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./.venv/bin/python -m pytest tests/test_grouping.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'rep2struct.grouping'`

- [ ] **Step 3: Write minimal implementation**

Modify `src/rep2struct/schema.py`, replacing the FoldJob dataclass body:

```python
@dataclass
class FoldJob:
    clonotype_id: str
    construct_fasta: str         # A..E chains
    msa_ref: Optional[str] = None
    status: str = "pending"      # pending, done, failed
    model_paths: list[str] = field(default_factory=list)
    mhc_class: int = 1           # 1 or 2
    has_tcr: bool = True
    species: str = "human"
    output_needed: str = "structure"   # structure | binding_score
    tool: Optional[str] = None   # tool the strategist assigned
    group_id: Optional[str] = None
```

Create `src/rep2struct/grouping.py`:

```python
from __future__ import annotations
from .schema import FoldJob


def group_key(job: FoldJob) -> str:
    tcr = "tcr" if job.has_tcr else "notcr"
    return f"c{job.mhc_class}_{tcr}_{job.species}_{job.output_needed}"


def partition(jobs: list[FoldJob]) -> dict[str, list[FoldJob]]:
    groups: dict[str, list[FoldJob]] = {}
    for j in jobs:
        k = group_key(j)
        j.group_id = k
        groups.setdefault(k, []).append(j)
    return groups
```

- [ ] **Step 4: Run test to verify it passes**

Run: `./.venv/bin/python -m pytest tests/test_grouping.py tests/test_foldprep.py tests/test_fold.py -q`
Expected: PASS (new grouping tests pass; existing foldprep/fold tests still pass because new FoldJob fields have defaults)

- [ ] **Step 5: Commit**

```bash
git add src/rep2struct/schema.py src/rep2struct/grouping.py tests/test_grouping.py
git commit -m "feat: tag FoldJob with routing metadata and partition into homogeneous groups

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: MSA prep artifact with fallback chain

**Files:**
- Create: `src/rep2struct/msa.py`
- Test: `tests/test_msa.py`

**Interfaces:**
- Consumes: `FoldJob`.
- Produces:
  - `build_msa(job: FoldJob, run_dir, local_runner=None, colab_runner=None) -> tuple[str, str]` returns `(msa_ref, basis)` where `basis` is `"local"`, `"colab_cpu"`, or `"none"`. `msa_ref` is a path string, or `""` when basis is `"none"`.
  - A runner is a callable `runner(fasta: str) -> str` that returns a3m text, or raises to signal failure. `local_runner` is tried first, then `colab_runner`, then MSA-free fallback.

**Interface note:** runners are injected exactly like `agent_tools.configure(sim_fn=..., assign_fn=...)` injects offline fakes, so this task needs no mmseqs2 binary to be tested.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_msa.py
from rep2struct.schema import FoldJob
from rep2struct import msa


def _job():
    return FoldJob(clonotype_id="c1", construct_fasta=">A\nAAAA")


def test_local_runner_wins_and_caches_a3m(tmp_path):
    ref, basis = msa.build_msa(_job(), tmp_path, local_runner=lambda f: "A3M-LOCAL")
    assert basis == "local"
    assert open(ref).read() == "A3M-LOCAL"


def test_falls_back_to_colab_when_local_fails(tmp_path):
    def boom(f): raise RuntimeError("no local DB")
    ref, basis = msa.build_msa(_job(), tmp_path, local_runner=boom,
                               colab_runner=lambda f: "A3M-COLAB")
    assert basis == "colab_cpu"
    assert open(ref).read() == "A3M-COLAB"


def test_falls_back_to_msa_free_when_both_fail(tmp_path):
    def boom(f): raise RuntimeError("down")
    ref, basis = msa.build_msa(_job(), tmp_path, local_runner=boom, colab_runner=boom)
    assert basis == "none" and ref == ""
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./.venv/bin/python -m pytest tests/test_msa.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'rep2struct.msa'`

- [ ] **Step 3: Write minimal implementation**

```python
# src/rep2struct/msa.py
from __future__ import annotations
from pathlib import Path
from .schema import FoldJob


def build_msa(job: FoldJob, run_dir, local_runner=None, colab_runner=None) -> tuple[str, str]:
    """Compute an MSA artifact OUTSIDE the fold runtime and cache it.

    Tries local mmseqs2 first, then a Colab CPU step, then falls back to
    MSA-free. Returns (msa_ref, basis). Removing the MSA from the fold
    runtime is what kills the remote-MSA-server throttle failure.
    """
    for runner, basis in ((local_runner, "local"), (colab_runner, "colab_cpu")):
        if runner is None:
            continue
        try:
            a3m = runner(job.construct_fasta)
        except Exception:  # noqa: BLE001  -- any runner failure degrades to the next path
            continue
        out = Path(run_dir) / f"msa_{job.clonotype_id}.a3m"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(a3m)
        return str(out), basis
    return "", "none"
```

- [ ] **Step 4: Run test to verify it passes**

Run: `./.venv/bin/python -m pytest tests/test_msa.py -q`
Expected: PASS (3 passed)

- [ ] **Step 5: Commit**

```bash
git add src/rep2struct/msa.py tests/test_msa.py
git commit -m "feat: MSA prep artifact with local -> colab-cpu -> msa-free fallback

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: list_structure_tools MCP tool and tag emission in prep_and_select

**Files:**
- Modify: `src/rep2struct/agent_tools.py` (imports near line 8; `prep_and_select` at 59-79; add a new `@tool`; register it in `build_server` at 138-142)
- Test: `tests/test_agent_tools.py` (append)

**Interfaces:**
- Consumes: `structure_tools.as_dicts`, `grouping.partition` (Tasks 1, 2), `msa.build_msa` (Task 3).
- Produces:
  - New MCP tool `list_structure_tools(run_dir)` whose `structuredContent` is `{"tools": structure_tools.as_dicts()}`.
  - `prep_and_select` now calls `grouping.partition(jobs)` so each job carries a `group_id`, calls `msa.build_msa(job, run_dir)` per job to stamp `job.msa_ref` (with no runners injected this returns `""` = MSA-free, reproducing today's behaviour; production injects the mmseqs2 runners), then persists the jobs. The v1 class I / structure tags come from the FoldJob dataclass defaults.

**Interface note:** v1 class II construct geometry is out of scope (see spec); prep still emits class I structure jobs, but they now flow through partitioning so the strategist and QC see `group_id`. The tags default to class I / has_tcr / human / structure, which reproduces today's behaviour.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_agent_tools.py  (append)
import asyncio
from rep2struct import agent_tools as at
from rep2struct.runstate import RunState
from rep2struct.schema import Clonotype, Annotation


def _run(coro):
    return asyncio.get_event_loop().run_until_complete(coro)


def test_list_structure_tools_returns_registry():
    res = _run(at.list_structure_tools.handler({"run_dir": "/tmp/whatever"}))
    names = {t["name"] for t in res["structuredContent"]["tools"]}
    assert names == {"protenix", "af3", "mhcfine", "tcrdock", "affinetune"}


def test_prep_and_select_stamps_group_id(tmp_path):
    rd = str(tmp_path / "run")
    clon = Clonotype(id="c1", trav="TRAV1", cdr3a="CAA", trbv="TRBV2", cdr3b="CAB",
                     size=5, traj="TRAJ1", trbj="TRBJ1")
    ann = Annotation(clonotype_id="c1", annotatable=True, confidence_tier="high",
                     epitope="SIINFEKL", hla="A*02:01")
    RunState(rd).write_stage("ingest", [clon])
    RunState(rd).write_stage("annotate", [ann])
    at.configure(assign_fn=lambda c: c)  # no allele network call
    try:
        _run(at.prep_and_select.handler({"run_dir": rd, "top_n": 5}))
    finally:
        at.configure()
    jobs = RunState(rd).read_stage("foldjobs")
    assert jobs and all(j["group_id"] == "c1_tcr_human_structure" for j in jobs)
```

**Note for the implementer:** `build_tcr_seqs` / `build_mhc_seqs` reach the network for real sequences. If this test cannot resolve sequences offline in your environment, wrap the two calls in `prep_and_select` so a monkeypatched fake can be injected the same way `_CFG` injects `sim_fn`/`assign_fn`, and inject fakes in the test. Keep the group-id assertion regardless.

- [ ] **Step 2: Run test to verify it fails**

Run: `./.venv/bin/python -m pytest tests/test_agent_tools.py -q -k "structure_tools or group_id"`
Expected: FAIL with `AttributeError: module 'rep2struct.agent_tools' has no attribute 'list_structure_tools'`

- [ ] **Step 3: Write minimal implementation**

In `agent_tools.py`, add to the imports block:

```python
from . import structure_tools
from .grouping import partition
from .msa import build_msa
```

Add the new tool (near `list_fold_jobs`):

```python
@tool("list_structure_tools", "List the structure tools and their validity domains for the strategist.",
      {"run_dir": str})
async def list_structure_tools(args):
    r = _txt("structure tool registry")
    r["structuredContent"] = {"tools": structure_tools.as_dicts()}
    return r
```

In `prep_and_select`, replace the persist line
`RunState(args["run_dir"]).write_stage("foldjobs", jobs)` with:

```python
    partition(jobs)  # stamps group_id on each job in place
    for j in jobs:   # MSA is a pre-fold artifact; no runners here = MSA-free default
        j.msa_ref, _ = build_msa(j, args["run_dir"])
    RunState(args["run_dir"]).write_stage("foldjobs", jobs)
```

The `test_prep_and_select_stamps_group_id` test already exercises this path;
with no runners injected, `build_msa` returns `""` so `msa_ref` stays empty and
the assertion on `group_id` is unaffected.

Register the tool in `build_server`:

```python
    return create_sdk_mcp_server(name="rep2struct", version="0.1.0", tools=[
        ingest_repertoire, annotate_specificity, prep_and_select, list_fold_jobs,
        list_structure_tools, record_fold_result, qc_structure, render_final_report,
    ])
```

- [ ] **Step 4: Run test to verify it passes**

Run: `./.venv/bin/python -m pytest tests/test_agent_tools.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/rep2struct/agent_tools.py tests/test_agent_tools.py
git commit -m "feat: list_structure_tools tool and group_id stamping in prep_and_select

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 5: QCResult provenance fields and binding-score verdict path

**Files:**
- Modify: `src/rep2struct/schema.py:37-44` (QCResult)
- Modify: `src/rep2struct/qc.py` (add `verdict_binding`)
- Test: `tests/test_qc.py` (append)

**Interfaces:**
- Consumes: `QCResult`.
- Produces:
  - `QCResult` gains `tool: Optional[str] = None` and `calibration_basis: Optional[str] = None`.
  - `qc.verdict_binding(score: float, threshold: float, clonotype_id: str, tool: str) -> QCResult`: returns verdict `"presented"` when `score > threshold` else `"not_presented"`, reason text that says "predicted presentation", never "fold"/"structure". Sets `tool` and `calibration_basis="binding_score_null"`. Leaves geometric fields (`cdr3_pep_atoms`, `dockq`, `crossing_angle`) None.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_qc.py  (append)
from rep2struct.qc import verdict_binding


def test_binding_verdict_presented_and_not_presented():
    hi = verdict_binding(0.9, 0.5, "c1", tool="affinetune")
    lo = verdict_binding(0.3, 0.5, "c2", tool="affinetune")
    assert hi.qc_verdict == "presented" and lo.qc_verdict == "not_presented"


def test_binding_verdict_is_honest_not_a_fold():
    r = verdict_binding(0.9, 0.5, "c1", tool="affinetune")
    assert "presentation" in r.reason.lower()
    assert "fold" not in r.reason.lower() and "structure" not in r.reason.lower()
    assert r.tool == "affinetune" and r.calibration_basis == "binding_score_null"
    assert r.cdr3_pep_atoms is None and r.dockq is None
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./.venv/bin/python -m pytest tests/test_qc.py -q -k binding`
Expected: FAIL with `ImportError: cannot import name 'verdict_binding'`

- [ ] **Step 3: Write minimal implementation**

In `schema.py`, replace the QCResult dataclass body:

```python
@dataclass
class QCResult:
    clonotype_id: str
    qc_verdict: str              # reliable, suspect, qc_failed, presented, not_presented
    reason: str
    dockq: Optional[float] = None
    cdr3_pep_atoms: Optional[float] = None
    crossing_angle: Optional[float] = None
    tool: Optional[str] = None
    calibration_basis: Optional[str] = None
```

In `qc.py`, add:

```python
def verdict_binding(score: float, threshold: float, clonotype_id: str, tool: str) -> QCResult:
    presented = score > threshold
    return QCResult(
        clonotype_id,
        "presented" if presented else "not_presented",
        ("predicted presentation above the score null" if presented
         else "predicted presentation not above the score null"),
        tool=tool,
        calibration_basis="binding_score_null",
    )
```

- [ ] **Step 4: Run test to verify it passes**

Run: `./.venv/bin/python -m pytest tests/test_qc.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/rep2struct/schema.py src/rep2struct/qc.py tests/test_qc.py
git commit -m "feat: QCResult provenance fields and honest binding-score verdict path

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 6: Record the tool on fold results and make qc_structure per-group and output-aware

**Files:**
- Modify: `src/rep2struct/agent_tools.py` (`record_fold_result` 92-99; `qc_structure` 102-120)
- Test: `tests/test_agent_tools.py` (append)

**Interfaces:**
- Consumes: `qc.verdict_binding` (Task 5), `structure_tools` (Task 1), `qc.score_model`/`qc.verdict`.
- Produces:
  - `record_fold_result(run_dir, clonotype_id, model_paths, tool="protenix")` persists `{clonotype_id: {"paths": [...], "tool": tool}}`. Reading tolerates the old list shape by treating a bare list as `{"paths": list, "tool": "protenix"}`.
  - `qc_structure(run_dir, clonotype_id, scramble_threshold, output_type="structure", tool="protenix")`: when `output_type == "binding_score"`, treat the first recorded path as a float score file and call `verdict_binding`; otherwise the existing structure path. The stored/returned `QCResult` carries `tool`.

**Interface note:** the `scramble_threshold` is supplied PER CALL by the qc-agent from that group's own calibration; there is no global threshold. The binding-score path reads a single float from the recorded artifact (the executor writes the score there).

- [ ] **Step 1: Write the failing test**

```python
# tests/test_agent_tools.py  (append)
def test_record_fold_result_keeps_tool_and_is_back_compatible(tmp_path):
    rd = str(tmp_path / "run")
    _run(at.record_fold_result.handler(
        {"run_dir": rd, "clonotype_id": "c1", "model_paths": ["c1.cif"], "tool": "tcrdock"}))
    done = RunState(rd).read_stage("folds")
    assert done["c1"]["tool"] == "tcrdock" and done["c1"]["paths"] == ["c1.cif"]


def test_qc_structure_binding_path_uses_binding_verdict(tmp_path):
    rd = str(tmp_path / "run")
    score_file = tmp_path / "c1.score"
    score_file.write_text("0.9")
    _run(at.record_fold_result.handler(
        {"run_dir": rd, "clonotype_id": "c1", "model_paths": [str(score_file)],
         "tool": "affinetune"}))
    res = _run(at.qc_structure.handler(
        {"run_dir": rd, "clonotype_id": "c1", "scramble_threshold": 0.5,
         "output_type": "binding_score", "tool": "affinetune"}))
    assert res["structuredContent"]["qc_verdict"] == "presented"
    stored = RunState(rd).read_stage("qc")
    assert stored[0]["tool"] == "affinetune"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./.venv/bin/python -m pytest tests/test_agent_tools.py -q -k "tool or binding_path"`
Expected: FAIL (record_fold_result ignores `tool`; qc_structure has no binding path)

- [ ] **Step 3: Write minimal implementation**

Replace `record_fold_result`:

```python
@tool("record_fold_result", "Record the model paths a fold produced for one clonotype, with the tool used.",
      {"run_dir": str, "clonotype_id": str, "model_paths": list, "tool": str})
async def record_fold_result(args):
    rs = RunState(args["run_dir"])
    done = rs.read_stage("folds") if rs.stage_done("folds") else {}
    done[args["clonotype_id"]] = {"paths": args["model_paths"],
                                  "tool": args.get("tool", "protenix")}
    rs.write_stage("folds", done)
    return _txt(f"recorded {len(args['model_paths'])} models for {args['clonotype_id']} via {args.get('tool', 'protenix')}")
```

Replace `qc_structure` with an output-aware version:

```python
@tool("qc_structure", "Score a fold (per-group threshold) and return a skeptical verdict; output-type aware.",
      {"run_dir": str, "clonotype_id": str, "scramble_threshold": float,
       "output_type": str, "tool": str})
async def qc_structure(args):
    from .qc import verdict_binding
    rs = RunState(args["run_dir"])
    done = rs.read_stage("folds") if rs.stage_done("folds") else {}
    rec = done.get(args["clonotype_id"], {})
    if isinstance(rec, list):                      # back-compat with the old list shape
        rec = {"paths": rec, "tool": "protenix"}
    paths = rec.get("paths", [])
    tool = args.get("tool", rec.get("tool", "protenix"))
    output_type = args.get("output_type", "structure")
    if not paths:
        res = QCResult(args["clonotype_id"], "qc_failed", "no model recorded", tool=tool)
    elif output_type == "binding_score":
        score = float(Path(paths[0]).read_text().strip())
        res = verdict_binding(score, args["scramble_threshold"], args["clonotype_id"], tool=tool)
    else:
        s = score_model(paths[0])
        s["clonotype_id"] = args["clonotype_id"]
        res = verdict(s, args["scramble_threshold"])
        res.tool = tool
        res.calibration_basis = "scramble_null"
    qcs = rs.read_stage("qc") if rs.stage_done("qc") else []
    qcs = [q for q in qcs if q["clonotype_id"] != res.clonotype_id]
    qcs.append(asdict(res))
    rs.write_stage("qc", qcs)
    r = _txt(f"{res.clonotype_id}: {res.qc_verdict} ({res.reason})")
    r["structuredContent"] = {"qc_verdict": res.qc_verdict, "reason": res.reason}
    return r
```

- [ ] **Step 4: Run test to verify it passes**

Run: `./.venv/bin/python -m pytest tests/test_agent_tools.py -q`
Expected: PASS (existing structure QC test still passes via back-compat; new tests pass)

- [ ] **Step 5: Commit**

```bash
git add src/rep2struct/agent_tools.py tests/test_agent_tools.py
git commit -m "feat: record fold tool and make qc_structure per-group and output-aware

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 7: Report honesty guards

**Files:**
- Modify: `src/rep2struct/report.py` (`render_report` signature and template data)
- Test: `tests/test_report.py` (append)

**Interfaces:**
- Consumes: `QCResult` with `tool`, `qc_verdict` in `{reliable, suspect, presented, not_presented, qc_failed}`.
- Produces: `render_report(clonotypes, annotations, qcs)` HTML that:
  - labels any `presented`/`not_presented` row as "predicted presentation", never "fold" or "structure";
  - shows each row's `tool`;
  - contains no sort of rows by `cdr3_pep_atoms` (rows stay in input order so raw distances are never presented as a ranking).

**Interface note:** read `report.py` first and follow its existing jinja2 template structure. Add a per-row `tool` column and an evidence-type label derived from `qc_verdict`; do not restructure the template.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_report.py  (append)
from rep2struct.report import render_report
from rep2struct.schema import Clonotype, Annotation, QCResult


def _fixtures():
    c = Clonotype(id="c1", trav="TRAV1", cdr3a="CAA", trbv="TRBV2", cdr3b="CAB", size=5)
    a = Annotation(clonotype_id="c1", annotatable=True, confidence_tier="high",
                   epitope="SIINFEKL", hla="A*02:01")
    return [c], [a]


def test_binding_row_is_labelled_predicted_presentation():
    c, a = _fixtures()
    q = QCResult("c1", "presented", "predicted presentation above the score null",
                 tool="affinetune", calibration_basis="binding_score_null")
    html = render_report(c, a, [q])
    assert "predicted presentation" in html.lower()
    assert "affinetune" in html


def test_binding_row_not_called_a_structure_or_fold():
    c, a = _fixtures()
    q = QCResult("c1", "presented", "predicted presentation above the score null",
                 tool="affinetune")
    html = render_report(c, a, [q])
    row = [ln for ln in html.splitlines() if "affinetune" in ln]
    assert row and all("fold" not in ln.lower() and "structure" not in ln.lower() for ln in row)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./.venv/bin/python -m pytest tests/test_report.py -q -k "binding or presentation"`
Expected: FAIL (report has no tool column / presentation label)

- [ ] **Step 3: Write minimal implementation**

Read `report.py`, then in the row-building logic add, per QCResult:

```python
    EVIDENCE = {
        "reliable": "structure (reliable)",
        "suspect": "structure (suspect)",
        "qc_failed": "structure (qc failed)",
        "presented": "predicted presentation",
        "not_presented": "predicted presentation",
    }
    # for each qc result q, expose to the template:
    #   q.tool  and  evidence = EVIDENCE.get(q.qc_verdict, "structure")
```

Add a `Tool` and an `Evidence` cell to the results table template (jinja2), binding `evidence` and `q.tool`. Do not sort rows by any numeric field; keep the clonotype input order.

- [ ] **Step 4: Run test to verify it passes**

Run: `./.venv/bin/python -m pytest tests/test_report.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/rep2struct/report.py tests/test_report.py
git commit -m "feat: report labels binding-score rows as predicted presentation, shows tool, no raw ranking

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 8: Wire the structure-strategist and per-tool executor agents

**Files:**
- Modify: `src/rep2struct/agents.py` (`build_agents` 6-37; `build_options` 40-56; `orchestrator_prompt` 59-66)
- Test: `tests/test_agents_config.py` (replace assertions)

**Interfaces:**
- Consumes: everything above; `structure_tools` names for the executor set.
- Produces:
  - `build_agents()` returns, in addition to `qc-agent` and `report-agent`: a `structure-strategist` agent and one executor per tool: `protenix-agent`, `af3-agent`, `mhcfine-agent`, `tcrdock-agent`, `affinetune-agent`.
  - The strategist has tools `mcp__rep2struct__list_structure_tools`, `mcp__rep2struct__list_fold_jobs`, and `Agent` (to delegate); model `opus`.
  - Each executor has `mcp__rep2struct__list_fold_jobs`, `mcp__rep2struct__record_fold_result`, `mcp__playwright__*`; model `sonnet`.
  - `build_options` lists the new tool `mcp__rep2struct__list_structure_tools` and all agents in `allowed_tools`/`agents`.
  - `orchestrator_prompt(csv_path, run_dir, top_n, question=None)` gains an optional `question`; when present it is passed to the strategist as the steering question; when absent the strategist routes from construct tags. The prompt still names ingest/annotate/fold/qc/report stages.

**Interface note:** only `protenix-agent` has a working Colab notebook today. The others are wired so the strategist can select and delegate, but their executor prompt must state honestly that if their Colab adapter is not yet available they must report the job as not-run rather than fabricate a result. This preserves the "no fabricated evidence" rule while the notebooks are built.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_agents_config.py  (replace file contents)
from rep2struct.agents import build_agents, build_options, orchestrator_prompt


def test_strategist_and_executors_present():
    agents = build_agents()
    assert "structure-strategist" in agents
    assert {"protenix-agent", "af3-agent", "mhcfine-agent",
            "tcrdock-agent", "affinetune-agent"} <= set(agents)
    assert {"qc-agent", "report-agent"} <= set(agents)


def test_strategist_reads_registry_and_can_delegate():
    a = build_agents()["structure-strategist"]
    assert any("list_structure_tools" in t for t in a.tools)
    assert "Agent" in a.tools


def test_executors_can_drive_browser_and_record():
    ex = build_agents()["protenix-agent"]
    assert any("playwright" in t for t in ex.tools)
    assert any("record_fold_result" in t for t in ex.tools)


def test_options_wire_new_tool_and_agents(tmp_path):
    opts = build_options(str(tmp_path / "run"))
    assert any("list_structure_tools" in t for t in opts.allowed_tools)
    assert "structure-strategist" in opts.agents
    assert opts.permission_mode == "bypassPermissions"


def test_prompt_supports_optional_question(tmp_path):
    base = orchestrator_prompt("x.csv", str(tmp_path), 8)
    for kw in ["ingest", "annotate", "fold", "qc", "report"]:
        assert kw in base.lower()
    steered = orchestrator_prompt("x.csv", str(tmp_path), 8, question="which clones are presented in class II")
    assert "class ii" in steered.lower()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./.venv/bin/python -m pytest tests/test_agents_config.py -q`
Expected: FAIL (no `structure-strategist`; `orchestrator_prompt` takes no `question`)

- [ ] **Step 3: Write minimal implementation**

Rewrite `agents.py`:

```python
from __future__ import annotations
from claude_agent_sdk import ClaudeAgentOptions, AgentDefinition
from .agent_tools import build_server
from . import structure_tools

_EXEC_TOOLS = ["mcp__rep2struct__list_fold_jobs", "mcp__rep2struct__record_fold_result",
               "mcp__playwright__*"]


def _executor(name, tool):
    return AgentDefinition(
        description=f"Folds the {tool} group by driving its Colab notebook through the browser.",
        prompt=(
            f"You run the {tool} structure tool for the jobs assigned to your group. "
            f"Call list_fold_jobs, and for each job whose tool is '{tool}', drive the "
            f"{tool} Colab notebook with the mcp__playwright tools: open it, submit the "
            f"construct (MSA already embedded, run MSA-free at runtime), wait, download "
            f"the model or score, then call record_fold_result with tool='{tool}'. The "
            f"loop is resumable; skip jobs already recorded. If the {tool} Colab adapter "
            f"is not available in this environment, report the job as not-run. Never "
            f"fabricate a model or a score."),
        tools=_EXEC_TOOLS,
        model="sonnet",
    )


def build_agents():
    agents = {
        "structure-strategist": AgentDefinition(
            description="Reasons over the tool registry and construct tags; routes each group to a tool.",
            prompt=(
                "You choose structure tools. Call list_structure_tools to read each tool's "
                "validity domain, and list_fold_jobs to see the jobs (each carries group_id, "
                "mhc_class, has_tcr, species, output_needed). For each homogeneous group pick "
                "ONE tool: Protenix is the default workhorse; switch to a specialized tool "
                "only when the group justifies it (af3 if the user has weights and it helps, "
                "affinetune for is-it-presented in class I or II, mhcfine for a precise class I "
                "pose, tcrdock for the TCR interface). Never Boltz. If no tool's validity "
                "domain covers a group, fall back to Protenix and state plainly that an "
                "un-wired tool would fit better. Justify each choice in one sentence, then "
                "delegate the group to that tool's executor agent (protenix-agent, af3-agent, "
                "mhcfine-agent, tcrdock-agent, affinetune-agent)."),
            tools=["mcp__rep2struct__list_structure_tools", "mcp__rep2struct__list_fold_jobs", "Agent"],
            model="opus",
        ),
        "qc-agent": AgentDefinition(
            description="Skeptical QC per group; calibration is per tool, output-type aware.",
            prompt=(
                "You are a skeptical structural referee. For each folded clonotype call "
                "qc_structure with that group's OWN scramble_threshold (never a global one) "
                "and its output_type and tool. A clean fold does NOT confirm specificity. For "
                "binding_score tools you judge predicted presentation, not geometry. Report "
                "reliable only when the CDR3 to peptide contact beats the group's scramble "
                "calibration; otherwise suspect. Never upgrade a verdict to please the caller."),
            tools=["mcp__rep2struct__qc_structure"],
            model="opus",
        ),
        "report-agent": AgentDefinition(
            description="Renders the final HTML report tying clonotype to specificity to structure to QC.",
            prompt="Call render_final_report and return the report path. Add no unsupported claims.",
            tools=["mcp__rep2struct__render_final_report"],
            model="sonnet",
        ),
    }
    for t in structure_tools.REGISTRY:
        agents[f"{t.name}-agent"] = _executor(f"{t.name}-agent", t.name)
    return agents


def build_options(run_dir):
    return ClaudeAgentOptions(
        mcp_servers={
            "rep2struct": build_server(),
            "playwright": {"command": "npx", "args": ["@playwright/mcp@latest"]},
        },
        agents=build_agents(),
        allowed_tools=[
            "Agent",
            "mcp__rep2struct__ingest_repertoire", "mcp__rep2struct__annotate_specificity",
            "mcp__rep2struct__prep_and_select", "mcp__rep2struct__list_structure_tools",
            "mcp__rep2struct__list_fold_jobs", "mcp__rep2struct__record_fold_result",
            "mcp__rep2struct__qc_structure", "mcp__rep2struct__render_final_report",
            "mcp__playwright__*",
        ],
        permission_mode="bypassPermissions",
    )


def orchestrator_prompt(csv_path, run_dir, top_n, question=None):
    steer = (f"\nUser question steering the routing: {question}\n" if question else
             "\nNo user question: route each group from its construct tags.\n")
    return (
        f"Run the repertoire to structure pipeline on {csv_path} with run_dir {run_dir}.{steer}"
        f"1. Call ingest_repertoire, then annotate_specificity (honest annotation, keep unannotatable as is).\n"
        f"2. Call prep_and_select with top_n {top_n}.\n"
        f"3. Delegate to the structure-strategist to route each group to a tool and drive its executor to fold.\n"
        f"4. Delegate to the qc-agent to QC each folded clonotype with that group's calibration.\n"
        f"5. Delegate to the report-agent to render the final HTML report, and return its path.")
```

- [ ] **Step 4: Run the full suite to verify nothing regressed**

Run: `./.venv/bin/python -m pytest -q`
Expected: PASS (all prior tests plus the new ones; strategist/executor config asserted)

- [ ] **Step 5: Commit**

```bash
git add src/rep2struct/agents.py tests/test_agents_config.py
git commit -m "feat: structure-strategist routes homogeneous groups to per-tool executors

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Notes for the implementer

- The MSA runners (Task 3) and the non-Protenix Colab adapters (Task 8) are the only
  parts that touch external resources. Everything else is fully unit-tested offline.
- v1 does NOT build class II construct geometry (`build_construct` stays class I). The
  routing, registry, grouping, QC branching, and report honesty all support class II
  and binding_score already, so wiring class II constructs later is additive.
- Wiring a real non-Protenix tool later = write its Colab notebook adapter + confirm the
  executor's Playwright sequence against it. The executor prompt already forbids
  fabricating a result when its adapter is absent.
- After Task 8, verify end to end with the verify skill against a tiny fixture CSV
  before any real fold run, and do not disturb the Protenix session running elsewhere.
- KNOWN v1 SIMPLIFICATION vs the spec: the spec asks that an out-of-domain group
  routed to fallback Protenix carry a reservation VISIBLE IN THE REPORT. This plan
  surfaces that reservation through the strategist's one-sentence justification in the
  run narration, not as a structured report cell. Encoding it as a `FoldJob.reservation`
  field plus a report column is additive and deferred; call it out if the run produces
  an out-of-domain group so it is not silently lost.
