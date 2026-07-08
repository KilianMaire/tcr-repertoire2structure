# Wire tcrdock, mhcfine, affinetune (offline layer) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Deliver the offline-testable layer that wires tcrdock, mhcfine, and affinetune into the structure-strategist pipeline (class I): per-tool input builders, a two-layer QC (common validity gate plus tool-specific skeptical control), a `qc_metric` registry field, notebook scaffolds, and report honesty for pose/validity, leaving the live notebook cells and real scramble calibration to be validated together per pipeline afterward.

**Architecture:** QC metrics operate on pure `dict[str, np.ndarray]` chain-coordinate maps (unit-testable with synthetic arrays), with thin CIF wrappers reusing the existing `_heavy_by_chain` parser. Tools carry a `qc_metric` in the registry so the QC control is chosen by the tool, not guessed. Input builders mirror `scripts/build_protenix_inputs.py`, reusing shared construct helpers. Notebook builders mirror the Protenix notebook pattern with scaffolded live cells.

**Tech Stack:** Python 3.11, numpy, Biopython (MMCIFParser), jinja2, dataclasses, pytest. No new dependencies.

## Global Constraints

- Python 3.11; run tests with `./.venv/bin/python -m pytest -q` from the repo root.
- Protenix is the default tool; NEVER Boltz. af3 stays deferred (non-redistributable weights).
- Class I only this round; MHC II construct geometry and Der p 1 are out of scope.
- Honesty is enforced in schemas/output: mhcfine output is a POSE, never a "fold" or a proof of TCR recognition; a binding_score row is never "fold"/"structure"; layer-2 skeptical distances are never cross-compared across tools; a clean fold never confirms specificity.
- Canonical chain IDs: A=TCR alpha, B=TCR beta, C=MHC heavy, D=beta2-microglobulin, E=peptide. mhcfine emits C, D, E only.
- Docs/report human-facing copy: no dash-as-punctuation (use commas, periods, parentheses).
- This repo KEEPS the commit trailer `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.
- Do NOT modify the Protenix fold code path (`scripts/build_protenix_inputs.py` logic, `scripts/build_colab_notebook.py`) except the DRY extraction explicitly called for in Task 5. Do NOT stage anything under `src/rep2struct/data/`.
- Full pytest suite must stay green (currently 62/62). Targeted `git add` per task.

---

### Task 1: registry qc_metric field

**Files:**
- Modify: `src/rep2struct/structure_tools.py` (StructureTool dataclass; each REGISTRY entry; `as_dicts`)
- Test: `tests/test_structure_tools.py` (append)

**Interfaces:**
- Consumes: existing REGISTRY.
- Produces:
  - `StructureTool` gains `qc_metric: str = "cdr3_peptide"`.
  - Entries: protenix `cdr3_peptide`, tcrdock `cdr3_peptide`, mhcfine `peptide_groove`, affinetune `binding_score`, af3 `cdr3_peptide`.
  - `qc_metric_for(name: str) -> str` returns the tool's qc_metric, default `"cdr3_peptide"` for unknown.
  - `as_dicts()` includes `"qc_metric"`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_structure_tools.py (append)
def test_qc_metric_per_tool():
    from rep2struct import structure_tools as st
    by = {t.name: t.qc_metric for t in st.REGISTRY}
    assert by == {"protenix": "cdr3_peptide", "tcrdock": "cdr3_peptide",
                  "mhcfine": "peptide_groove", "affinetune": "binding_score",
                  "af3": "cdr3_peptide"}


def test_qc_metric_for_defaults():
    from rep2struct import structure_tools as st
    assert st.qc_metric_for("mhcfine") == "peptide_groove"
    assert st.qc_metric_for("unknown") == "cdr3_peptide"


def test_as_dicts_exposes_qc_metric():
    from rep2struct import structure_tools as st
    assert all("qc_metric" in d for d in st.as_dicts())
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./.venv/bin/python -m pytest tests/test_structure_tools.py -q -k qc_metric`
Expected: FAIL with `AttributeError: 'StructureTool' object has no attribute 'qc_metric'`

- [ ] **Step 3: Write minimal implementation**

Add the field to the dataclass (after `colab_adapter`, before `is_default`):

```python
    qc_metric: str = "cdr3_peptide"     # cdr3_peptide | peptide_groove | binding_score
```

Add `qc_metric=` to each REGISTRY entry with the values above (protenix/tcrdock/af3 `"cdr3_peptide"`, mhcfine `"peptide_groove"`, affinetune `"binding_score"`).

Add the helper near `output_type_for`:

```python
def qc_metric_for(name: str) -> str:
    for t in REGISTRY:
        if t.name == name:
            return t.qc_metric
    return "cdr3_peptide"
```

In `as_dicts`, add `"qc_metric": t.qc_metric` to each emitted dict.

- [ ] **Step 4: Run test to verify it passes**

Run: `./.venv/bin/python -m pytest tests/test_structure_tools.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/rep2struct/structure_tools.py tests/test_structure_tools.py
git commit -m "feat: qc_metric registry field so the QC control is chosen by the tool

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: common QC layer (validity gate)

**Files:**
- Modify: `src/rep2struct/qc.py` (add `load_chains`, `common_checks`)
- Test: `tests/test_qc.py` (append)

**Interfaces:**
- Consumes: `_heavy_by_chain` (existing).
- Produces:
  - `load_chains(cif_path) -> dict[str, np.ndarray]` = thin alias of `_heavy_by_chain` (public name for reuse by Task 4).
  - `common_checks(chains: dict, expected: set[str]) -> dict` PURE over a chain-coordinate map. Returns `{"ok": bool, "issues": list[str], "n_chains": int, "has_peptide": bool, "min_interatomic": float|None}`. Rules: all `expected` chains present (else issue "missing chains ..."); all coords finite (else "non-finite coords"); peptide chain "E" present and non-empty sets `has_peptide`; `min_interatomic` is the smallest nonzero distance between atoms of DIFFERENT chains, and a value `< 0.5` (Angstrom) adds issue "severe steric clash". `ok` is True iff `issues` is empty.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_qc.py (append)
import numpy as np
from rep2struct.qc import common_checks


def _chain(*xyz):
    return np.array(xyz, dtype=float)


def test_common_checks_passes_a_sane_two_chain_model():
    chains = {"C": _chain([0, 0, 0], [10, 0, 0]), "E": _chain([5, 0, 0], [6, 0, 0])}
    r = common_checks(chains, expected={"C", "E"})
    assert r["ok"] and r["issues"] == [] and r["has_peptide"] and r["n_chains"] == 2


def test_common_checks_flags_missing_chain():
    chains = {"C": _chain([0, 0, 0])}
    r = common_checks(chains, expected={"C", "E"})
    assert not r["ok"] and any("missing" in i for i in r["issues"]) and not r["has_peptide"]


def test_common_checks_flags_nonfinite_coords():
    chains = {"C": _chain([0, 0, 0]), "E": _chain([np.nan, 0, 0])}
    r = common_checks(chains, expected={"C", "E"})
    assert not r["ok"] and any("non-finite" in i for i in r["issues"])


def test_common_checks_flags_severe_clash():
    # two different chains with atoms 0.1A apart
    chains = {"C": _chain([0, 0, 0]), "E": _chain([0.1, 0, 0])}
    r = common_checks(chains, expected={"C", "E"})
    assert not r["ok"] and any("clash" in i for i in r["issues"])
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./.venv/bin/python -m pytest tests/test_qc.py -q -k common_checks`
Expected: FAIL with `ImportError: cannot import name 'common_checks'`

- [ ] **Step 3: Write minimal implementation**

```python
# in src/rep2struct/qc.py
def load_chains(cif_path):
    return _heavy_by_chain(cif_path)


def common_checks(chains: dict, expected: set) -> dict:
    issues = []
    missing = expected - set(chains)
    if missing:
        issues.append(f"missing chains {sorted(missing)}")
    finite = all(np.isfinite(a).all() for a in chains.values())
    if not finite:
        issues.append("non-finite coords")
    has_peptide = "E" in chains and len(chains["E"]) > 0
    min_inter = None
    if finite and len(chains) >= 2:
        ids = list(chains)
        best = np.inf
        for i in range(len(ids)):
            for j in range(i + 1, len(ids)):
                a, b = chains[ids[i]], chains[ids[j]]
                d = np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(-1))
                best = min(best, float(d.min()))
        min_inter = best
        if min_inter < 0.5:
            issues.append("severe steric clash")
    return {"ok": not issues, "issues": issues, "n_chains": len(chains),
            "has_peptide": has_peptide, "min_interatomic": min_inter}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `./.venv/bin/python -m pytest tests/test_qc.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/rep2struct/qc.py tests/test_qc.py
git commit -m "feat: common QC validity gate over a chain-coordinate map

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: peptide_groove metric (mhcfine, layer 2)

**Files:**
- Modify: `src/rep2struct/qc.py` (add `score_pose`, `mean_confidence`, `verdict_groove`)
- Test: `tests/test_qc.py` (append)

**Interfaces:**
- Consumes: nothing new.
- Produces:
  - `score_pose(chains: dict) -> float|None` PURE: count of peptide (chain E) heavy atoms within 4.5 A of MHC heavy chain (chain C) atoms; returns None if C or E absent.
  - `mean_confidence(bfactors: list[float]|None) -> float|None`: mean of the list, or None if empty/None. This is the offline-testable helper for the LIVE confidence path: the executor will read mhcfine's pLDDT-like B-factors and pass the mean into `verdict_groove`. It is intentionally not wired into `qc_structure` this round (no confidence source until live validation), same status as the notebook live cells. Kept and tested now so the live wiring has a verified helper.
  - `verdict_groove(pose_atoms: float|None, threshold: float, clonotype_id: str, tool: str, confidence: float|None = None) -> QCResult`: returns `qc_failed` if pose_atoms is None; else `pose_reliable` when `pose_atoms > threshold`, `pose_suspect` otherwise. The geometric QCResult fields (`cdr3_pep_atoms`, `dockq`, `crossing_angle`) are left None to avoid overloading them with a non-CDR3 quantity; the pose contact count and the confidence are carried in the reason string. Reason says "peptide in groove" and "pose", never "fold"/"structure"/"recognition". Sets `tool` and `calibration_basis="groove_scramble_null"`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_qc.py (append)
from rep2struct.qc import score_pose, mean_confidence, verdict_groove


def test_score_pose_counts_peptide_mhc_contacts():
    chains = {"C": _chain([0, 0, 0], [100, 0, 0]), "E": _chain([1, 0, 0], [50, 0, 0])}
    # only the E atom at (1,0,0) is within 4.5A of a C atom
    assert score_pose(chains) == 1.0
    assert score_pose({"E": _chain([0, 0, 0])}) is None  # no MHC heavy chain


def test_mean_confidence():
    assert mean_confidence([80.0, 90.0]) == 85.0
    assert mean_confidence(None) is None
    assert mean_confidence([]) is None


def test_verdict_groove_is_a_pose_not_a_fold():
    hi = verdict_groove(20.0, 10.0, "c1", tool="mhcfine", confidence=88.0)
    lo = verdict_groove(5.0, 10.0, "c2", tool="mhcfine")
    assert hi.qc_verdict == "pose_reliable" and lo.qc_verdict == "pose_suspect"
    for r in (hi, lo):
        assert "pose" in r.reason.lower()
        for bad in ("fold", "structure", "recognition"):
            assert bad not in r.reason.lower()
    assert hi.tool == "mhcfine" and hi.calibration_basis == "groove_scramble_null"
    assert verdict_groove(None, 10.0, "c3", tool="mhcfine").qc_verdict == "qc_failed"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./.venv/bin/python -m pytest tests/test_qc.py -q -k "pose or groove or confidence"`
Expected: FAIL with `ImportError: cannot import name 'score_pose'`

- [ ] **Step 3: Write minimal implementation**

```python
# in src/rep2struct/qc.py
def score_pose(chains: dict):
    if "C" not in chains or "E" not in chains:
        return None
    mhc, pep = chains["C"], chains["E"]
    d = np.sqrt(((pep[:, None, :] - mhc[None, :, :]) ** 2).sum(-1))
    return float((d < 4.5).sum())


def mean_confidence(bfactors):
    if not bfactors:
        return None
    return float(sum(bfactors) / len(bfactors))


def verdict_groove(pose_atoms, threshold: float, clonotype_id: str, tool: str,
                   confidence=None) -> QCResult:
    if pose_atoms is None:
        return QCResult(clonotype_id, "qc_failed", "no MHC-peptide pose to score", tool=tool)
    ok = pose_atoms > threshold
    reason = (f"peptide in groove contact {pose_atoms:.0f} "
              f"{'beats' if ok else 'not above'} scramble null; "
              f"model confidence {confidence}")
    return QCResult(clonotype_id, "pose_reliable" if ok else "pose_suspect", reason,
                    tool=tool, calibration_basis="groove_scramble_null")
```

- [ ] **Step 4: Run test to verify it passes**

Run: `./.venv/bin/python -m pytest tests/test_qc.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/rep2struct/qc.py tests/test_qc.py
git commit -m "feat: peptide_groove pose metric for mhcfine (honest pose, not recognition)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: qc_structure two-layer wiring

**Files:**
- Modify: `src/rep2struct/agent_tools.py` (`qc_structure`, 116-147)
- Test: `tests/test_agent_tools.py` (append)

**Interfaces:**
- Consumes: `structure_tools.qc_metric_for`, `qc.load_chains`, `qc.common_checks`, `qc.score_model`, `qc.verdict`, `qc.score_pose`, `qc.verdict_groove`, `qc.verdict_binding`.
- Produces: `qc_structure` now: for `binding_score` tools, unchanged (float + verdict_binding). For structure tools (`cdr3_peptide`, `peptide_groove`), it (1) `load_chains(paths[0])`, (2) runs `common_checks(chains, expected)` where expected is `{"A","B","C","D","E"}` for `cdr3_peptide` and `{"C","D","E"}` for `peptide_groove`; if not `ok`, returns `qc_failed` with the joined issues; (3) else dispatches: `cdr3_peptide` -> existing `score_model`/`verdict`; `peptide_groove` -> `score_pose`/`verdict_groove`. `res.tool` and `res.calibration_basis` set as today (scramble_null for cdr3_peptide; verdict_groove sets its own).

**Interface note:** keep the binding_score branch first (it reads a float, not a CIF). Derive the metric via `structure_tools.qc_metric_for(tool)`. The existing `test_qc_tool_flags_scramble` (protenix, cdr3_peptide, 5-chain fixture) must still pass.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_agent_tools.py (append)
def test_qc_structure_common_gate_fails_closed_on_missing_chains(tmp_path):
    # a mhcfine pose recorded but the CIF is a fixture missing chain C/E -> qc_failed
    rd = str(tmp_path / "run")
    bad = str(_ROOT / "tests" / "fixtures" / "threechain_min.cif")  # implementer: pick a fixture that lacks the expected mhcfine chains
    _run(at.record_fold_result.handler(
        {"run_dir": rd, "clonotype_id": "c1", "model_paths": [bad], "tool": "mhcfine"}))
    res = _run(at.qc_structure.handler(
        {"run_dir": rd, "clonotype_id": "c1", "scramble_threshold": 1.0,
         "output_type": "structure", "tool": "mhcfine"}))
    assert res["structuredContent"]["qc_verdict"] == "qc_failed"
```

**Note for the implementer:** define `_ROOT = Path(__file__).resolve().parents[1]` if not already present, and READ the three fixtures in tests/fixtures to choose one that genuinely lacks the mhcfine-expected chains (C,D,E) OR construct a tiny CIF in tmp_path with only chain C so the gate fails. The point of the test is that the common gate fails closed before the pose metric runs. Also add a positive test: a synthetic CIF (or a fixture) with valid C,D,E where mhcfine yields pose_reliable/pose_suspect. If building a CIF inline is impractical, unit-test the pose/groove path directly in test_qc.py (already done in Task 3) and keep this test focused on the fail-closed gate.

- [ ] **Step 2: Run test to verify it fails**

Run: `./.venv/bin/python -m pytest tests/test_agent_tools.py -q -k common_gate`
Expected: FAIL (mhcfine currently routed through the cdr3 structure path, no gate)

- [ ] **Step 3: Write minimal implementation**

Replace the structure branch of `qc_structure` so it gates then dispatches on the metric:

```python
    from .qc import (verdict_binding, load_chains, common_checks,
                     score_model, verdict, score_pose, verdict_groove)
    ...
    tool = args.get("tool", rec.get("tool", "protenix"))
    metric = structure_tools.qc_metric_for(tool)
    if not paths:
        res = QCResult(args["clonotype_id"], "qc_failed", "no model recorded", tool=tool)
    elif metric == "binding_score":
        score = float(Path(paths[0]).read_text().strip())
        res = verdict_binding(score, args["scramble_threshold"], args["clonotype_id"], tool=tool)
    else:
        expected = {"A", "B", "C", "D", "E"} if metric == "cdr3_peptide" else {"C", "D", "E"}
        chains = load_chains(paths[0])
        cc = common_checks(chains, expected)
        if not cc["ok"]:
            res = QCResult(args["clonotype_id"], "qc_failed",
                           "; ".join(cc["issues"]), tool=tool)
        elif metric == "peptide_groove":
            res = verdict_groove(score_pose(chains), args["scramble_threshold"],
                                 args["clonotype_id"], tool=tool)
        else:  # cdr3_peptide
            s = score_model(paths[0])
            s["clonotype_id"] = args["clonotype_id"]
            res = verdict(s, args["scramble_threshold"])
            res.tool = tool
            res.calibration_basis = "scramble_null"
```

(Keep the qc-stage persistence and return exactly as they are now.)

- [ ] **Step 4: Run test to verify it passes**

Run: `./.venv/bin/python -m pytest tests/test_agent_tools.py -q`
Expected: PASS (new fail-closed test passes; existing test_qc_tool_flags_scramble still passes via cdr3_peptide)

- [ ] **Step 5: Commit**

```bash
git add src/rep2struct/agent_tools.py tests/test_agent_tools.py
git commit -m "feat: qc_structure runs the common gate then dispatches on qc_metric

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 5: shared construct IO helpers (DRY)

**Files:**
- Create: `src/rep2struct/tools/__init__.py` (empty), `src/rep2struct/tools/construct_io.py`
- Modify: `scripts/build_protenix_inputs.py` (import the shared helpers instead of its local copies)
- Test: `tests/test_construct_io.py`

**Interfaces:**
- Produces:
  - `parse_fasta(text) -> dict[str,str]` (moved from build_protenix_inputs.py).
  - `scramble_peptide(pep) -> str` (moved from build_protenix_inputs.py, same deterministic reverse-then-rotate logic).
  - `pmhc_only(chains: dict) -> dict` returns a copy with chains A and B removed (keeps C, D, E) for the no-TCR tools.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_construct_io.py
from rep2struct.tools.construct_io import parse_fasta, scramble_peptide, pmhc_only


def test_parse_fasta_roundtrip():
    c = parse_fasta(">A\nAAA\n>E\nSII\n")
    assert c == {"A": "AAA", "E": "SII"}


def test_scramble_is_deterministic_and_non_identity():
    assert scramble_peptide("SIINFEKL") == scramble_peptide("SIINFEKL")
    assert scramble_peptide("SIINFEKL") != "SIINFEKL"
    assert sorted(scramble_peptide("SIINFEKL")) == sorted("SIINFEKL")  # same composition


def test_pmhc_only_drops_tcr_chains():
    chains = {"A": "a", "B": "b", "C": "c", "D": "d", "E": "e"}
    assert pmhc_only(chains) == {"C": "c", "D": "d", "E": "e"}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./.venv/bin/python -m pytest tests/test_construct_io.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'rep2struct.tools'`

- [ ] **Step 3: Write minimal implementation**

Create `src/rep2struct/tools/__init__.py` (empty). Create `src/rep2struct/tools/construct_io.py`:

```python
from __future__ import annotations


def parse_fasta(text):
    chains, cur = {}, None
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            cur = line[1:].strip()
            chains[cur] = ""
        elif cur:
            chains[cur] += line
    return chains


def scramble_peptide(pep):
    # deterministic non-identity shuffle: reverse then rotate by 1 (reproducible).
    s = pep[::-1]
    s = s[1:] + s[:1]
    return s if s != pep else pep[1:] + pep[:1]


def pmhc_only(chains: dict) -> dict:
    return {k: v for k, v in chains.items() if k not in ("A", "B")}
```

In `scripts/build_protenix_inputs.py`, delete the local `parse_fasta` and `scramble_peptide` and import them: add `sys.path.insert` is already implicit for scripts; import via `from rep2struct.tools.construct_io import parse_fasta, scramble_peptide`. Verify the script still runs (it is not covered by the suite, so just import-check).

- [ ] **Step 4: Run test to verify it passes**

Run: `./.venv/bin/python -m pytest tests/test_construct_io.py -q && ./.venv/bin/python -c "import ast; ast.parse(open('scripts/build_protenix_inputs.py').read())"`
Expected: PASS and no syntax error.

- [ ] **Step 5: Commit**

```bash
git add src/rep2struct/tools/__init__.py src/rep2struct/tools/construct_io.py scripts/build_protenix_inputs.py tests/test_construct_io.py
git commit -m "refactor: shared construct IO helpers (parse_fasta, scramble, pmhc_only)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 6: per-tool input builders

**Files:**
- Create: `src/rep2struct/tools/tcrdock_inputs.py`, `src/rep2struct/tools/mhcfine_inputs.py`, `src/rep2struct/tools/affinetune_inputs.py`
- Test: `tests/test_tool_inputs.py`

**Interfaces:**
- Consumes: `construct_io.parse_fasta`, `scramble_peptide`, `pmhc_only`.
- Produces, each module exposes `build(construct_fasta: str) -> dict` returning `{"cognate": <inputs>, "scramble": <inputs>}`:
  - tcrdock: full A-E chains; input is `{"chains": {id: seq}, ...}` for all five; scramble shuffles E.
  - mhcfine: `pmhc_only` (C,D,E); input carries only those; scramble shuffles E.
  - affinetune: class I pMHC; input is `{"mhc": chains["C"], "b2m": chains["D"], "peptide": chains["E"]}`; scramble shuffles the peptide.
  Each also exposes `TOOL = "<name>"` for traceability. The exact serialization to each tool's real format is refined during live validation; the builders lock the chain selection, the cognate/scramble pair, and the deterministic scramble.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_tool_inputs.py
from rep2struct.tools import tcrdock_inputs, mhcfine_inputs, affinetune_inputs

FASTA = ">A\nAAAA\n>B\nBBBB\n>C\nCCCC\n>D\nDDDD\n>E\nSIINFEKL\n"


def test_tcrdock_keeps_all_five_chains_and_scrambles_peptide():
    out = tcrdock_inputs.build(FASTA)
    assert set(out["cognate"]["chains"]) == {"A", "B", "C", "D", "E"}
    assert out["scramble"]["chains"]["E"] != "SIINFEKL"
    assert sorted(out["scramble"]["chains"]["E"]) == sorted("SIINFEKL")


def test_mhcfine_drops_tcr_chains():
    out = mhcfine_inputs.build(FASTA)
    assert set(out["cognate"]["chains"]) == {"C", "D", "E"}
    assert "A" not in out["cognate"]["chains"]


def test_affinetune_maps_class_i_fields_and_scrambles():
    out = affinetune_inputs.build(FASTA)
    assert out["cognate"]["mhc"] == "CCCC" and out["cognate"]["b2m"] == "DDDD"
    assert out["cognate"]["peptide"] == "SIINFEKL"
    assert out["scramble"]["peptide"] != "SIINFEKL"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./.venv/bin/python -m pytest tests/test_tool_inputs.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'rep2struct.tools.tcrdock_inputs'`

- [ ] **Step 3: Write minimal implementation**

`src/rep2struct/tools/tcrdock_inputs.py`:

```python
from __future__ import annotations
from .construct_io import parse_fasta, scramble_peptide

TOOL = "tcrdock"


def build(construct_fasta: str) -> dict:
    chains = parse_fasta(construct_fasta)
    sc = dict(chains)
    sc["E"] = scramble_peptide(chains["E"])
    return {"cognate": {"chains": chains}, "scramble": {"chains": sc}}
```

`src/rep2struct/tools/mhcfine_inputs.py`:

```python
from __future__ import annotations
from .construct_io import parse_fasta, scramble_peptide, pmhc_only

TOOL = "mhcfine"


def build(construct_fasta: str) -> dict:
    chains = pmhc_only(parse_fasta(construct_fasta))
    sc = dict(chains)
    sc["E"] = scramble_peptide(chains["E"])
    return {"cognate": {"chains": chains}, "scramble": {"chains": sc}}
```

`src/rep2struct/tools/affinetune_inputs.py`:

```python
from __future__ import annotations
from .construct_io import parse_fasta, scramble_peptide

TOOL = "affinetune"


def build(construct_fasta: str) -> dict:
    chains = parse_fasta(construct_fasta)
    def rec(pep):
        return {"mhc": chains["C"], "b2m": chains["D"], "peptide": pep}
    return {"cognate": rec(chains["E"]),
            "scramble": rec(scramble_peptide(chains["E"]))}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `./.venv/bin/python -m pytest tests/test_tool_inputs.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/rep2struct/tools/tcrdock_inputs.py src/rep2struct/tools/mhcfine_inputs.py src/rep2struct/tools/affinetune_inputs.py tests/test_tool_inputs.py
git commit -m "feat: per-tool class I input builders (cognate + scramble, chain selection)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 7: notebook scaffolds

**Files:**
- Create: `src/rep2struct/tools/notebook.py`
- Test: `tests/test_tool_notebook.py`

**Interfaces:**
- Consumes: nothing (pure notebook JSON assembly).
- Produces:
  - `build_notebook(tool: str, inputs: dict) -> dict` returns a minimal, valid Jupyter notebook (nbformat 4) as a Python dict with three cells: (1) a code cell embedding `inputs` as JSON, (2) a code cell containing a clearly marked live TODO (`# TODO(live): invoke {tool} here ...`) and a `raise NotImplementedError` so an un-validated notebook cannot silently produce a fake result, (3) a code cell writing the output path. The notebook must be JSON-serializable and carry `nbformat == 4`.

**Interface note:** this is the scaffold only. The real invocation cell is filled during live validation. The `raise NotImplementedError` in the TODO cell is the code-level guarantee that an unwired notebook fails loudly instead of fabricating output, mirroring the executor's "report not-run, never fabricate" rule.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_tool_notebook.py
import json
from rep2struct.tools.notebook import build_notebook


def test_notebook_is_valid_and_embeds_inputs():
    nb = build_notebook("tcrdock", {"cognate": {"chains": {"E": "SII"}}})
    assert nb["nbformat"] == 4
    json.dumps(nb)  # serializable
    src = "".join(src_ for cell in nb["cells"] for src_ in cell["source"])
    assert "SII" in src                       # inputs embedded
    assert "TODO(live)" in src                 # live marker present


def test_notebook_scaffold_fails_loud_not_fake():
    nb = build_notebook("mhcfine", {})
    src = "".join(s for cell in nb["cells"] for s in cell["source"])
    assert "NotImplementedError" in src        # unwired notebook cannot fake a result
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./.venv/bin/python -m pytest tests/test_tool_notebook.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'rep2struct.tools.notebook'`

- [ ] **Step 3: Write minimal implementation**

```python
# src/rep2struct/tools/notebook.py
from __future__ import annotations
import json


def _code(*lines):
    return {"cell_type": "code", "metadata": {}, "execution_count": None,
            "outputs": [], "source": list(lines)}


def build_notebook(tool: str, inputs: dict) -> dict:
    return {
        "nbformat": 4, "nbformat_minor": 5, "metadata": {},
        "cells": [
            _code(f"# {tool} inputs (embedded, MSA-free at runtime)\n",
                  "INPUTS = ", json.dumps(inputs)),
            _code(f"# TODO(live): invoke {tool} here against INPUTS and write the model/score.\n",
                  "raise NotImplementedError('live cell not yet validated')\n"),
            _code("# TODO(live): save the output (CIF for structure, float for binding) to OUT_PATH\n"),
        ],
    }
```

- [ ] **Step 4: Run test to verify it passes**

Run: `./.venv/bin/python -m pytest tests/test_tool_notebook.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/rep2struct/tools/notebook.py tests/test_tool_notebook.py
git commit -m "feat: notebook scaffold generator (embeds inputs, fails loud until live cell added)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 8: report validity column and pose honesty labels

**Files:**
- Modify: `src/rep2struct/report.py` (EVIDENCE map, row build, render), `src/rep2struct/templates/report.html.j2`
- Test: `tests/test_report.py` (append)

**Interfaces:**
- Consumes: QCResult with new verdicts `pose_reliable`/`pose_suspect`; optional per-clonotype common-validity descriptors.
- Produces:
  - `EVIDENCE` gains `"pose_reliable": "pose (peptide in groove)"`, `"pose_suspect": "pose (peptide in groove)"`. A pose row must never read "fold"/"structure".
  - `render_report(clonotypes, annotations, qc_results, metrics=None, msa_basis=None, validity=None)` gains an optional `validity` map `{clonotype_id: str}` (the layer-1 summary, e.g. "valid" or the joined issues). Existing call sites keep working (validity optional, defaults empty). Each row exposes `validity` for a new "validity (common)" column.

**Interface note:** READ report.py and report.html.j2 first; add the validity cell the same minimal way the tool/evidence cells were added. Do not sort rows. No dash-as-punctuation in added copy.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_report.py (append)
def test_pose_row_labelled_pose_not_fold_or_structure():
    c, a = _fixtures()
    q = QCResult("c1", "pose_reliable", "peptide in groove contact 20 beats scramble null; model confidence 88.0",
                 tool="mhcfine", calibration_basis="groove_scramble_null")
    html = render_report(c, a, [q])
    row = [ln for ln in html.splitlines() if "mhcfine" in ln]
    assert row
    assert any("pose" in ln.lower() for ln in row)
    # the evidence label for this row must not call it a fold or a structure
    ev = [ln for ln in html.splitlines() if "peptide in groove" in ln.lower() or "pose (" in ln.lower()]
    assert ev and all("fold" not in ln.lower() and "structure" not in ln.lower() for ln in ev)


def test_validity_column_renders_when_supplied():
    c, a = _fixtures()
    q = QCResult("c1", "reliable", "ok", tool="protenix")
    html = render_report(c, a, [q], validity={"c1": "valid"})
    assert "valid" in html.lower()


def test_render_report_validity_is_optional():
    c, a = _fixtures()
    q = QCResult("c1", "reliable", "ok", tool="protenix")
    render_report(c, a, [q])  # must not raise without validity
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./.venv/bin/python -m pytest tests/test_report.py -q -k "pose or validity"`
Expected: FAIL (no pose label / no validity param)

- [ ] **Step 3: Write minimal implementation**

Add to `EVIDENCE`:

```python
    "pose_reliable": "pose (peptide in groove)",
    "pose_suspect": "pose (peptide in groove)",
```

Change `render_report` signature to add `validity=None`, default `validity = validity or {}`, and add `"validity": validity.get(c.id, "n/a")` to each row dict. In `report.html.j2`, add a "Validity (common)" header and a `<td>{{ row.validity }}</td>` cell, following the existing Tool/Evidence cell pattern. Do not reorder rows.

- [ ] **Step 4: Run test to verify it passes**

Run: `./.venv/bin/python -m pytest tests/test_report.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/rep2struct/report.py src/rep2struct/templates/report.html.j2 tests/test_report.py
git commit -m "feat: report labels mhcfine as pose and shows the common validity column

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Notes for the implementer

- The QC metric functions are pure over `dict[str, np.ndarray]` chain maps, so they are
  unit-tested with synthetic coordinates and do not depend on fixture internals. Only
  Task 4's fail-closed test touches a real fixture; inspect tests/fixtures first.
- The three tool notebooks' real invocation cells and the executor Playwright drives are
  NOT in this plan. They are validated live, one pipeline at a time, after this merges.
  Until then each tool stays "not-run" (executor rule) and the scaffold raises
  NotImplementedError, so nothing fabricates a result.
- Do not wire the input builders or notebooks into prep_and_select/executors yet: that
  binding happens during live validation per pipeline. This plan delivers the tested
  building blocks.
- After all tasks, run the full suite and confirm the Protenix path is unchanged
  (test_qc_tool_flags_scramble, agents/agent_tools tests green).
