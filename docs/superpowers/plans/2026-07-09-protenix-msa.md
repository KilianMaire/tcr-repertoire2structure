# Protenix pre-fold MSA (v1) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give the Protenix fold a real unpaired MSA per protein chain, computed in a Colab CPU cell before the GPU fold, so the TCR actually docks and the QC verdicts become meaningful.

**Architecture:** `_protenix_notebook` gains one cell (inserted before the write-inputs cell) that deduplicates the unique protein chain sequences, queries the ColabFold MMseqs2 API for each, writes an a3m, injects `unpairedMsaPath` into each `proteinChain` of the embedded INPUTS in place, and writes a per-clonotype MSA manifest under `out/`. The manifest travels back in the existing repatriation zip, and `render_report` uses it to state the actual per-clonotype MSA basis.

**Tech Stack:** Python 3.11, the `colabfold` package (`run_mmseqs2`) on Colab, Protenix input JSON `unpairedMsaPath`, jinja2 report, pytest.

## Global Constraints

- Strict TDD: write the failing test first, watch it fail, then implement. `./.venv/bin/python -m pytest -q` must stay green after every task (currently 132 passing).
- Protenix only. Do not touch tcrdock, affinetune, mhcfine, af3.
- Unpaired MSA only. No paired MSA. Peptide chains (len < 20) get no MSA.
- Notebook cells are self-contained strings that run on Colab with no access to `rep2struct`; the MSA cell logic is inline. Notebook-content tests assert markers and Python validity, not live execution.
- Commit messages keep the trailer `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`. This repo keeps the Claude co-author trailer (hackathon repo).
- Docs use no dashes as punctuation (periods, commas, parentheses instead).

---

## File Structure

- `src/rep2struct/tools/notebook.py` — `_protenix_notebook` gains the MSA cell (Task 1).
- `src/rep2struct/report.py` — `msa_basis_from_manifest` helper + `_msa_note` extension (Task 2).
- `src/rep2struct/agent_tools.py` — `render_final_report` reads the repatriated manifests (Task 2).
- `tests/test_tool_notebook.py` — MSA cell markers (Task 1).
- `tests/test_report.py` — manifest -> basis -> note (Task 2).
- `tests/test_agent_tools.py` — `render_final_report` picks up the manifest (Task 2).
- Live validation drives the real notebook on Colab (Task 3), no code file.

---

## Task 1: MSA cell in the Protenix notebook

**Files:**
- Modify: `src/rep2struct/tools/notebook.py` (`_protenix_notebook`, insert a cell before the write-inputs cell)
- Test: `tests/test_tool_notebook.py`

**Interfaces:**
- Consumes: `build_notebook("protenix", inputs)` where `inputs = {key: <protenix JSON list>}`, key = `{cid}_cognate` / `{cid}_scramble`, each `proteinChain` has `sequence`, `count`, `id`.
- Produces: the notebook, when run on Colab, mutates each `proteinChain` to add `unpairedMsaPath` (for `len(sequence) >= 20`) and writes `out/{cid}_msa_manifest.json` = `{chain_id: {"got_msa": bool}, ...}` (protein chains only, peptide excluded). Task 2 consumes that manifest.

- [ ] **Step 1: Write the failing test**

Add to `tests/test_tool_notebook.py`:

```python
def test_protenix_notebook_computes_pre_fold_msa():
    # v1 MSA: a Colab-CPU cell computes one unpaired a3m per unique protein chain via ColabFold
    # MMseqs2, injects it as unpairedMsaPath, BEFORE the fold. Peptides (len<20) get none, and a
    # failed search folds that chain MSA-free (never fatal).
    inputs = {
        "c1_cognate": [{"name": "cognate", "sequences": [
            {"proteinChain": {"sequence": "M" * 100, "count": 1, "id": ["C"]}},
            {"proteinChain": {"sequence": "GILGFVFTL", "count": 1, "id": ["E"]}}],
            "covalent_bonds": []}],
        "c1_scramble": [{"name": "scramble", "sequences": [
            {"proteinChain": {"sequence": "M" * 100, "count": 1, "id": ["C"]}},
            {"proteinChain": {"sequence": "TFVFGLIGL", "count": 1, "id": ["E"]}}],
            "covalent_bonds": []}],
    }
    src = "".join(s for cell in build_notebook("protenix", inputs)["cells"] for s in cell["source"])
    for marker in ("pip install -q colabfold", "run_mmseqs2", "use_pairing=False",
                   "/content/msa", "unpairedMsaPath", "len(s) >= 20", "_msa_manifest.json",
                   "MSA_FAIL"):
        assert marker in src, f"missing MSA marker: {marker}"
    # the MSA cell must run BEFORE the write-inputs cell (the fold consumes the injected paths)
    assert src.index("run_mmseqs2") < src.index("write each embedded record")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./.venv/bin/python -m pytest tests/test_tool_notebook.py::test_protenix_notebook_computes_pre_fold_msa -q`
Expected: FAIL (markers like `run_mmseqs2` absent from the current notebook).

- [ ] **Step 3: Insert the MSA cell**

In `src/rep2struct/tools/notebook.py`, inside `_protenix_notebook`, insert this cell between the install cell (the one containing `# 1. install Protenix`) and the write-inputs cell (the one containing `# 2. write each embedded record`):

```python
            _code("# 2. MSA (Colab CPU): one unpaired a3m per UNIQUE protein chain via ColabFold\n",
                  "# MMseqs2, injected as unpairedMsaPath. Runs BEFORE the GPU fold so a slow MSA\n",
                  "# server cannot wedge it. Peptides (len<20) get none; a failed search folds that\n",
                  "# chain MSA-free (never fatal).\n",
                  "import os, json, hashlib, subprocess\n",
                  "subprocess.run('pip install -q colabfold', shell=True, check=True)\n",
                  "from colabfold.batch import run_mmseqs2\n",
                  "os.makedirs('/content/msa', exist_ok=True)\n",
                  "uniq = {}\n",
                  "for rec in INPUTS.values():\n",
                  "    for ch in rec[0]['sequences']:\n",
                  "        s = ch['proteinChain']['sequence']\n",
                  "        if len(s) >= 20:\n",
                  "            uniq.setdefault(s, hashlib.sha1(s.encode()).hexdigest()[:12])\n",
                  "seq2path = {}\n",
                  "for s, h in uniq.items():\n",
                  "    try:\n",
                  "        a3m = run_mmseqs2(s, prefix=f'/content/msa/{h}', use_pairing=False)\n",
                  "        a3m = a3m[0] if isinstance(a3m, list) else a3m\n",
                  "        p = f'/content/msa/{h}.a3m'; open(p, 'w').write(a3m)\n",
                  "        seq2path[s] = p; print('MSA_OK', h, flush=True)\n",
                  "    except Exception as e:\n",
                  "        print('MSA_FAIL', h, type(e).__name__, e, flush=True)\n",
                  "os.makedirs('out', exist_ok=True)\n",
                  "manifests = {}\n",
                  "for key, rec in INPUTS.items():\n",
                  "    cid = key.rsplit('_', 1)[0]   # {cid}_cognate / {cid}_scramble -> cid\n",
                  "    man = manifests.setdefault(cid, {})\n",
                  "    for ch in rec[0]['sequences']:\n",
                  "        pc = ch['proteinChain']; s = pc['sequence']\n",
                  "        if s in seq2path:\n",
                  "            pc['unpairedMsaPath'] = seq2path[s]\n",
                  "            man[pc['id'][0]] = {'got_msa': True}\n",
                  "        elif len(s) >= 20:\n",
                  "            man[pc['id'][0]] = {'got_msa': False}\n",
                  "for cid, man in manifests.items():\n",
                  "    json.dump(man, open(f'out/{cid}_msa_manifest.json', 'w'))\n",
                  "print('MSA done; a3m:', len(seq2path), 'manifests:', list(manifests), flush=True)\n"),
```

The existing write-inputs cell needs no change: it dumps `INPUTS`, which this cell has mutated in place to carry `unpairedMsaPath`.

- [ ] **Step 4: Run test to verify it passes**

Run: `./.venv/bin/python -m pytest tests/test_tool_notebook.py::test_protenix_notebook_computes_pre_fold_msa -q`
Expected: PASS.

- [ ] **Step 5: Run the full suite**

Run: `./.venv/bin/python -m pytest -q`
Expected: PASS (133 tests). The existing `test_protenix_notebook_is_wired_not_a_stub` and `test_inputs_cell_is_executable_python_with_bool_and_none` still pass (the INPUTS cell is unchanged and still a valid Python literal).

- [ ] **Step 6: Commit**

```bash
git add src/rep2struct/tools/notebook.py tests/test_tool_notebook.py
git commit -m "feat(msa): pre-fold ColabFold MSA cell in the Protenix notebook

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: Report the actual MSA basis from the manifest

**Files:**
- Modify: `src/rep2struct/report.py` (add `msa_basis_from_manifest`, extend `_msa_note`)
- Modify: `src/rep2struct/agent_tools.py` (`render_final_report` reads `*_msa_manifest.json`)
- Test: `tests/test_report.py`, `tests/test_agent_tools.py`

**Interfaces:**
- Consumes: the per-clonotype manifest from Task 1, `{chain_id: {"got_msa": bool}, ...}`, repatriated into the run dir as `{cid}_msa_manifest.json` (anywhere under it; the executor unzips the fold zip there).
- Produces: `msa_basis_from_manifest(manifest: dict) -> str` returning `"colab_cpu:k/n"` or `"none"`; `_msa_note(basis)` renders `"colab_cpu:3/4"` as `"MSA colab_cpu (3/4 chains)"`. `render_final_report` feeds `render_report(msa_basis=...)` with the manifest-derived basis per clonotype, falling back to the foldjob's prep-time `msa_basis` when no manifest is present.

- [ ] **Step 1: Write the failing tests**

Add to `tests/test_report.py`:

```python
def test_msa_basis_from_manifest():
    from rep2struct.report import msa_basis_from_manifest
    assert msa_basis_from_manifest(
        {"A": {"got_msa": True}, "B": {"got_msa": True},
         "C": {"got_msa": True}, "D": {"got_msa": True}}) == "colab_cpu:4/4"
    assert msa_basis_from_manifest(
        {"A": {"got_msa": True}, "B": {"got_msa": False},
         "C": {"got_msa": True}, "D": {"got_msa": False}}) == "colab_cpu:2/4"
    # nothing got an MSA -> honestly MSA-free
    assert msa_basis_from_manifest({"A": {"got_msa": False}, "B": {"got_msa": False}}) == "none"


def test_msa_note_renders_manifest_basis():
    from rep2struct.report import _msa_note
    assert _msa_note("colab_cpu:3/4") == "MSA colab_cpu (3/4 chains)"
    assert _msa_note("colab_cpu") == "MSA colab_cpu"
    assert _msa_note("none") == "MSA-free (reduced confidence)"
    assert _msa_note(None) == "MSA-free (reduced confidence)"


def test_render_report_shows_manifest_msa_basis():
    from rep2struct.report import render_report
    from rep2struct.schema import Clonotype, Annotation, QCResult
    c = Clonotype(id="c1", trav="", cdr3a="", trbv="", cdr3b="", size=1)
    a = Annotation(clonotype_id="c1", annotatable=True, confidence_tier="high", epitope="X")
    q = QCResult("c1", "reliable", "ok", tool="protenix")
    html = render_report([c], [a], [q], msa_basis={"c1": "colab_cpu:3/4"})
    assert "MSA colab_cpu (3/4 chains)" in html
```

- [ ] **Step 2: Run to verify they fail**

Run: `./.venv/bin/python -m pytest tests/test_report.py -q -k "msa_basis or msa_note or manifest_msa"`
Expected: FAIL (`msa_basis_from_manifest` does not exist; `_msa_note` does not parse `colab_cpu:3/4`).

- [ ] **Step 3: Implement the report helpers**

In `src/rep2struct/report.py`, replace `_msa_note` and add `msa_basis_from_manifest`:

```python
def msa_basis_from_manifest(manifest: dict) -> str:
    """Per-clonotype MSA manifest {chain_id: {"got_msa": bool}, ...} -> basis token.
    "colab_cpu:k/n" when any chain got an MSA, else "none" (honestly MSA-free)."""
    chains = [v for v in manifest.values() if isinstance(v, dict) and "got_msa" in v]
    n = len(chains)
    k = sum(1 for v in chains if v["got_msa"])
    return f"colab_cpu:{k}/{n}" if k else "none"


def _msa_note(basis) -> str:
    if basis and basis.startswith("colab_cpu"):
        _, _, cnt = basis.partition(":")
        return f"MSA colab_cpu ({cnt} chains)" if cnt else "MSA colab_cpu"
    if basis == "local":
        return "MSA local"
    return "MSA-free (reduced confidence)"
```

- [ ] **Step 4: Run to verify the report tests pass**

Run: `./.venv/bin/python -m pytest tests/test_report.py -q -k "msa_basis or msa_note or manifest_msa"`
Expected: PASS.

- [ ] **Step 5: Write the failing render_final_report test**

Add to `tests/test_agent_tools.py`:

```python
def test_render_final_report_reads_msa_manifest(tmp_path):
    # The report states the ACTUAL per-clonotype MSA basis from the repatriated manifest,
    # overriding the prep-time msa_basis on the foldjob.
    import json
    from rep2struct.runstate import RunState
    from rep2struct.schema import Clonotype, Annotation, QCResult, FoldJob
    at.configure()
    rd = str(tmp_path / "run")
    rs = RunState(rd)
    rs.write_stage("ingest", [Clonotype(id="c1", trav="", cdr3a="", trbv="", cdr3b="", size=1)])
    rs.write_stage("annotate", [Annotation(clonotype_id="c1", annotatable=True,
                                           confidence_tier="high", epitope="X")])
    rs.write_stage("foldjobs", [FoldJob(clonotype_id="c1", construct_fasta=">E\nX", msa_basis="none")])
    rs.write_stage("qc", [__import__("dataclasses").asdict(QCResult("c1", "reliable", "ok", tool="protenix"))])
    # a repatriated manifest sits somewhere under the run dir
    md = tmp_path / "run" / "folds_cif"; md.mkdir(parents=True)
    (md / "c1_msa_manifest.json").write_text(json.dumps(
        {"A": {"got_msa": True}, "B": {"got_msa": True}, "C": {"got_msa": True}, "D": {"got_msa": False}}))
    out = _run(at.render_final_report.handler({"run_dir": rd}))
    html = Path(out["structuredContent"]["report_path"]).read_text()
    assert "MSA colab_cpu (3/4 chains)" in html   # not "MSA-free", the foldjob prep value
```

- [ ] **Step 6: Run to verify it fails**

Run: `./.venv/bin/python -m pytest tests/test_agent_tools.py::test_render_final_report_reads_msa_manifest -q`
Expected: FAIL (report shows "MSA-free"; the manifest is not read yet).

- [ ] **Step 7: Wire the manifest into render_final_report**

In `src/rep2struct/agent_tools.py`, in `render_final_report`, replace the `msa_basis` line

```python
    msa_basis = {j["clonotype_id"]: j.get("msa_basis") for j in fjs}
```

with (also import `msa_basis_from_manifest` at the top of the function body):

```python
    import json as _json
    from .report import msa_basis_from_manifest
    manifests = {}
    for mp in Path(args["run_dir"]).rglob("*_msa_manifest.json"):
        cid = mp.name[:-len("_msa_manifest.json")]
        try:
            manifests[cid] = _json.loads(mp.read_text())
        except (ValueError, OSError):
            pass
    msa_basis = {
        j["clonotype_id"]: (msa_basis_from_manifest(manifests[j["clonotype_id"]])
                            if j["clonotype_id"] in manifests else j.get("msa_basis"))
        for j in fjs
    }
```

- [ ] **Step 8: Run the full suite**

Run: `./.venv/bin/python -m pytest -q`
Expected: PASS (137 tests: 133 after Task 1, plus 4 here).

- [ ] **Step 9: Commit**

```bash
git add src/rep2struct/report.py src/rep2struct/agent_tools.py tests/test_report.py tests/test_agent_tools.py
git commit -m "feat(msa): report the real per-clonotype MSA basis from the manifest

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: Live validation (the success criterion)

No code file. This is the proof the slice worked, run on the coquetlab.sic A100 Colab as in the loop-closing session. It also settles the one residual interface unknown.

- [ ] **Step 1: Build the MSA notebook for the flu M1 clonotype**

Reuse the existing run dir (it already has ingest/annotate/foldjobs for `data/validation_tcrdock_classI.csv`, clonotype `9ab6b3bfa998`). Regenerate its Protenix notebook so it now carries the MSA cell:

```bash
./.venv/bin/python - <<'PY'
import asyncio, sys; sys.path.insert(0, "src")
from rep2struct import agent_tools as at
RD = "/private/tmp/claude-501/-Users-fzd181/675f8148-b48a-46bd-81f3-c0d57eae74fe/scratchpad/live_run"
print(asyncio.run(at.build_fold_notebook.handler(
    {"run_dir": RD, "clonotype_id": "9ab6b3bfa998", "tool": "protenix"}))["structuredContent"])
PY
```

Expected: a notebook path under `<RD>/notebooks/9ab6b3bfa998_protenix.ipynb`.

- [ ] **Step 2: Drive it on the A100 Colab**

Upload and run it on the authenticated coquetlab.sic Colab (same CDP-attached Chrome path used before). Watch the MSA cell print `MSA_OK` for the 4 protein chains (A, B, C, D) and the fold cell print `FOLDED ... 5 models` for cognate and scramble.

- [ ] **Step 3: Verify --use_msa false honored the provided MSA**

Confirm the Protenix run actually used the a3m (the fold log references the MSA / does not fall back to single-sequence). If the fold ignored `unpairedMsaPath` because of `--use_msa false`, edit the fold cell in `_protenix_notebook` to drop `--use_msa false` (Protenix then uses provided paths and searches none for chains that have them), add a test asserting the flag is gone, re-run the suite, commit, and re-drive.

- [ ] **Step 4: Repatriate and run QC**

Download `protenix_folds.zip`, unzip into the run dir, then:

```bash
./.venv/bin/python /private/tmp/claude-501/-Users-fzd181/675f8148-b48a-46bd-81f3-c0d57eae74fe/scratchpad/real_qc.py
```

- [ ] **Step 5: Confirm the success criterion**

The cognate must now dock: median CDR3-peptide contact > 0 and > its scramble, versus today's MSA-free run where 4 of 5 poses had zero contact. Confirm the report shows `MSA colab_cpu (4/4 chains)` for the clonotype. Record the numbers in `docs/fold_qc_results.md`.

---

## Self-Review

**Spec coverage:**
- In-notebook Colab CPU MSA cell (spec 2, 3): Task 1.
- Unpaired, per protein chain, dedup, peptide exclusion by len<20 (spec 3, scope): Task 1.
- `unpairedMsaPath` injection into the JSON (spec interface, 4): Task 1 (mutates INPUTS in place; write cell unchanged).
- Graceful per-chain fallback, never fatal (spec 3): Task 1 (`MSA_FAIL` log, no path).
- `--use_msa false` kept, verified live, drop-flag fallback (spec 4, risks): Task 3 Step 3.
- Manifest produced (spec 3) and report honesty per clonotype (spec, report honesty, chosen for v1): Task 1 (writes `{cid}_msa_manifest.json`) + Task 2 (reads and renders it).
- msa.py unchanged (spec msa.py): honored, no task touches it.
- Live success criterion, flu M1 docks (spec testing): Task 3.

**Placeholder scan:** none. Every code step shows complete code; the one interface unknown (`--use_msa false` honoring paths) is a live verification step with a concrete contingency, not a deferred TODO.

**Type consistency:** manifest shape `{chain_id: {"got_msa": bool}}` is written in Task 1 and read by `msa_basis_from_manifest` in Task 2. The basis token `"colab_cpu:k/n"` is produced by `msa_basis_from_manifest` and parsed by `_msa_note`, both in Task 2. `render_final_report` feeds `render_report(msa_basis=...)`, whose `_msa_note` consumes the token. Consistent.
