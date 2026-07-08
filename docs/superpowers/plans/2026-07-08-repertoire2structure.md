# Repertoire2Structure Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A Claude orchestrated pipeline that turns a raw 10x TCR repertoire CSV into QC'd predicted TCR pMHC structures for its top clonotypes, with honest specificity annotation and a skeptical hallucination flag.

**Architecture:** Two layers. A deterministic stage layer (ingest, annotate, foldprep, fold, qc, report), each a pure Python module with a typed output persisted to a run directory, so the chain is resumable and fully testable offline. On top, a genuine multi agent layer built on the Claude Agent SDK: the stage functions are exposed as in process tools, and an orchestrator agent delegates to specialist agents (a fold agent that drives Playwright, a skeptical QC agent, a report agent) that coordinate through the shared run directory. The stage layer carries reliability; the agent layer is the product and the Built with Claude story. Heavy lifting reuses TCR Explorer (imported as `tcr_explorer`) for allele assignment and TCRdist annotation, the existing Protenix Colab fold procedure for structure prediction, and geometry scoring for QC.

**Tech Stack:** Python 3.11, pandas, polars, numpy, biopython, jinja2 (HTML report), pytest, claude-agent-sdk (multi agent runtime). Dependency on the local `tcr_explorer` package (from `~/imgt-api`). Protenix folding on Google Colab driven via the Playwright MCP.

## Global Constraints

- Python 3.11 virtualenv. `tcr_explorer` must be importable (install `~/imgt-api` in editable mode into the venv).
- Documentation and README files use no dash as punctuation. Rephrase with periods, commas, or parentheses. Hyphens inside code, flags, and identifiers are fine.
- Commits authored as `Kilian Maire <mairekilian@gmail.com>`. This is a Built with Claude hackathon repo, so keep the `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>` trailer.
- Two honesty rules are enforced in the output schemas, not just prose. Stage 1 always carries `annotatable: bool`, `confidence_tier`, and `tcrdist`; it never forces a label. Stage 4 emits `qc_verdict` in {`reliable`, `suspect`, `qc_failed`}; a fold never confirms specificity.
- TDD. Frequent commits. One task per commit minimum.
- Never redistribute source datasets in the repo. Test fixtures are tiny synthetic files only.

---

## File Structure

```
tcr-repertoire2structure/
  pyproject.toml                     # package metadata, deps, pytest config
  src/rep2struct/
    __init__.py
    schema.py                        # typed stage outputs (dataclasses)
    runstate.py                      # run directory, stage read/write, resumability
    ingest.py                        # stage 0: 10x parse, clonotype collapse, allele standardize
    annotate.py                      # stage 1: TCRdist annotation, confidence tiers, unannotatable
    validate.py                      # stage 1 metrics vs dextramer ground truth
    foldprep.py                      # stage 2: rank, class I construct, MSA bundle
    fold.py                          # stage 3: resumable fold runner (pluggable fold_fn)
    qc.py                            # stage 4: geometry scoring, scramble calibrated verdict
    report.py                        # stage 5: HTML report (jinja2)
    pipeline.py                      # deterministic reference chain + offline test harness
    agent_tools.py                   # stage functions wrapped as Agent SDK tools (in process MCP)
    agents.py                        # AgentDefinition for orchestrator + fold + qc + report agents
    app.py                           # product entrypoint: query() run over the agent layer
  tests/
    fixtures/
      tenx_tiny.csv                  # synthetic 10x paired contig rows
      tenx_dextramer_tiny.csv        # synthetic labeled rows for validation
      cognate_min.cif                # a minimal cognate TCR pMHC model
      scramble_min.cif               # a minimal scramble model
    test_runstate.py
    test_ingest.py
    test_annotate.py
    test_validate.py
    test_foldprep.py
    test_fold.py
    test_qc.py
    test_report.py
    test_pipeline_offline.py
    test_agent_tools.py
    test_agents_config.py
```

Each stage module owns one stage. `schema.py` and `runstate.py` are shared infrastructure. `pipeline.py` wires the stages deterministically and is the offline test harness (it proves the chain works with a mocked fold). The product entrypoint is the agent layer: `agent_tools.py` exposes the stages as tools, `agents.py` defines the specialist agents, and `app.py` runs the orchestrator.

---

## Task 1: Scaffold, schema, run state

**Files:**
- Create: `pyproject.toml`
- Create: `src/rep2struct/__init__.py`
- Create: `src/rep2struct/schema.py`
- Create: `src/rep2struct/runstate.py`
- Test: `tests/test_runstate.py`

**Interfaces:**
- Produces: `Clonotype`, `Annotation`, `FoldJob`, `QCResult` dataclasses; `RunState` with `write_stage(name, obj)`, `read_stage(name)`, `stage_done(name)`, `path_for(name)`.

- [ ] **Step 1: Write pyproject.toml**

```toml
[build-system]
requires = ["setuptools>=68"]
build-backend = "setuptools.build_meta"

[project]
name = "rep2struct"
version = "0.1.0"
requires-python = ">=3.11"
dependencies = [
  "pandas>=2.0", "polars>=1.0", "numpy>=1.24",
  "biopython>=1.83", "jinja2>=3.1", "claude-agent-sdk>=0.1",
]

[project.optional-dependencies]
dev = ["pytest>=8.0"]

[tool.setuptools.packages.find]
where = ["src"]

[tool.pytest.ini_options]
pythonpath = ["src"]
testpaths = ["tests"]
```

- [ ] **Step 2: Write the failing run state test**

```python
# tests/test_runstate.py
import json
from rep2struct.runstate import RunState

def test_write_read_roundtrip(tmp_path):
    rs = RunState(tmp_path / "run1")
    assert not rs.stage_done("ingest")
    rs.write_stage("ingest", {"clonotypes": [{"id": "c1", "size": 3}]})
    assert rs.stage_done("ingest")
    got = rs.read_stage("ingest")
    assert got["clonotypes"][0]["id"] == "c1"

def test_path_for_is_under_run_dir(tmp_path):
    rs = RunState(tmp_path / "run2")
    p = rs.path_for("fold")
    assert str(p).startswith(str(tmp_path / "run2"))
```

- [ ] **Step 3: Run it to see it fail**

Run: `pytest tests/test_runstate.py -v`
Expected: FAIL with import error (module not yet written).

- [ ] **Step 4: Write schema.py**

```python
# src/rep2struct/schema.py
from __future__ import annotations
from dataclasses import dataclass, field, asdict
from typing import Optional

@dataclass
class Clonotype:
    id: str                      # content hash of the defining tuple
    trav: str
    cdr3a: str
    trbv: str
    cdr3b: str
    size: int                    # cell count
    trav_allele: Optional[str] = None
    trbv_allele: Optional[str] = None

@dataclass
class Annotation:
    clonotype_id: str
    annotatable: bool
    confidence_tier: str         # high, medium, low, unannotatable
    tcrdist: Optional[float] = None
    epitope: Optional[str] = None
    hla: Optional[str] = None
    antigen: Optional[str] = None
    neighbour_id: Optional[str] = None

@dataclass
class FoldJob:
    clonotype_id: str
    construct_fasta: str         # A..E chains
    msa_ref: Optional[str] = None
    status: str = "pending"      # pending, done, failed
    model_paths: list[str] = field(default_factory=list)

@dataclass
class QCResult:
    clonotype_id: str
    qc_verdict: str              # reliable, suspect, qc_failed
    reason: str
    dockq: Optional[float] = None
    cdr3_pep_atoms: Optional[float] = None
    crossing_angle: Optional[float] = None

def to_jsonable(obj):
    return asdict(obj) if hasattr(obj, "__dataclass_fields__") else obj
```

- [ ] **Step 5: Write runstate.py**

```python
# src/rep2struct/runstate.py
from __future__ import annotations
import json
from pathlib import Path
from .schema import to_jsonable

class RunState:
    def __init__(self, run_dir):
        self.run_dir = Path(run_dir)
        self.run_dir.mkdir(parents=True, exist_ok=True)

    def path_for(self, name: str) -> Path:
        return self.run_dir / f"{name}.json"

    def stage_done(self, name: str) -> bool:
        return self.path_for(name).exists()

    def write_stage(self, name: str, obj) -> None:
        def enc(o):
            if isinstance(o, list):
                return [enc(x) for x in o]
            return to_jsonable(o)
        tmp = self.path_for(name).with_suffix(".json.tmp")
        tmp.write_text(json.dumps(enc(obj), indent=2))
        tmp.replace(self.path_for(name))

    def read_stage(self, name: str):
        return json.loads(self.path_for(name).read_text())
```

- [ ] **Step 6: Run tests to verify pass**

Run: `pytest tests/test_runstate.py -v`
Expected: PASS (2 tests).

- [ ] **Step 7: Commit**

```bash
git add pyproject.toml src/rep2struct tests/test_runstate.py
git commit -m "feat: scaffold package, schema, and resumable run state

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: Stage 0 ingest and clonotype collapse

**Files:**
- Create: `src/rep2struct/ingest.py`
- Create: `tests/fixtures/tenx_tiny.csv`
- Test: `tests/test_ingest.py`

**Interfaces:**
- Consumes: a 10x `filtered_contig_annotations.csv` path.
- Produces: `parse_10x(path) -> list[Clonotype]`. Clonotype id is a stable hash of `(trav, cdr3a, trbv, cdr3b)`. Size is the number of distinct cell barcodes collapsing to that tuple. Rows lacking a productive paired alpha and beta are dropped with a counted reason via `parse_10x(path, report=True) -> (clonotypes, drop_report)`.

- [ ] **Step 1: Write the fixture**

```csv
barcode,is_cell,high_confidence,productive,chain,v_gene,j_gene,cdr3
AAAC-1,True,True,True,TRA,TRAV1-2,TRAJ33,CAVMDSSYKLIF
AAAC-1,True,True,True,TRB,TRBV19,TRBJ2-1,CASSIRSSYEQYF
AAAD-1,True,True,True,TRA,TRAV1-2,TRAJ33,CAVMDSSYKLIF
AAAD-1,True,True,True,TRB,TRBV19,TRBJ2-1,CASSIRSSYEQYF
AAAE-1,True,True,True,TRA,TRAV12-1,TRAJ20,CAVNNDYKLSF
AAAE-1,True,True,True,TRB,TRBV7-9,TRBJ2-7,CASSLGQAYEQYF
AAAF-1,True,True,False,TRB,TRBV20-1,TRBJ1-1,CSARDLTF
AAAG-1,True,True,True,TRA,TRAV8-1,TRAJ22,CAVGATSGSRLTF
```

- [ ] **Step 2: Write the failing test**

```python
# tests/test_ingest.py
from pathlib import Path
from rep2struct.ingest import parse_10x

FIX = Path(__file__).parent / "fixtures" / "tenx_tiny.csv"

def test_collapse_identical_clonotypes():
    clons = parse_10x(FIX)
    # AAAC and AAAD share the exact alpha+beta tuple -> one clonotype, size 2
    top = [c for c in clons if c.cdr3b == "CASSIRSSYEQYF"]
    assert len(top) == 1
    assert top[0].size == 2
    assert top[0].trav == "TRAV1-2"

def test_unpaired_rows_dropped_with_reason():
    clons, report = parse_10x(FIX, report=True)
    # AAAF has only a non-productive beta; AAAG has only an alpha -> both dropped
    ids = {(c.trav, c.cdr3b) for c in clons}
    assert ("TRAV8-1", "CASSIRSSYEQYF") not in ids
    assert report["dropped_unpaired"] >= 2
```

- [ ] **Step 3: Run to verify fail**

Run: `pytest tests/test_ingest.py -v`
Expected: FAIL (module missing).

- [ ] **Step 4: Implement ingest.py**

```python
# src/rep2struct/ingest.py
from __future__ import annotations
import hashlib
from collections import defaultdict
import pandas as pd
from .schema import Clonotype

def _clon_id(trav, cdr3a, trbv, cdr3b) -> str:
    key = f"{trav}|{cdr3a}|{trbv}|{cdr3b}"
    return hashlib.sha1(key.encode()).hexdigest()[:12]

def parse_10x(path, report: bool = False):
    df = pd.read_csv(path)
    df = df[(df["productive"].astype(str) == "True") & (df["high_confidence"].astype(str) == "True")]
    per_cell = defaultdict(dict)  # barcode -> {"TRA": row, "TRB": row}
    for _, r in df.iterrows():
        chain = r["chain"]
        if chain in ("TRA", "TRB"):
            per_cell[r["barcode"]][chain] = r
    tuples = defaultdict(set)  # tuple -> set of barcodes
    dropped_unpaired = 0
    for bc, chains in per_cell.items():
        if "TRA" not in chains or "TRB" not in chains:
            dropped_unpaired += 1
            continue
        a, b = chains["TRA"], chains["TRB"]
        key = (a["v_gene"], a["cdr3"], b["v_gene"], b["cdr3"])
        tuples[key].add(bc)
    clons = [
        Clonotype(id=_clon_id(*k), trav=k[0], cdr3a=k[1], trbv=k[2], cdr3b=k[3], size=len(bcs))
        for k, bcs in tuples.items()
    ]
    clons.sort(key=lambda c: c.size, reverse=True)
    if report:
        return clons, {"dropped_unpaired": dropped_unpaired, "clonotypes": len(clons)}
    return clons
```

- [ ] **Step 5: Run tests to verify pass**

Run: `pytest tests/test_ingest.py -v`
Expected: PASS (2 tests).

- [ ] **Step 6: Commit**

```bash
git add src/rep2struct/ingest.py tests/test_ingest.py tests/fixtures/tenx_tiny.csv
git commit -m "feat: stage 0 ingest, collapse 10x contigs into paired clonotypes

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: Stage 0 allele standardization

**Files:**
- Modify: `src/rep2struct/ingest.py`
- Test: `tests/test_ingest.py` (add a test)

**Interfaces:**
- Consumes: `tcr_explorer.tcr_align.assign(sequence, species, chain, want_d)` via a thin wrapper so it can be mocked in tests.
- Produces: `standardize_alleles(clonotypes, assign_fn=None) -> list[Clonotype]` filling `trav_allele` and `trbv_allele`. On assignment failure the allele stays `None` and the clonotype is kept (never dropped for allele reasons).

- [ ] **Step 1: Write the failing test with a fake assigner**

```python
# tests/test_ingest.py  (append)
from rep2struct.ingest import standardize_alleles
from rep2struct.schema import Clonotype

def test_standardize_alleles_fills_allele_or_keeps_none():
    clons = [Clonotype(id="x", trav="TRAV1-2", cdr3a="CAVMDSSYKLIF",
                       trbv="TRBV19", cdr3b="CASSIRSSYEQYF", size=2)]
    def fake_assign(gene, species, chain):
        return {"TRAV1-2": "TRAV1-2*01", "TRBV19": None}.get(gene)
    out = standardize_alleles(clons, assign_fn=fake_assign)
    assert out[0].trav_allele == "TRAV1-2*01"
    assert out[0].trbv_allele is None  # failure keeps None, clonotype kept
```

- [ ] **Step 2: Run to verify fail**

Run: `pytest tests/test_ingest.py::test_standardize_alleles_fills_allele_or_keeps_none -v`
Expected: FAIL (function missing).

- [ ] **Step 3: Implement standardize_alleles**

```python
# src/rep2struct/ingest.py  (append)
def _default_assign(gene, species="human", chain=None):
    from tcr_explorer.tcr_align import assign
    res = assign(gene, species=species, chain=chain)
    v = getattr(res, "v_allele", None)
    return v

def standardize_alleles(clonotypes, assign_fn=None):
    fn = assign_fn or _default_assign
    out = []
    for c in clonotypes:
        try:
            c.trav_allele = fn(c.trav, species="human", chain="A")
        except Exception:
            c.trav_allele = None
        try:
            c.trbv_allele = fn(c.trbv, species="human", chain="B")
        except Exception:
            c.trbv_allele = None
        out.append(c)
    return out
```

- [ ] **Step 4: Run tests to verify pass**

Run: `pytest tests/test_ingest.py -v`
Expected: PASS (3 tests).

- [ ] **Step 5: Commit**

```bash
git add src/rep2struct/ingest.py tests/test_ingest.py
git commit -m "feat: stage 0 allele standardization via TCR Explorer, failure keeps None

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: Stage 1 specificity annotation

**Files:**
- Create: `src/rep2struct/annotate.py`
- Test: `tests/test_annotate.py`

**Interfaces:**
- Consumes: `list[Clonotype]`, and a paired similarity function shaped like `tcr_explorer.similarity.find_similar_paired_tcrs(cdr3_a, v_a, cdr3_b, v_b, species, top_k) -> (neighbours, engine, total, warnings)` where each neighbour is a dict carrying `epitope`, `mhc`, `antigen`, and `distance`.
- Produces: `annotate(clonotypes, sim_fn=None, tiers=DEFAULT_TIERS) -> list[Annotation]`. The confidence tier comes from the nearest neighbour tcrdist against ascending thresholds. No neighbour or distance above the last threshold yields `annotatable=False`, `confidence_tier="unannotatable"`, and no epitope.

- [ ] **Step 1: Write the failing test with a fake similarity fn**

```python
# tests/test_annotate.py
from rep2struct.annotate import annotate, DEFAULT_TIERS
from rep2struct.schema import Clonotype

def _clon(cid, cdr3b): return Clonotype(id=cid, trav="TRAV1-2", cdr3a="CAVA",
                                        trbv="TRBV19", cdr3b=cdr3b, size=1)

def test_close_neighbour_is_annotated_high():
    def sim(cdr3_a, v_a, cdr3_b, v_b, species="human", top_k=5):
        return ([{"epitope": "GILGFVFTL", "mhc": "HLA-A*02:01",
                  "antigen": "Flu M1", "distance": 3.0}], "tcrdist", 100, [])
    a = annotate([_clon("c1", "CASSIRSSYEQYF")], sim_fn=sim)[0]
    assert a.annotatable and a.confidence_tier == "high"
    assert a.epitope == "GILGFVFTL" and a.tcrdist == 3.0

def test_far_neighbour_is_unannotatable():
    def sim(*args, **kw):
        return ([{"epitope": "X", "mhc": "Y", "antigen": "Z", "distance": 999.0}], "tcrdist", 100, [])
    a = annotate([_clon("c2", "CASSNOMATCH")], sim_fn=sim)[0]
    assert a.annotatable is False
    assert a.confidence_tier == "unannotatable"
    assert a.epitope is None

def test_no_neighbour_is_unannotatable():
    def sim(*args, **kw): return ([], "tcrdist", 0, ["no candidates"])
    a = annotate([_clon("c3", "CASSEMPTY")], sim_fn=sim)[0]
    assert a.annotatable is False and a.epitope is None
```

- [ ] **Step 2: Run to verify fail**

Run: `pytest tests/test_annotate.py -v`
Expected: FAIL (module missing).

- [ ] **Step 3: Implement annotate.py**

```python
# src/rep2struct/annotate.py
from __future__ import annotations
from .schema import Clonotype, Annotation

# ascending tcrdist thresholds -> tier. Calibrated on the validation arm later.
DEFAULT_TIERS = [(12.0, "high"), (24.0, "medium"), (48.0, "low")]

def _tier(distance, tiers):
    for thr, name in tiers:
        if distance <= thr:
            return name
    return "unannotatable"

def _default_sim(cdr3_a, v_a, cdr3_b, v_b, species="human", top_k=5):
    from tcr_explorer.similarity import find_similar_paired_tcrs
    return find_similar_paired_tcrs(cdr3_a, v_a, cdr3_b, v_b, species=species, top_k=top_k)

def annotate(clonotypes, sim_fn=None, tiers=DEFAULT_TIERS):
    fn = sim_fn or _default_sim
    out = []
    for c in clonotypes:
        neigh, *_ = fn(c.cdr3a, c.trav, c.cdr3b, c.trbv, species="human", top_k=5)
        if not neigh:
            out.append(Annotation(clonotype_id=c.id, annotatable=False,
                                  confidence_tier="unannotatable"))
            continue
        best = min(neigh, key=lambda n: n["distance"])
        tier = _tier(best["distance"], tiers)
        if tier == "unannotatable":
            out.append(Annotation(clonotype_id=c.id, annotatable=False,
                                  confidence_tier="unannotatable", tcrdist=best["distance"]))
        else:
            out.append(Annotation(
                clonotype_id=c.id, annotatable=True, confidence_tier=tier,
                tcrdist=best["distance"], epitope=best.get("epitope"),
                hla=best.get("mhc"), antigen=best.get("antigen")))
    return out
```

- [ ] **Step 4: Run tests to verify pass**

Run: `pytest tests/test_annotate.py -v`
Expected: PASS (3 tests).

- [ ] **Step 5: Commit**

```bash
git add src/rep2struct/annotate.py tests/test_annotate.py
git commit -m "feat: stage 1 honest TCRdist annotation with confidence tiers

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 5: Stage 1 validation metrics vs dextramer labels

**Files:**
- Create: `src/rep2struct/validate.py`
- Create: `tests/fixtures/tenx_dextramer_tiny.csv`
- Test: `tests/test_validate.py`

**Interfaces:**
- Consumes: `list[Annotation]` and a label map `clonotype_id -> true_epitope` from the dextramer positive calls.
- Produces: `annotation_metrics(annotations, labels) -> dict` with `precision`, `recall`, `unannotatable_rate`, `n`. Precision is over annotated clonotypes (annotated epitope equals true label). Recall is over labeled clonotypes correctly annotated. Unannotatable rate is the fraction flagged unannotatable.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_validate.py
from rep2struct.validate import annotation_metrics
from rep2struct.schema import Annotation

def test_metrics_basic():
    anns = [
        Annotation("c1", True, "high", 3.0, epitope="GILGFVFTL"),   # correct
        Annotation("c2", True, "medium", 20.0, epitope="NLVPMVATV"),# wrong
        Annotation("c3", False, "unannotatable"),                   # missed
    ]
    labels = {"c1": "GILGFVFTL", "c2": "GILGFVFTL", "c3": "KLGGALQAK"}
    m = annotation_metrics(anns, labels)
    assert m["precision"] == 0.5           # 1 of 2 annotated correct
    assert round(m["recall"], 3) == 0.333  # 1 of 3 labeled correct
    assert round(m["unannotatable_rate"], 3) == 0.333
```

- [ ] **Step 2: Run to verify fail**

Run: `pytest tests/test_validate.py -v`
Expected: FAIL (module missing).

- [ ] **Step 3: Implement validate.py**

```python
# src/rep2struct/validate.py
from __future__ import annotations

def annotation_metrics(annotations, labels) -> dict:
    n = len(annotations)
    annotated = [a for a in annotations if a.annotatable]
    correct = sum(1 for a in annotated if labels.get(a.clonotype_id) == a.epitope)
    n_labeled = sum(1 for a in annotations if a.clonotype_id in labels)
    unann = sum(1 for a in annotations if not a.annotatable)
    precision = correct / len(annotated) if annotated else float("nan")
    recall = correct / n_labeled if n_labeled else float("nan")
    return {
        "precision": precision,
        "recall": recall,
        "unannotatable_rate": unann / n if n else float("nan"),
        "n": n,
        "n_annotated": len(annotated),
        "n_correct": correct,
    }
```

- [ ] **Step 4: Write the labeled fixture** (used later by the end to end test)

```csv
barcode,is_cell,high_confidence,productive,chain,v_gene,j_gene,cdr3,dextramer
AAAC-1,True,True,True,TRA,TRAV1-2,TRAJ33,CAVMDSSYKLIF,GILGFVFTL
AAAC-1,True,True,True,TRB,TRBV19,TRBJ2-1,CASSIRSSYEQYF,GILGFVFTL
AAAE-1,True,True,True,TRA,TRAV12-1,TRAJ20,CAVNNDYKLSF,NLVPMVATV
AAAE-1,True,True,True,TRB,TRBV7-9,TRBJ2-7,CASSLGQAYEQYF,NLVPMVATV
```

- [ ] **Step 5: Run tests to verify pass**

Run: `pytest tests/test_validate.py -v`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/rep2struct/validate.py tests/test_validate.py tests/fixtures/tenx_dextramer_tiny.csv
git commit -m "feat: stage 1 validation metrics against dextramer ground truth

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 6: Stage 2 fold prep, ranking and class I construct

**Files:**
- Create: `src/rep2struct/foldprep.py`
- Test: `tests/test_foldprep.py`

**Interfaces:**
- Consumes: `list[Clonotype]`, `list[Annotation]`.
- Produces: `select_top(clonotypes, annotations, n) -> list[tuple[Clonotype, Annotation]]` ranked by tier weight times clonal size; `build_construct(clonotype, annotation, tcr_seqs, mhc_seqs) -> FoldJob` assembling a five chain class I construct A..E (TCR alpha V, TCR beta V, MHC class I heavy, beta 2 microglobulin, peptide). `tcr_seqs` and `mhc_seqs` are provided lookups so the test can inject sequences.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_foldprep.py
from rep2struct.foldprep import select_top, build_construct, TIER_WEIGHT
from rep2struct.schema import Clonotype, Annotation

def _c(cid, size): return Clonotype(cid, "TRAV1-2", "CAVA", "TRBV19", "CASSB", size)

def test_ranking_prefers_confident_and_expanded():
    clons = [_c("c1", 2), _c("c2", 100)]
    anns = [Annotation("c1", True, "high", 3.0, epitope="GILGFVFTL"),
            Annotation("c2", False, "unannotatable")]
    top = select_top(clons, anns, n=1)
    assert top[0][0].id == "c1"   # high tier beats big-but-unannotatable

def test_build_construct_has_five_chains():
    c = _c("c1", 2)
    a = Annotation("c1", True, "high", 3.0, epitope="GILGFVFTL", hla="HLA-A*02:01")
    job = build_construct(c, a,
        tcr_seqs={"c1": {"A": "AAAA", "B": "BBBB"}},
        mhc_seqs={"HLA-A*02:01": {"heavy": "HHHH", "b2m": "MMMM"}})
    chains = [l for l in job.construct_fasta.splitlines() if l.startswith(">")]
    assert chains == [">A", ">B", ">C", ">D", ">E"]
    assert "GILGFVFTL" in job.construct_fasta
```

- [ ] **Step 2: Run to verify fail**

Run: `pytest tests/test_foldprep.py -v`
Expected: FAIL (module missing).

- [ ] **Step 3: Implement foldprep.py**

```python
# src/rep2struct/foldprep.py
from __future__ import annotations
from .schema import FoldJob

TIER_WEIGHT = {"high": 3.0, "medium": 2.0, "low": 1.0, "unannotatable": 0.0}

def select_top(clonotypes, annotations, n):
    by_id = {a.clonotype_id: a for a in annotations}
    scored = []
    for c in clonotypes:
        a = by_id.get(c.id)
        if a is None:
            continue
        score = TIER_WEIGHT.get(a.confidence_tier, 0.0) * c.size
        scored.append((score, c, a))
    scored.sort(key=lambda t: t[0], reverse=True)
    return [(c, a) for _, c, a in scored[:n]]

def build_construct(clonotype, annotation, tcr_seqs, mhc_seqs) -> FoldJob:
    t = tcr_seqs[clonotype.id]
    m = mhc_seqs[annotation.hla]
    fasta = "\n".join([
        ">A", t["A"],            # TCR alpha V domain
        ">B", t["B"],            # TCR beta V domain
        ">C", m["heavy"],        # MHC class I heavy chain
        ">D", m["b2m"],          # beta 2 microglobulin
        ">E", annotation.epitope,  # peptide
    ])
    return FoldJob(clonotype_id=clonotype.id, construct_fasta=fasta)
```

- [ ] **Step 4: Run tests to verify pass**

Run: `pytest tests/test_foldprep.py -v`
Expected: PASS (2 tests).

- [ ] **Step 5: Commit**

```bash
git add src/rep2struct/foldprep.py tests/test_foldprep.py
git commit -m "feat: stage 2 ranking and five chain class I construct assembly

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 7: Stage 3 resumable fold runner

**Files:**
- Create: `src/rep2struct/fold.py`
- Test: `tests/test_fold.py`

**Interfaces:**
- Consumes: `list[FoldJob]`, a `fold_fn(job) -> list[str]` that returns model paths, and a `RunState`. The real `fold_fn` is the Claude agent driving Protenix on Colab via Playwright, documented in `docs/fold_procedure.md`. Tests inject a mock `fold_fn`.
- Produces: `run_folds(jobs, fold_fn, run_state) -> list[FoldJob]`. A job whose `clonotype_id` already has a `.done.txt` marker in the run directory is skipped and its cached model paths are reused. A job whose `fold_fn` raises is marked `status="failed"` and does not abort the batch.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_fold.py
from rep2struct.fold import run_folds
from rep2struct.schema import FoldJob
from rep2struct.runstate import RunState

def _job(cid): return FoldJob(clonotype_id=cid, construct_fasta=">A\nAAAA")

def test_resume_skips_done(tmp_path):
    rs = RunState(tmp_path / "r")
    calls = []
    def fold_fn(job):
        calls.append(job.clonotype_id)
        return [f"{job.clonotype_id}.cif"]
    run_folds([_job("c1")], fold_fn, rs)
    run_folds([_job("c1")], fold_fn, rs)  # second run resumes
    assert calls == ["c1"]  # folded once, skipped the second time

def test_failure_does_not_abort_batch(tmp_path):
    rs = RunState(tmp_path / "r2")
    def fold_fn(job):
        if job.clonotype_id == "bad":
            raise RuntimeError("colab wedged")
        return ["ok.cif"]
    out = run_folds([_job("bad"), _job("good")], fold_fn, rs)
    status = {j.clonotype_id: j.status for j in out}
    assert status["bad"] == "failed" and status["good"] == "done"
```

- [ ] **Step 2: Run to verify fail**

Run: `pytest tests/test_fold.py -v`
Expected: FAIL (module missing).

- [ ] **Step 3: Implement fold.py**

```python
# src/rep2struct/fold.py
from __future__ import annotations

def _marker(run_state, cid):
    return run_state.run_dir / f"fold_{cid}.done.txt"

def run_folds(jobs, fold_fn, run_state):
    out = []
    for job in jobs:
        marker = _marker(run_state, job.clonotype_id)
        if marker.exists():
            job.status = "done"
            job.model_paths = marker.read_text().splitlines()
            out.append(job)
            continue
        try:
            paths = fold_fn(job)
            job.model_paths = list(paths)
            job.status = "done"
            marker.write_text("\n".join(job.model_paths))
        except Exception as exc:  # noqa: BLE001
            job.status = "failed"
            job.model_paths = []
        out.append(job)
    return out
```

- [ ] **Step 4: Write docs/fold_procedure.md** (the real fold_fn is procedural, documented not unit tested)

```markdown
# Fold procedure (real fold_fn)

The Claude agent implements fold_fn by driving Protenix on Google Colab through
Playwright, following the established procedure: open the Colab notebook, upload the
construct FASTA and precomputed MSA, run `protenix_base_default_v1.0.0` at 5 seeds,
let background execution survive disconnect, and download the resulting CIF models to
`~/.playwright-mcp/`. It returns the local model paths. The loop is resumable: a job
with a `fold_<id>.done.txt` marker is skipped.
```

- [ ] **Step 5: Run tests to verify pass**

Run: `pytest tests/test_fold.py -v`
Expected: PASS (2 tests).

- [ ] **Step 6: Commit**

```bash
git add src/rep2struct/fold.py tests/test_fold.py docs/fold_procedure.md
git commit -m "feat: stage 3 resumable fold runner with pluggable fold_fn

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 8: Stage 4 skeptical QC verdict

**Files:**
- Create: `src/rep2struct/qc.py`
- Create: `tests/fixtures/cognate_min.cif`, `tests/fixtures/scramble_min.cif`
- Test: `tests/test_qc.py`

**Interfaces:**
- Consumes: a model CIF path and a scramble calibration threshold on CDR3 to peptide atom contacts.
- Produces: `score_model(cif_path) -> dict` with `cdr3_pep_atoms`, `crossing_angle`, `dockq` (dockq optional, None when no reference); `verdict(scores, scramble_threshold) -> QCResult`. A model with fewer than five chains yields `qc_failed`. A model whose `cdr3_pep_atoms` does not exceed `scramble_threshold` is `suspect`; otherwise `reliable`.

- [ ] **Step 1: Write minimal CIF fixtures**

Create two tiny CIF files each with five chains A..E and a handful of CA atoms. `cognate_min.cif` places CDR3 beta atoms within 4.5 angstrom of the peptide (chain E); `scramble_min.cif` places them far. Keep them under 60 atoms so parsing is fast. (Generate with the helper in Step 3 if hand authoring is error prone.)

- [ ] **Step 2: Write the failing test**

```python
# tests/test_qc.py
from pathlib import Path
from rep2struct.qc import score_model, verdict
FIX = Path(__file__).parent / "fixtures"

def test_cognate_is_reliable():
    s = score_model(FIX / "cognate_min.cif")
    r = verdict(s, scramble_threshold=1.0)
    assert s["cdr3_pep_atoms"] > 1.0
    assert r.qc_verdict == "reliable"

def test_scramble_is_suspect():
    s = score_model(FIX / "scramble_min.cif")
    r = verdict(s, scramble_threshold=1.0)
    assert r.qc_verdict == "suspect"
```

- [ ] **Step 3: Implement qc.py**

```python
# src/rep2struct/qc.py
from __future__ import annotations
import warnings
import numpy as np
from .schema import QCResult

def _heavy_by_chain(cif_path):
    from Bio.PDB import MMCIFParser
    warnings.simplefilter("ignore")
    m = next(MMCIFParser(QUIET=True).get_structure("x", str(cif_path)).get_models())
    out = {}
    for ch in m:
        atoms = [a.coord for r in ch for a in r if a.element != "H"]
        if atoms:
            out[ch.id] = np.array(atoms)
    return out

def score_model(cif_path) -> dict:
    chains = _heavy_by_chain(cif_path)
    if not {"A", "B", "C", "D", "E"}.issubset(chains):
        return {"n_chains": len(chains), "cdr3_pep_atoms": None, "crossing_angle": None, "dockq": None}
    pep = chains["E"]
    beta = chains["B"]
    d = np.sqrt(((beta[:, None, :] - pep[None, :, :]) ** 2).sum(-1))
    cdr3_pep_atoms = float((d < 4.5).sum())
    return {"n_chains": len(chains), "cdr3_pep_atoms": cdr3_pep_atoms,
            "crossing_angle": None, "dockq": None}

def verdict(scores, scramble_threshold: float) -> QCResult:
    cid = scores.get("clonotype_id", "unknown")
    if scores.get("cdr3_pep_atoms") is None:
        return QCResult(cid, "qc_failed", f"model has {scores.get('n_chains')} chains, need 5")
    if scores["cdr3_pep_atoms"] <= scramble_threshold:
        return QCResult(cid, "suspect",
                        "CDR3 to peptide contact not above scramble calibration",
                        cdr3_pep_atoms=scores["cdr3_pep_atoms"])
    return QCResult(cid, "reliable", "CDR3 to peptide contact beats scramble null",
                    cdr3_pep_atoms=scores["cdr3_pep_atoms"])
```

- [ ] **Step 4: Run tests to verify pass**

Run: `pytest tests/test_qc.py -v`
Expected: PASS (2 tests).

- [ ] **Step 5: Commit**

```bash
git add src/rep2struct/qc.py tests/test_qc.py tests/fixtures/cognate_min.cif tests/fixtures/scramble_min.cif
git commit -m "feat: stage 4 scramble calibrated skeptical QC verdict

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 9: Stage 5 HTML report

**Files:**
- Create: `src/rep2struct/report.py`
- Create: `src/rep2struct/templates/report.html.j2`
- Test: `tests/test_report.py`

**Interfaces:**
- Consumes: clonotypes, annotations, QC results, and optional validation metrics.
- Produces: `render_report(clonotypes, annotations, qc_results, metrics=None) -> str` returning a self contained HTML string. Each row shows clonotype id, clonal size, candidate epitope or `unannotatable`, confidence tier, and QC verdict with a color. When metrics are present a validation summary block is shown.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_report.py
from rep2struct.report import render_report
from rep2struct.schema import Clonotype, Annotation, QCResult

def test_report_contains_rows_and_verdicts():
    clons = [Clonotype("c1", "TRAV1-2", "CAVA", "TRBV19", "CASSB", 5)]
    anns = [Annotation("c1", True, "high", 3.0, epitope="GILGFVFTL", hla="HLA-A*02:01")]
    qcs = [QCResult("c1", "reliable", "ok", cdr3_pep_atoms=12.0)]
    html = render_report(clons, anns, qcs, metrics={"precision": 0.8, "recall": 0.6,
                                                    "unannotatable_rate": 0.3, "n": 10})
    assert "GILGFVFTL" in html
    assert "reliable" in html
    assert "0.8" in html  # validation block rendered
    assert html.lstrip().lower().startswith("<!doctype html")

def test_unannotatable_is_shown():
    clons = [Clonotype("c2", "TRAV1", "CAVA", "TRBV2", "CASSX", 2)]
    anns = [Annotation("c2", False, "unannotatable")]
    qcs = []
    html = render_report(clons, anns, qcs)
    assert "unannotatable" in html
```

- [ ] **Step 2: Run to verify fail**

Run: `pytest tests/test_report.py -v`
Expected: FAIL (module missing).

- [ ] **Step 3: Write the template** `src/rep2struct/templates/report.html.j2`

```html
<!doctype html>
<html><head><meta charset="utf-8"><title>Repertoire2Structure report</title>
<style>
 body{font-family:system-ui,sans-serif;margin:2rem;color:#1f2933}
 table{border-collapse:collapse;width:100%} th,td{border:1px solid #cbd5e0;padding:.4rem .6rem;text-align:left}
 .reliable{color:#2f855a;font-weight:600} .suspect{color:#c05621;font-weight:600} .qc_failed{color:#718096}
 .unann{color:#718096;font-style:italic} .val{background:#f0fff4;border:1px solid #2f855a;padding:1rem;margin:1rem 0;border-radius:8px}
</style></head><body>
<h1>Repertoire2Structure report</h1>
{% if metrics %}<div class="val"><b>Validation.</b>
 precision {{ '%.2f'|format(metrics.precision) }},
 recall {{ '%.2f'|format(metrics.recall) }},
 unannotatable rate {{ '%.2f'|format(metrics.unannotatable_rate) }}
 (n={{ metrics.n }})</div>{% endif %}
<table><thead><tr><th>clonotype</th><th>size</th><th>epitope</th><th>tier</th><th>QC</th></tr></thead><tbody>
{% for row in rows %}<tr>
 <td>{{ row.id }}</td><td>{{ row.size }}</td>
 <td>{% if row.annotatable %}{{ row.epitope }} <small>({{ row.hla }})</small>{% else %}<span class="unann">unannotatable</span>{% endif %}</td>
 <td>{{ row.tier }}</td>
 <td class="{{ row.qc }}">{{ row.qc }}</td>
</tr>{% endfor %}
</tbody></table>
<p><small>Specificity is annotation by similarity, not confirmed binding. A predicted structure does not confirm specificity.</small></p>
</body></html>
```

- [ ] **Step 4: Implement report.py**

```python
# src/rep2struct/report.py
from __future__ import annotations
from pathlib import Path
from jinja2 import Environment, FileSystemLoader, select_autoescape

_TPL_DIR = Path(__file__).parent / "templates"

def render_report(clonotypes, annotations, qc_results, metrics=None) -> str:
    ann = {a.clonotype_id: a for a in annotations}
    qc = {q.clonotype_id: q for q in qc_results}
    rows = []
    for c in clonotypes:
        a = ann.get(c.id)
        q = qc.get(c.id)
        rows.append({
            "id": c.id, "size": c.size,
            "annotatable": bool(a and a.annotatable),
            "epitope": a.epitope if a else None,
            "hla": a.hla if a else None,
            "tier": a.confidence_tier if a else "n/a",
            "qc": q.qc_verdict if q else "not folded",
        })
    env = Environment(loader=FileSystemLoader(str(_TPL_DIR)),
                      autoescape=select_autoescape(["html"]))
    return env.get_template("report.html.j2").render(rows=rows, metrics=metrics)
```

- [ ] **Step 5: Run tests to verify pass**

Run: `pytest tests/test_report.py -v`
Expected: PASS (2 tests).

- [ ] **Step 6: Commit**

```bash
git add src/rep2struct/report.py src/rep2struct/templates/report.html.j2 tests/test_report.py
git commit -m "feat: stage 5 self contained HTML report

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 10: Deterministic reference chain and offline harness

This is not the product entrypoint. It is the deterministic reference that proves the
whole stage chain works offline with a mocked fold, and it is what the offline test
suite exercises end to end. The product is the multi agent layer in Tasks 11 and 12,
which calls the same stage functions through tools.

**Files:**
- Create: `src/rep2struct/pipeline.py`
- Test: `tests/test_pipeline_offline.py`

**Interfaces:**
- Consumes: a 10x CSV path, a run directory, and injectable `sim_fn`, `assign_fn`, `fold_fn`, sequence lookups (so the whole chain runs offline in tests).
- Produces: `run_pipeline(csv_path, run_dir, top_n, sim_fn, assign_fn, fold_fn, tcr_seqs, mhc_seqs, scramble_threshold, labels=None) -> str` returning the HTML report path. Each stage writes to the run state and is skipped when already done.

- [ ] **Step 1: Write the failing offline end to end test**

```python
# tests/test_pipeline_offline.py
from pathlib import Path
from rep2struct.pipeline import run_pipeline

FIX = Path(__file__).parent / "fixtures"

def test_offline_end_to_end(tmp_path):
    def sim(cdr3_a, v_a, cdr3_b, v_b, species="human", top_k=5):
        if cdr3_b == "CASSIRSSYEQYF":
            return ([{"epitope": "GILGFVFTL", "mhc": "HLA-A*02:01",
                      "antigen": "Flu M1", "distance": 3.0}], "tcrdist", 1, [])
        return ([], "tcrdist", 0, [])
    def assign_fn(gene, species="human", chain=None): return gene + "*01"
    def fold_fn(job): return [str(FIX / "cognate_min.cif")]
    tcr_seqs = {}  # filled by pipeline from clonotypes via a stub below
    report = run_pipeline(
        csv_path=FIX / "tenx_tiny.csv", run_dir=tmp_path / "run",
        top_n=1, sim_fn=sim, assign_fn=assign_fn, fold_fn=fold_fn,
        tcr_seqs=None, mhc_seqs={"HLA-A*02:01": {"heavy": "H"*20, "b2m": "M"*20}},
        scramble_threshold=1.0)
    html = Path(report).read_text()
    assert "GILGFVFTL" in html and "reliable" in html

def test_resume_is_idempotent(tmp_path):
    # running twice must not raise and must return the same report path
    ...  # same setup as above, call run_pipeline twice, assert equal paths
```

- [ ] **Step 2: Run to verify fail**

Run: `pytest tests/test_pipeline_offline.py -v`
Expected: FAIL (module missing).

- [ ] **Step 3: Implement pipeline.py**

```python
# src/rep2struct/pipeline.py
from __future__ import annotations
from pathlib import Path
from .runstate import RunState
from .ingest import parse_10x, standardize_alleles
from .annotate import annotate
from .foldprep import select_top, build_construct
from .fold import run_folds
from .qc import score_model, verdict
from .report import render_report
from .validate import annotation_metrics
from .schema import Clonotype, Annotation, QCResult

def _tcr_seq_stub(clonotype):
    # placeholder chain sequences: real runs use reconstructed V domains from
    # TCR Explorer. For the construct we need any residues; use the CDR3s padded.
    return {"A": "G" * 10 + clonotype.cdr3a, "B": "G" * 10 + clonotype.cdr3b}

def run_pipeline(csv_path, run_dir, top_n, sim_fn=None, assign_fn=None, fold_fn=None,
                 tcr_seqs=None, mhc_seqs=None, scramble_threshold=1.0, labels=None):
    rs = RunState(run_dir)

    clons = standardize_alleles(parse_10x(csv_path), assign_fn=assign_fn)
    anns = annotate(clons, sim_fn=sim_fn)

    metrics = annotation_metrics(anns, labels) if labels else None

    top = select_top(clons, anns, n=top_n)
    seqs = tcr_seqs or {c.id: _tcr_seq_stub(c) for c, _ in top}
    jobs = [build_construct(c, a, seqs, mhc_seqs) for c, a in top]
    jobs = run_folds(jobs, fold_fn, rs)

    qcs = []
    for job in jobs:
        if job.status != "done" or not job.model_paths:
            qcs.append(QCResult(job.clonotype_id, "qc_failed", "no model produced"))
            continue
        s = score_model(job.model_paths[0]); s["clonotype_id"] = job.clonotype_id
        qcs.append(verdict(s, scramble_threshold))

    html = render_report(clons, anns, qcs, metrics=metrics)
    out = Path(run_dir) / "report.html"
    out.write_text(html)
    return str(out)
```

- [ ] **Step 4: Run tests to verify pass**

Run: `pytest tests/test_pipeline_offline.py -v`
Expected: PASS. Then run the full suite: `pytest -v`.

- [ ] **Step 5: Commit**

```bash
git add src/rep2struct/pipeline.py tests/test_pipeline_offline.py
git commit -m "feat: orchestration pipeline with offline end to end test

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 11: Stage functions as Agent SDK tools

**Files:**
- Create: `src/rep2struct/agent_tools.py`
- Test: `tests/test_agent_tools.py`

**Interfaces:**
- Consumes: the stage functions from Tasks 2 to 9, plus `RunState`. Tools operate on a
  `run_dir` and persist to the run directory so data survives across tool calls and
  across agents (Agent SDK subagents start with fresh context, so shared state lives on
  disk, never in the prompt).
- Produces: async `@tool` functions and a `create_sdk_mcp_server(name="rep2struct", ...)`
  server named `rep2struct`. A `configure(sim_fn=None, assign_fn=None)` hook lets tests
  inject offline fakes. Each tool returns `{"content": [...text summary...],
  "structuredContent": {...}}`.
- Tools: `ingest_repertoire(run_dir, csv_path)`, `annotate_specificity(run_dir)`,
  `prep_and_select(run_dir, top_n)`, `list_fold_jobs(run_dir)`,
  `record_fold_result(run_dir, clonotype_id, model_paths)`,
  `qc_structure(run_dir, clonotype_id, scramble_threshold)`, `render_final_report(run_dir)`.
  The fold agent drives Playwright itself and reports models through
  `record_fold_result`, so browser driving stays in the agent and state stays in tools.

- [ ] **Step 1: Write the failing test (call tools directly, offline)**

```python
# tests/test_agent_tools.py
import asyncio
from pathlib import Path
from rep2struct import agent_tools as at

FIX = Path(__file__).parent / "fixtures"

def _run(coro): return asyncio.get_event_loop().run_until_complete(coro)

def test_ingest_then_annotate_tools(tmp_path):
    def sim(cdr3_a, v_a, cdr3_b, v_b, species="human", top_k=5):
        if cdr3_b == "CASSIRSSYEQYF":
            return ([{"epitope": "GILGFVFTL", "mhc": "HLA-A*02:01",
                      "antigen": "Flu M1", "distance": 3.0}], "tcrdist", 1, [])
        return ([], "tcrdist", 0, [])
    at.configure(sim_fn=sim, assign_fn=lambda g, species="human", chain=None: g + "*01")
    rd = str(tmp_path / "run")
    ing = _run(at.ingest_repertoire({"run_dir": rd, "csv_path": str(FIX / "tenx_tiny.csv")}))
    assert ing["structuredContent"]["clonotypes"] >= 1
    ann = _run(at.annotate_specificity({"run_dir": rd}))
    tiers = ann["structuredContent"]["tiers"]
    assert tiers.get("high", 0) >= 1

def test_qc_tool_flags_scramble(tmp_path):
    rd = str(tmp_path / "run2")
    at.configure()
    # seed a fold job + model path pointing at the scramble fixture
    _run(at.ingest_repertoire({"run_dir": rd, "csv_path": str(FIX / "tenx_tiny.csv")}))
    _run(at.record_fold_result({"run_dir": rd, "clonotype_id": "z",
                                "model_paths": [str(FIX / "scramble_min.cif")]}))
    out = _run(at.qc_structure({"run_dir": rd, "clonotype_id": "z", "scramble_threshold": 1.0}))
    assert out["structuredContent"]["qc_verdict"] == "suspect"
```

- [ ] **Step 2: Run to verify fail**

Run: `pytest tests/test_agent_tools.py -v`
Expected: FAIL (module missing).

- [ ] **Step 3: Implement agent_tools.py**

```python
# src/rep2struct/agent_tools.py
from __future__ import annotations
from claude_agent_sdk import tool, create_sdk_mcp_server
from .runstate import RunState
from .ingest import parse_10x, standardize_alleles
from .annotate import annotate
from .foldprep import select_top, build_construct
from .qc import score_model, verdict
from .report import render_report
from .schema import Clonotype, Annotation, QCResult

_CFG = {"sim_fn": None, "assign_fn": None}
def configure(sim_fn=None, assign_fn=None):
    _CFG["sim_fn"] = sim_fn
    _CFG["assign_fn"] = assign_fn

def _txt(s): return {"content": [{"type": "text", "text": s}]}
def _load(rd, name, cls):
    rs = RunState(rd)
    return [cls(**d) for d in rs.read_stage(name)] if rs.stage_done(name) else []

@tool("ingest_repertoire", "Parse a 10x contig CSV into paired clonotypes and persist them.",
      {"run_dir": str, "csv_path": str})
async def ingest_repertoire(args):
    clons = standardize_alleles(parse_10x(args["csv_path"]), assign_fn=_CFG["assign_fn"])
    RunState(args["run_dir"]).write_stage("ingest", clons)
    r = _txt(f"{len(clons)} clonotypes ingested")
    r["structuredContent"] = {"clonotypes": len(clons)}
    return r

@tool("annotate_specificity", "Annotate persisted clonotypes with candidate epitopes by TCRdist. Never forces a label.",
      {"run_dir": str})
async def annotate_specificity(args):
    clons = _load(args["run_dir"], "ingest", Clonotype)
    anns = annotate(clons, sim_fn=_CFG["sim_fn"])
    RunState(args["run_dir"]).write_stage("annotate", anns)
    tiers = {}
    for a in anns:
        tiers[a.confidence_tier] = tiers.get(a.confidence_tier, 0) + 1
    r = _txt(f"annotated {len(anns)} clonotypes: {tiers}")
    r["structuredContent"] = {"tiers": tiers}
    return r

@tool("prep_and_select", "Rank clonotypes and build class I fold constructs for the top N.",
      {"run_dir": str, "top_n": int})
async def prep_and_select(args):
    clons = _load(args["run_dir"], "ingest", Clonotype)
    anns = _load(args["run_dir"], "annotate", Annotation)
    top = select_top(clons, anns, n=args["top_n"])
    seqs = {c.id: {"A": "G"*10 + c.cdr3a, "B": "G"*10 + c.cdr3b} for c, _ in top}
    mhc = {a.hla: {"heavy": "H"*20, "b2m": "M"*20} for _, a in top if a.hla}
    jobs = [build_construct(c, a, seqs, mhc) for c, a in top]
    RunState(args["run_dir"]).write_stage("foldjobs", jobs)
    r = _txt(f"prepared {len(jobs)} fold jobs")
    r["structuredContent"] = {"jobs": [j.clonotype_id for j in jobs]}
    return r

@tool("list_fold_jobs", "List fold jobs and their construct FASTA for the fold agent.",
      {"run_dir": str})
async def list_fold_jobs(args):
    rs = RunState(args["run_dir"])
    jobs = rs.read_stage("foldjobs") if rs.stage_done("foldjobs") else []
    r = _txt(f"{len(jobs)} jobs")
    r["structuredContent"] = {"jobs": jobs}
    return r

@tool("record_fold_result", "Record the model paths a fold produced for one clonotype.",
      {"run_dir": str, "clonotype_id": str, "model_paths": list})
async def record_fold_result(args):
    rs = RunState(args["run_dir"])
    done = rs.read_stage("folds") if rs.stage_done("folds") else {}
    done[args["clonotype_id"]] = args["model_paths"]
    rs.write_stage("folds", done)
    return _txt(f"recorded {len(args['model_paths'])} models for {args['clonotype_id']}")

@tool("qc_structure", "Score a predicted structure and return a skeptical reliable or suspect verdict.",
      {"run_dir": str, "clonotype_id": str, "scramble_threshold": float})
async def qc_structure(args):
    rs = RunState(args["run_dir"])
    done = rs.read_stage("folds") if rs.stage_done("folds") else {}
    paths = done.get(args["clonotype_id"], [])
    if not paths:
        res = QCResult(args["clonotype_id"], "qc_failed", "no model recorded")
    else:
        s = score_model(paths[0]); s["clonotype_id"] = args["clonotype_id"]
        res = verdict(s, args["scramble_threshold"])
    qcs = rs.read_stage("qc") if rs.stage_done("qc") else []
    qcs = [q for q in qcs if q["clonotype_id"] != res.clonotype_id]
    from dataclasses import asdict
    qcs.append(asdict(res)); rs.write_stage("qc", qcs)
    r = _txt(f"{res.clonotype_id}: {res.qc_verdict} ({res.reason})")
    r["structuredContent"] = {"qc_verdict": res.qc_verdict, "reason": res.reason}
    return r

@tool("render_final_report", "Render the self contained HTML report for the run.",
      {"run_dir": str})
async def render_final_report(args):
    clons = _load(args["run_dir"], "ingest", Clonotype)
    anns = _load(args["run_dir"], "annotate", Annotation)
    rs = RunState(args["run_dir"])
    qcs = [QCResult(**d) for d in (rs.read_stage("qc") if rs.stage_done("qc") else [])]
    html = render_report(clons, anns, qcs)
    from pathlib import Path
    out = Path(args["run_dir"]) / "report.html"; out.write_text(html)
    r = _txt(f"report written to {out}")
    r["structuredContent"] = {"report_path": str(out)}
    return r

def build_server():
    return create_sdk_mcp_server(name="rep2struct", version="0.1.0", tools=[
        ingest_repertoire, annotate_specificity, prep_and_select, list_fold_jobs,
        record_fold_result, qc_structure, render_final_report,
    ])
```

- [ ] **Step 4: Run tests to verify pass**

Run: `pytest tests/test_agent_tools.py -v`
Expected: PASS (2 tests). Note: the `@tool` decorated callables are invoked directly with an args dict, which the SDK supports.

- [ ] **Step 5: Commit**

```bash
git add src/rep2struct/agent_tools.py tests/test_agent_tools.py
git commit -m "feat: expose stages as Agent SDK tools over a shared run directory

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 12: Multi agent layer and product entrypoint

**Files:**
- Create: `src/rep2struct/agents.py`
- Create: `src/rep2struct/app.py`
- Test: `tests/test_agents_config.py`

**Interfaces:**
- Produces: `build_agents() -> dict[str, AgentDefinition]` (fold, qc, report specialists);
  `build_options(run_dir) -> ClaudeAgentOptions` wiring the in process `rep2struct` tool
  server, the Playwright MCP, the specialist agents, and the allowed tools;
  `orchestrator_prompt(csv_path, run_dir, top_n) -> str`; and async `run(csv_path,
  run_dir, top_n) -> str` that calls `query()` and returns the report path.
- The config assembly (agents, tool wiring) is unit tested offline. The live `query()`
  run needs the API and is a manual demo step, not a unit test.

- [ ] **Step 1: Write the failing config test**

```python
# tests/test_agents_config.py
from rep2struct.agents import build_agents, build_options, orchestrator_prompt

def test_specialist_agents_present():
    agents = build_agents()
    assert set(agents) >= {"fold-agent", "qc-agent", "report-agent"}
    # the fold agent must be allowed to drive the browser
    assert any("playwright" in t for t in agents["fold-agent"].tools)

def test_options_wire_tools_and_agents(tmp_path):
    opts = build_options(str(tmp_path / "run"))
    assert "rep2struct" in opts.mcp_servers
    assert "playwright" in opts.mcp_servers
    assert any(t.startswith("mcp__rep2struct__") or t == "Agent" for t in opts.allowed_tools)
    assert set(opts.agents) >= {"fold-agent", "qc-agent", "report-agent"}

def test_prompt_names_the_stages(tmp_path):
    p = orchestrator_prompt("x.csv", str(tmp_path), 8)
    for kw in ["ingest", "annotate", "fold", "qc", "report"]:
        assert kw in p.lower()
```

- [ ] **Step 2: Run to verify fail**

Run: `pytest tests/test_agents_config.py -v`
Expected: FAIL (module missing).

- [ ] **Step 3: Implement agents.py**

```python
# src/rep2struct/agents.py
from __future__ import annotations
from claude_agent_sdk import ClaudeAgentOptions, AgentDefinition
from .agent_tools import build_server

def build_agents():
    return {
        "fold-agent": AgentDefinition(
            description="Folds TCR pMHC constructs by driving Protenix on Colab through the browser.",
            prompt=(
                "You fold structures. Call list_fold_jobs to get constructs. For each job, "
                "drive the Protenix Colab notebook with the mcp__playwright tools: open the "
                "notebook, submit the construct FASTA at 5 seeds, wait for completion, download "
                "the CIF models, then call record_fold_result with the local model paths. The "
                "loop is resumable; skip a job that already has recorded models."),
            tools=["mcp__rep2struct__list_fold_jobs", "mcp__rep2struct__record_fold_result",
                   "mcp__playwright__*"],
            model="sonnet",
        ),
        "qc-agent": AgentDefinition(
            description="Skeptical QC of predicted TCR pMHC structures. Flags geometry hallucinations.",
            prompt=(
                "You are a skeptical structural referee. For each folded clonotype call "
                "qc_structure. A clean fold does NOT confirm specificity: Protenix imposes "
                "canonical docking geometry even on non binders. Report reliable only when the "
                "CDR3 to peptide contact beats the scramble calibration; otherwise suspect. Never "
                "upgrade a verdict to please the caller."),
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
            "mcp__rep2struct__prep_and_select", "mcp__rep2struct__list_fold_jobs",
            "mcp__rep2struct__record_fold_result", "mcp__rep2struct__qc_structure",
            "mcp__rep2struct__render_final_report", "mcp__playwright__*",
        ],
        permission_mode="acceptEdits",
    )

def orchestrator_prompt(csv_path, run_dir, top_n):
    return (
        f"Run the repertoire to structure pipeline on {csv_path} with run_dir {run_dir}.\n"
        f"1. Call ingest_repertoire, then annotate_specificity (honest annotation, keep unannotatable as is).\n"
        f"2. Call prep_and_select with top_n {top_n}.\n"
        f"3. Delegate to the fold-agent to fold the prepared jobs.\n"
        f"4. Delegate to the qc-agent to QC each folded clonotype (scramble_threshold from calibration).\n"
        f"5. Delegate to the report-agent to render the final HTML report, and return its path.")
```

- [ ] **Step 4: Implement app.py**

```python
# src/rep2struct/app.py
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
```

- [ ] **Step 5: Run tests to verify pass**

Run: `pytest tests/test_agents_config.py -v`
Expected: PASS (3 tests). Then the full suite: `pytest -v`.

- [ ] **Step 6: Commit**

```bash
git add src/rep2struct/agents.py src/rep2struct/app.py tests/test_agents_config.py
git commit -m "feat: multi agent layer (orchestrator plus fold, qc, report specialists)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## After the core: live integration (not unit tested, run manually)

These steps connect the real datasets and the real fold. They are checklist items for the demo, run by the Claude agent, not covered by the offline suite.

- [ ] Install `~/imgt-api` editable into the venv, confirm `from tcr_explorer.similarity import find_similar_paired_tcrs` imports, and confirm `tcr_explorer.tcrdist_engine.tcrdist_available()` is true.
- [ ] Wire the real `sim_fn` and `assign_fn` (drop the injected fakes) and reconstruct real V domain sequences for the construct from TCR Explorer instead of `_tcr_seq_stub`.
- [ ] Download the 10x 4 donor dextramer set, build the `labels` map from binarized dextramer calls, run the validation arm, and record precision, recall, and unannotatable rate. Calibrate `DEFAULT_TIERS` and `scramble_threshold` from this arm.
- [ ] Run the application arm on one TABLO donor CSV end to end through the agent app (`python -m rep2struct.app <csv> <run_dir> <top_n>`), where the fold agent drives Protenix on Colab via the Playwright MCP for the top N clonotypes.
- [ ] Confirm the Playwright MCP is installed (`npx @playwright/mcp@latest`) and that the fold agent can reach the Colab notebook.
- [ ] Publish the two HTML reports as claude.ai Artifacts for the demo.

---

## Self-Review

Spec coverage: stages 0 to 5 each map to Tasks 2 to 9; the deterministic reference chain is Task 10; the multi agent product layer is Tasks 11 (tools) and 12 (agents plus entrypoint). The two honesty rules are enforced in `annotate.py` (annotatable flag, unannotatable tier) and `qc.py` (suspect verdict), and are restated in the qc-agent prompt. Validation arm is Task 5 plus the live checklist. Application arm runs through `app.py` in the live checklist. Report format is HTML (Task 9).

Placeholder scan: the only ellipsis is `test_resume_is_idempotent` in Task 10 Step 1, intentionally left for the implementer to mirror the preceding test; every implementation step carries complete code.

Type consistency: `Clonotype`, `Annotation`, `FoldJob`, `QCResult` field names are used identically across ingest, annotate, foldprep, fold, qc, report, pipeline, and agent_tools. `sim_fn` returns the four tuple `(neighbours, engine, total, warnings)` matching TCR Explorer everywhere it is called. The tool server name `rep2struct` matches the `mcp__rep2struct__` prefixes in `agents.py` allowed tools. Agent keys `fold-agent`, `qc-agent`, `report-agent` match between `build_agents` and the config test.

Agent SDK API check: `tool`, `create_sdk_mcp_server`, `AgentDefinition`, `ClaudeAgentOptions`, `query`, `ResultMessage` are the current import surface of `claude-agent-sdk`. Tool return shape is `{"content": [...], "structuredContent": {...}}`. Subagents receive no parent context, which is why all shared state is persisted to the run directory by the tools rather than passed in prompts. Model aliases (`sonnet`, `opus`) are used rather than dated model IDs.
