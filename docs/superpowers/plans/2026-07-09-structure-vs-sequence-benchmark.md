# Structure-vs-Sequence Retrieval Benchmark Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Measure whether the scramble-calibrated CDR3-peptide contact recovers the correct epitope from a candidate panel for sequence-novel TCRs, where sequence annotation fails.

**Architecture:** A new `rep2struct.benchmark` module of pure functions (panel/decoy selection, construct building, retrieval scoring, baselines, statistics) built on the released R2S modules (`annotate`, `foldprep`, `qc`, `seqs`, `tools.protenix_inputs`). A driver `scripts/run_benchmark_arm.py` selects TCRs from the 10x dextramer data, emits fold constructs and the Colab notebook (the user runs the folds), then scores the returned structures. The paper's Methods cite the released package.

**Tech Stack:** Python 3.11, numpy, pandas, Biopython (MMCIFParser), pytest. Existing R2S package installed editable (`pip install -e .`).

## Global Constraints

- Reuse released R2S modules; new code is additive (`src/rep2struct/benchmark.py`, `scripts/run_benchmark_arm.py`). Copied verbatim from spec.
- The contact metric is reused **unchanged** from `qc.ensemble_contact` / `qc.score_model` — do not redefine it. Its definition (chain B = TCRbeta, chain E = peptide, heavy-atom < 4.5 A, median across samples) is frozen before looking at retrieval results.
- Class I only. Protenix only. No new folding model. No class II.
- Novelty (leakage) guard: a TCR is **novel** iff its nearest reference TCRdist is `None` or `> 1.0`.
- Folds run through the existing user-driven Colab notebook: Claude builds/uploads the notebook, the user runs cells and pastes outputs. No Playwright.
- Keep the `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>` trailer on commits (R2S convention).
- Seed target: one HLA, ~10-15 TCRs, cognate + 3 decoys, 3 samples (~120-240 folds), by 2026-07-13.
- All statistics report bootstrap CIs and a permutation p-value; no claim beyond the pilot N.

---

### Task 1: Panel and novel-TCR pre-check

Gates the whole benchmark: confirm one HLA has enough de-leaked novel labeled TCRs. Reuses the validation-arm loader and the real annotation.

**Files:**
- Create: `src/rep2struct/benchmark.py`
- Create: `tests/test_benchmark.py`
- Test: `tests/test_benchmark.py`

**Interfaces:**
- Consumes: `scripts/run_validation_arm.py::labeled_clonotypes(contig_csv, bm_csv)` which returns **three** values `(clons: list[Clonotype], labels: dict[id,str], hlas: dict[id,str])`; the driver composes `truth[cid] = (labels[cid], hlas[cid])`. `benchmark` functions take that composed `truth` dict. Annotations come from the validation arm's cached path `nearest_cache(clons)` + `annotations_from_cache(clons, cache)` (avoids thousands of live similarity queries), not a direct `annotate()` over the full set.
- Produces:
  - `panel_epitopes(truth: dict[str, tuple[str,str]]) -> list[tuple[str,str]]` — sorted unique `(epitope, hla)` pairs.
  - `is_novel(tcrdist: float | None, leak_thr: float = 1.0) -> bool`.
  - `per_hla_novel_counts(clonotypes, truth, annotations) -> dict[str, dict]` — `{hla: {"n_total": int, "n_novel": int, "epitopes": {epitope: {"n": int, "n_novel": int}}}}`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_benchmark.py
from rep2struct import benchmark as bm
from rep2struct.schema import Clonotype, Annotation

def test_is_novel():
    assert bm.is_novel(None) is True
    assert bm.is_novel(0.0) is False
    assert bm.is_novel(1.0) is False
    assert bm.is_novel(1.5) is True

def test_panel_epitopes_sorted_unique():
    truth = {"c1": ("GILGFVFTL", "HLA-A*02:01"),
             "c2": ("GILGFVFTL", "HLA-A*02:01"),
             "c3": ("NLVPMVATV", "HLA-A*02:01")}
    assert bm.panel_epitopes(truth) == [
        ("GILGFVFTL", "HLA-A*02:01"), ("NLVPMVATV", "HLA-A*02:01")]

def test_per_hla_novel_counts():
    clonos = [Clonotype("c1","TRAV1","CAA","TRBV1","CBB",5),
              Clonotype("c2","TRAV1","CAC","TRBV1","CBD",3)]
    truth = {"c1": ("GILGFVFTL","HLA-A*02:01"),
             "c2": ("GILGFVFTL","HLA-A*02:01")}
    anns = [Annotation("c1", False, "unannotatable", tcrdist=None),
            Annotation("c2", True, "high", tcrdist=0.0)]
    out = bm.per_hla_novel_counts(clonos, truth, anns)
    assert out["HLA-A*02:01"]["n_total"] == 2
    assert out["HLA-A*02:01"]["n_novel"] == 1
    assert out["HLA-A*02:01"]["epitopes"]["GILGFVFTL"]["n_novel"] == 1
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./.venv/bin/python -m pytest tests/test_benchmark.py -v`
Expected: FAIL with `ModuleNotFoundError` or `AttributeError: module ... has no attribute 'is_novel'`.

- [ ] **Step 3: Write minimal implementation**

```python
# src/rep2struct/benchmark.py
from __future__ import annotations
from collections import defaultdict

def is_novel(tcrdist, leak_thr: float = 1.0) -> bool:
    return tcrdist is None or tcrdist > leak_thr

def panel_epitopes(truth):
    return sorted({(ep, hla) for (ep, hla) in truth.values()})

def per_hla_novel_counts(clonotypes, truth, annotations):
    dist = {a.clonotype_id: getattr(a, "tcrdist", None) for a in annotations}
    out = {}
    for c in clonotypes:
        if c.id not in truth:
            continue
        ep, hla = truth[c.id]
        novel = is_novel(dist.get(c.id))
        h = out.setdefault(hla, {"n_total": 0, "n_novel": 0, "epitopes": {}})
        h["n_total"] += 1
        h["n_novel"] += int(novel)
        e = h["epitopes"].setdefault(ep, {"n": 0, "n_novel": 0})
        e["n"] += 1
        e["n_novel"] += int(novel)
    return out
```

- [ ] **Step 4: Run test to verify it passes**

Run: `./.venv/bin/python -m pytest tests/test_benchmark.py -v`
Expected: PASS (3 tests).

- [ ] **Step 5: Commit**

```bash
git add src/rep2struct/benchmark.py tests/test_benchmark.py
git commit -m "feat(benchmark): panel epitopes + per-HLA novel-TCR pre-check

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: Decoy selection and panel construct building

Build the cognate + k decoy fold constructs for one TCR, reusing `foldprep.build_construct` with synthetic decoy annotations.

**Files:**
- Modify: `src/rep2struct/benchmark.py`
- Test: `tests/test_benchmark.py`

**Interfaces:**
- Consumes: `rep2struct.foldprep.build_construct(clonotype, annotation, tcr_seqs, mhc_seqs) -> FoldJob`; `rep2struct.schema.Annotation`, `rep2struct.schema.Clonotype`.
- Produces:
  - `decoys_for(cognate: str, hla: str, panel: list[tuple[str,str]], k: int) -> list[str]` — up to `k` decoy peptides, same HLA first (sorted, deterministic), then filled from other HLAs if short.
  - `build_panel_constructs(clonotype, cognate: str, hla: str, decoys: list[str], tcr_seqs, mhc_seqs) -> dict[str, "FoldJob"]` — keyed by epitope (cognate + each decoy).

- [ ] **Step 1: Write the failing test**

```python
# append to tests/test_benchmark.py
def test_decoys_same_hla_first():
    panel = [("GILGFVFTL","HLA-A*02:01"), ("NLVPMVATV","HLA-A*02:01"),
             ("GLCTLVAML","HLA-A*02:01"), ("KLGGALQAK","HLA-A*03:01")]
    d = bm.decoys_for("GILGFVFTL", "HLA-A*02:01", panel, k=2)
    assert d == ["GLCTLVAML", "NLVPMVATV"]      # sorted same-HLA, cognate excluded
    assert "GILGFVFTL" not in d

def test_decoys_fill_from_other_hla_when_short():
    panel = [("GILGFVFTL","HLA-A*02:01"), ("KLGGALQAK","HLA-A*03:01")]
    d = bm.decoys_for("GILGFVFTL", "HLA-A*02:01", panel, k=2)
    assert d == ["KLGGALQAK"]                    # only one other epitope exists

def test_build_panel_constructs_keys_and_peptides():
    clono = Clonotype("c1","TRAV1","CAA","TRBV1","CBB",5)
    tcr_seqs = {"c1": {"A": "AAAA", "B": "BBBB"}}
    mhc_seqs = {"HLA-A*02:01": {"heavy": "HHHH", "b2m": "MMMM"}}
    jobs = bm.build_panel_constructs(clono, "GILGFVFTL", "HLA-A*02:01",
                                     ["NLVPMVATV"], tcr_seqs, mhc_seqs)
    assert set(jobs) == {"GILGFVFTL", "NLVPMVATV"}
    assert ">E\nGILGFVFTL" in jobs["GILGFVFTL"].construct_fasta
    assert ">E\nNLVPMVATV" in jobs["NLVPMVATV"].construct_fasta
    assert jobs["GILGFVFTL"].clonotype_id == "c1"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./.venv/bin/python -m pytest tests/test_benchmark.py -k "decoys or build_panel" -v`
Expected: FAIL with `AttributeError: ... has no attribute 'decoys_for'`.

- [ ] **Step 3: Write minimal implementation**

```python
# append to src/rep2struct/benchmark.py
from .schema import Annotation
from .foldprep import build_construct

def decoys_for(cognate, hla, panel, k):
    same = sorted(ep for (ep, h) in panel if h == hla and ep != cognate)
    other = sorted(ep for (ep, h) in panel if h != hla and ep != cognate)
    return (same + other)[:k]

def build_panel_constructs(clonotype, cognate, hla, decoys, tcr_seqs, mhc_seqs):
    jobs = {}
    for ep in [cognate, *decoys]:
        ann = Annotation(clonotype_id=clonotype.id, annotatable=True,
                         confidence_tier="benchmark", epitope=ep, hla=hla)
        jobs[ep] = build_construct(clonotype, ann, tcr_seqs, mhc_seqs)
    return jobs
```

- [ ] **Step 4: Run test to verify it passes**

Run: `./.venv/bin/python -m pytest tests/test_benchmark.py -v`
Expected: PASS (all tests).

- [ ] **Step 5: Commit**

```bash
git add src/rep2struct/benchmark.py tests/test_benchmark.py
git commit -m "feat(benchmark): decoy selection + panel construct builder

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: Retrieval scoring on the contact metric

Rank epitopes by the reused contact metric and compute Top-1 and pairwise AUROC. Tested on the existing cif fixtures (no real folds needed).

**Files:**
- Modify: `src/rep2struct/benchmark.py`
- Test: `tests/test_benchmark.py`
- Uses fixtures: `tests/fixtures/cognate_min.cif`, `tests/fixtures/scramble_min.cif`

**Interfaces:**
- Consumes: `rep2struct.qc.ensemble_contact(paths) -> (median|None, n_models, n_valid)`.
- Produces:
  - `contact_by_epitope(paths_by_epitope: dict[str, list[str]]) -> dict[str, float | None]` — median contact per epitope via `qc.ensemble_contact`.
  - `retrieval_result(contacts: dict[str, float | None], cognate: str) -> dict` with keys `ranked` (list[str], high contact first, `None`->-1), `top1` (bool), `cognate_contact` (float|None).
  - `auroc(pairs: list[tuple[float, list[float]]]) -> float | None` — mean over all (cognate, decoy) pairs of `1[cognate>decoy] + 0.5*1[==]`; `None` if no pairs.

- [ ] **Step 1: Write the failing test**

```python
# append to tests/test_benchmark.py
from pathlib import Path
FIX = Path(__file__).parent / "fixtures"

def test_contact_and_retrieval_with_fixtures():
    paths = {"COGNATE": [str(FIX/"cognate_min.cif")],
             "DECOY":   [str(FIX/"scramble_min.cif")]}
    contacts = bm.contact_by_epitope(paths)
    assert contacts["COGNATE"] >= contacts["DECOY"]     # cognate fixture has more contact
    res = bm.retrieval_result(contacts, "COGNATE")
    assert res["ranked"][0] == "COGNATE"
    assert res["top1"] is True

def test_auroc_pairs():
    assert bm.auroc([(10.0, [1.0, 2.0])]) == 1.0        # cognate beats both
    assert bm.auroc([(1.0, [10.0])]) == 0.0
    assert bm.auroc([(5.0, [5.0])]) == 0.5              # tie
    assert bm.auroc([]) is None
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./.venv/bin/python -m pytest tests/test_benchmark.py -k "retrieval or auroc" -v`
Expected: FAIL with `AttributeError: ... 'contact_by_epitope'`.

- [ ] **Step 3: Write minimal implementation**

```python
# append to src/rep2struct/benchmark.py
from .qc import ensemble_contact

def contact_by_epitope(paths_by_epitope):
    return {ep: ensemble_contact(paths)[0] for ep, paths in paths_by_epitope.items()}

def retrieval_result(contacts, cognate):
    ranked = sorted(contacts, key=lambda e: (-1.0 if contacts[e] is None else contacts[e]),
                    reverse=True)
    return {"ranked": ranked, "top1": ranked[0] == cognate,
            "cognate_contact": contacts.get(cognate)}

def auroc(pairs):
    num = den = 0.0
    for cog, decoys in pairs:
        for d in decoys:
            den += 1
            num += 1.0 if cog > d else (0.5 if cog == d else 0.0)
    return None if den == 0 else num / den
```

- [ ] **Step 4: Run test to verify it passes**

Run: `./.venv/bin/python -m pytest tests/test_benchmark.py -v`
Expected: PASS. If the fixtures' relative contact is unexpectedly equal, inspect with
`./.venv/bin/python -c "from rep2struct.qc import score_model; print(score_model('tests/fixtures/cognate_min.cif'), score_model('tests/fixtures/scramble_min.cif'))"` and adjust the assertion to the observed direction (do NOT change the metric).

- [ ] **Step 5: Commit**

```bash
git add src/rep2struct/benchmark.py tests/test_benchmark.py
git commit -m "feat(benchmark): contact retrieval scoring (top-1, AUROC)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: Baselines — sequence (tcrdist) and peptide pLDDT

Two honest competitors on the same folds. B1 (sequence) is chance-on-novel by construction; B2 (peptide pLDDT) is the 2025 prior-art confidence signal (peptide-chain pLDDT proxy for the seed; true CDR3 residues deferred to preprint).

**Files:**
- Modify: `src/rep2struct/benchmark.py`
- Test: `tests/test_benchmark.py`

**Interfaces:**
- Consumes: `qc.mean_confidence(bfactors) -> float | None`; Biopython `MMCIFParser` (already a dependency, used in `qc._heavy_by_chain`).
- Produces:
  - `model_plddt(cif_path: str) -> float | None` — mean B-factor over chain E (peptide) atoms; `None` if chain E absent.
  - `plddt_by_epitope(paths_by_epitope) -> dict[str, float | None]` — median peptide pLDDT per epitope across samples.
  - `sequence_baseline_top1(annotation_epitope: str | None, cognate: str) -> bool` — whether the tcrdist annotation names the cognate (novel TCRs annotate to `None` -> False).

- [ ] **Step 1: Write the failing test**

```python
# append to tests/test_benchmark.py
def test_model_plddt_reads_peptide_bfactors():
    v = bm.model_plddt(str(FIX/"cognate_min.cif"))
    assert v is None or isinstance(v, float)            # fixture may lack chain E bfactors
    assert bm.model_plddt(str(FIX/"threechain_min.cif")) is None  # no chain E

def test_sequence_baseline_top1():
    assert bm.sequence_baseline_top1("GILGFVFTL", "GILGFVFTL") is True
    assert bm.sequence_baseline_top1("NLVPMVATV", "GILGFVFTL") is False
    assert bm.sequence_baseline_top1(None, "GILGFVFTL") is False
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./.venv/bin/python -m pytest tests/test_benchmark.py -k "plddt or sequence_baseline" -v`
Expected: FAIL with `AttributeError: ... 'model_plddt'`.

- [ ] **Step 3: Write minimal implementation**

```python
# append to src/rep2struct/benchmark.py
import warnings
import numpy as np
from .qc import mean_confidence

def model_plddt(cif_path):
    from Bio.PDB import MMCIFParser
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = next(MMCIFParser(QUIET=True).get_structure("x", str(cif_path)).get_models())
    for ch in m:
        if ch.id == "E":
            bf = [a.get_bfactor() for r in ch for a in r if a.element != "H"]
            return mean_confidence(bf)
    return None

def plddt_by_epitope(paths_by_epitope):
    out = {}
    for ep, paths in paths_by_epitope.items():
        vals = [v for v in (model_plddt(p) for p in paths) if v is not None]
        out[ep] = float(np.median(vals)) if vals else None
    return out

def sequence_baseline_top1(annotation_epitope, cognate):
    return annotation_epitope == cognate
```

- [ ] **Step 4: Run test to verify it passes**

Run: `./.venv/bin/python -m pytest tests/test_benchmark.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/rep2struct/benchmark.py tests/test_benchmark.py
git commit -m "feat(benchmark): tcrdist + peptide-pLDDT baselines

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 5: Statistics — bootstrap CI and permutation test

Turn per-TCR hits into an honest headline with uncertainty. Deterministic via seeded `random.Random`.

**Files:**
- Modify: `src/rep2struct/benchmark.py`
- Test: `tests/test_benchmark.py`

**Interfaces:**
- Produces:
  - `bootstrap_ci(hits: list[bool], n_boot: int = 2000, seed: int = 0) -> tuple[float, float, float]` — `(point_mean, lo_2.5pct, hi_97.5pct)`.
  - `permutation_p(hits: list[bool], chance: float, n_perm: int = 10000, seed: int = 0) -> float` — one-sided p that mean(hits) > chance under a Binomial(len,chance) null.

- [ ] **Step 1: Write the failing test**

```python
# append to tests/test_benchmark.py
def test_bootstrap_ci_all_hits():
    pt, lo, hi = bm.bootstrap_ci([True]*20, n_boot=500, seed=1)
    assert pt == 1.0
    assert lo == 1.0 and hi == 1.0

def test_bootstrap_ci_bounds_order():
    pt, lo, hi = bm.bootstrap_ci([True, False, True, False], n_boot=500, seed=1)
    assert 0.0 <= lo <= pt <= hi <= 1.0

def test_permutation_p_strong_signal():
    p = bm.permutation_p([True]*10, chance=0.25, n_perm=5000, seed=1)
    assert p < 0.01
    p2 = bm.permutation_p([False]*10, chance=0.25, n_perm=5000, seed=1)
    assert p2 > 0.5
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./.venv/bin/python -m pytest tests/test_benchmark.py -k "bootstrap or permutation" -v`
Expected: FAIL with `AttributeError: ... 'bootstrap_ci'`.

- [ ] **Step 3: Write minimal implementation**

```python
# append to src/rep2struct/benchmark.py
import random as _random

def bootstrap_ci(hits, n_boot: int = 2000, seed: int = 0):
    n = len(hits)
    pt = sum(hits) / n if n else 0.0
    if n == 0:
        return 0.0, 0.0, 0.0
    rng = _random.Random(seed)
    means = []
    for _ in range(n_boot):
        s = sum(hits[rng.randrange(n)] for _ in range(n)) / n
        means.append(s)
    means.sort()
    lo = means[int(0.025 * n_boot)]
    hi = means[min(int(0.975 * n_boot), n_boot - 1)]
    return pt, lo, hi

def permutation_p(hits, chance, n_perm: int = 10000, seed: int = 0):
    n = len(hits)
    obs = sum(hits)
    if n == 0:
        return 1.0
    rng = _random.Random(seed)
    ge = 0
    for _ in range(n_perm):
        draw = sum(1 for _ in range(n) if rng.random() < chance)
        if draw >= obs:
            ge += 1
    return (ge + 1) / (n_perm + 1)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `./.venv/bin/python -m pytest tests/test_benchmark.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/rep2struct/benchmark.py tests/test_benchmark.py
git commit -m "feat(benchmark): bootstrap CI + permutation test

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 6: Selection + fold-emission driver (pre-fold half)

Driver that loads the dextramer data, runs the pre-check, selects the seed TCRs, builds panel constructs, emits Protenix inputs and the Colab notebook for the user to run. **Ends at the MANUAL FOLD GATE.**

**Files:**
- Create: `scripts/run_benchmark_arm.py`
- Test: `tests/test_benchmark_driver.py`

**Interfaces:**
- Consumes: `benchmark.per_hla_novel_counts`, `benchmark.panel_epitopes`, `benchmark.decoys_for`, `benchmark.build_panel_constructs`; `rep2struct.seqs.build_tcr_seqs`, `rep2struct.seqs.build_mhc_seqs`; `rep2struct.tools.protenix_inputs.build`; `scripts/build_colab_notebook.py` (imported as a module or invoked via subprocess to produce the notebook); `scripts/run_validation_arm.py::labeled_clonotypes`.
- Produces:
  - `select_seed_tcrs(clonotypes, truth, annotations, hla: str, n: int, prefer_novel: bool = True) -> list[str]` — clonotype ids in the chosen HLA, novel-first, then by size, deterministic.
  - `emit_manifest(out_dir, selected, truth, annotations, panel, tcr_seqs, mhc_seqs, k, samples) -> dict` — writes `constructs/*.fasta` + `manifest.json` mapping `clonotype_id -> {cognate, hla, decoys, epitopes: {epitope: fasta_path}, novel: bool, tcrdist: float|None}`; returns the manifest dict.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_benchmark_driver.py
import json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import run_benchmark_arm as drv
from rep2struct.schema import Clonotype, Annotation

def test_select_seed_novel_first():
    clonos = [Clonotype("a","TRAV1","CAA","TRBV1","CBB",2),
              Clonotype("b","TRAV1","CAC","TRBV1","CBD",9)]
    truth = {"a": ("GILGFVFTL","HLA-A*02:01"), "b": ("GILGFVFTL","HLA-A*02:01")}
    anns = [Annotation("a", False, "unannotatable", tcrdist=None),   # novel
            Annotation("b", True, "high", tcrdist=0.0)]              # leaked
    sel = drv.select_seed_tcrs(clonos, truth, anns, "HLA-A*02:01", n=2)
    assert sel[0] == "a"                                             # novel first

def test_emit_manifest_writes_constructs(tmp_path):
    clonos = [Clonotype("a","TRAV1","CAA","TRBV1","CBB",2)]
    truth = {"a": ("GILGFVFTL","HLA-A*02:01")}
    anns = [Annotation("a", False, "unannotatable", tcrdist=None)]
    panel = [("GILGFVFTL","HLA-A*02:01"), ("NLVPMVATV","HLA-A*02:01")]
    tcr_seqs = {"a": {"A": "AAAA", "B": "BBBB"}}
    mhc_seqs = {"HLA-A*02:01": {"heavy": "HHHH", "b2m": "MMMM"}}
    man = drv.emit_manifest(tmp_path, ["a"], truth, anns, panel,
                            tcr_seqs, mhc_seqs, k=1, samples=3)
    assert man["a"]["novel"] is True
    assert set(man["a"]["epitopes"]) == {"GILGFVFTL", "NLVPMVATV"}
    assert (tmp_path / "manifest.json").exists()
    for p in man["a"]["epitopes"].values():
        assert Path(p).exists()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./.venv/bin/python -m pytest tests/test_benchmark_driver.py -v`
Expected: FAIL with `ModuleNotFoundError: run_benchmark_arm`.

- [ ] **Step 3: Write minimal implementation**

```python
# scripts/run_benchmark_arm.py
"""Structure-vs-sequence retrieval benchmark driver.

Pre-fold half: select seed TCRs in one HLA (novel-first), build cognate+decoy
constructs, emit Protenix inputs + manifest for the user-driven Colab notebook.
Post-fold half (score_manifest) is added in the scoring task.

Usage (pre-fold):
  python scripts/run_benchmark_arm.py emit <dextramer_dir> <out_dir> \
      --hla 'HLA-A*02:01' --n 12 --k 3 --samples 3
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from rep2struct import benchmark as bm
from rep2struct.seqs import build_tcr_seqs, build_mhc_seqs
from run_validation_arm import labeled_clonotypes

def select_seed_tcrs(clonotypes, truth, annotations, hla, n, prefer_novel=True):
    dist = {a.clonotype_id: getattr(a, "tcrdist", None) for a in annotations}
    cands = [c for c in clonotypes if truth.get(c.id, (None, None))[1] == hla]
    cands.sort(key=lambda c: (0 if bm.is_novel(dist.get(c.id)) and prefer_novel else 1,
                              -c.size, c.id))
    return [c.id for c in cands[:n]]

def emit_manifest(out_dir, selected, truth, annotations, panel,
                  tcr_seqs, mhc_seqs, k, samples):
    out_dir = Path(out_dir)
    (out_dir / "constructs").mkdir(parents=True, exist_ok=True)
    by_id = {c: c for c in selected}
    dist = {a.clonotype_id: getattr(a, "tcrdist", None) for a in annotations}
    from rep2struct.schema import Clonotype
    clono_by_id = {}
    manifest = {}
    for cid in selected:
        cognate, hla = truth[cid]
        decoys = bm.decoys_for(cognate, hla, panel, k)
        # a minimal Clonotype carrier is enough for build_construct (uses .id)
        clono = Clonotype(cid, "", "", "", "", 0)
        jobs = bm.build_panel_constructs(clono, cognate, hla, decoys,
                                         tcr_seqs, mhc_seqs)
        eps = {}
        for ep, job in jobs.items():
            fp = out_dir / "constructs" / f"{cid}__{ep}.fasta"
            fp.write_text(job.construct_fasta)
            eps[ep] = str(fp)
        manifest[cid] = {"cognate": cognate, "hla": hla, "decoys": decoys,
                         "epitopes": eps, "novel": bm.is_novel(dist.get(cid)),
                         "tcrdist": dist.get(cid), "samples": samples}
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest

from run_validation_arm import nearest_cache, annotations_from_cache

def _load_truth_and_anns(dextramer_dir, hla):
    """Returns (clonotypes, truth, annotations). Annotations restricted to the
    target-HLA labeled clonotypes so we do not fire thousands of live queries."""
    clons, labels, hlas = labeled_clonotypes(
        f"{dextramer_dir}/donor1_all_contig_annotations.csv",
        f"{dextramer_dir}/donor1_binarized_matrix.csv")
    truth = {cid: (labels[cid], hlas[cid]) for cid in labels}
    in_hla = [c for c in clons if truth.get(c.id, (None, None))[1] == hla]
    cache = nearest_cache(in_hla)
    anns = annotations_from_cache(in_hla, cache)
    return clons, truth, anns, in_hla

def _emit_cmd(args):
    clonotypes, truth, anns, in_hla = _load_truth_and_anns(args.dextramer_dir, args.hla)
    counts = bm.per_hla_novel_counts(in_hla, truth, anns)
    print(json.dumps(counts.get(args.hla, {}), indent=2))
    selected = select_seed_tcrs(in_hla, truth, anns, args.hla, args.n)
    panel = bm.panel_epitopes(truth)
    sel_clonos = [c for c in clonotypes if c.id in set(selected)]
    tcr_seqs = build_tcr_seqs(sel_clonos)
    mhc_seqs = build_mhc_seqs([args.hla])
    emit_manifest(args.out_dir, selected, truth, anns, panel,
                  tcr_seqs, mhc_seqs, args.k, args.samples)
    print(f"emitted {len(selected)} TCRs x (1+{args.k}) constructs to {args.out_dir}")

def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("emit")
    e.add_argument("dextramer_dir"); e.add_argument("out_dir")
    e.add_argument("--hla", required=True); e.add_argument("--n", type=int, default=12)
    e.add_argument("--k", type=int, default=3); e.add_argument("--samples", type=int, default=3)
    e.set_defaults(func=_emit_cmd)
    args = ap.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run test to verify it passes**

Run: `./.venv/bin/python -m pytest tests/test_benchmark_driver.py -v`
Expected: PASS (2 tests).

- [ ] **Step 5: Commit**

```bash
git add scripts/run_benchmark_arm.py tests/test_benchmark_driver.py
git commit -m "feat(benchmark): seed selection + construct/manifest emission driver

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### MANUAL FOLD GATE (between Task 6 and Task 7)

Not an automated step. After Task 6 emits `constructs/*.fasta` + `manifest.json`:
1. Claude converts constructs to Protenix inputs via `rep2struct.tools.protenix_inputs.build` and builds the Colab notebook with `scripts/build_colab_notebook.py` (existing, user-driven pattern).
2. The user runs the notebook cells on Colab (GPU) and returns the per-construct Protenix output `.cif` sample files into `<out_dir>/folds/<clonotype_id>__<epitope>/`.
3. Proceed to Task 7 scoring once folds are present.

Chunk the batch to respect Colab wall-time; keep the seed N small.

---

### Task 7: Scoring + report (post-fold half)

Read the folds, run retrieval + both baselines + statistics, write an honest HTML/markdown report. Tested with a synthetic folds tree built from the cif fixtures.

**Files:**
- Modify: `scripts/run_benchmark_arm.py`
- Modify: `src/rep2struct/benchmark.py` (add the aggregation function)
- Test: `tests/test_benchmark_driver.py`

**Interfaces:**
- Consumes: `benchmark.contact_by_epitope`, `benchmark.retrieval_result`, `benchmark.plddt_by_epitope`, `benchmark.sequence_baseline_top1`, `benchmark.auroc`, `benchmark.bootstrap_ci`, `benchmark.permutation_p`.
- Produces:
  - `benchmark.evaluate(manifest, folds_root, annotations) -> dict` — per-TCR `{contact_top1, plddt_top1, seq_top1, cognate_contact, decoy_contacts}`, plus stratified aggregates (`overall`, `novel`, `leaked`) each with `{n, chance, contact: {top1, ci, p, auroc}, plddt: {top1, ci, p}, seq: {top1}}`.
  - `scripts/run_benchmark_arm.py::score_manifest(out_dir, dextramer_dir) -> dict` and a `score` subcommand; writes `benchmark_report.md`.

- [ ] **Step 1: Write the failing test**

```python
# append to tests/test_benchmark_driver.py
from rep2struct import benchmark as bm2

def _make_folds(root, cid, cognate, decoy, fix):
    for ep, src in [(cognate, fix/"cognate_min.cif"), (decoy, fix/"scramble_min.cif")]:
        d = root / "folds" / f"{cid}__{ep}"
        d.mkdir(parents=True)
        (d / "sample_0.cif").write_bytes(src.read_bytes())

def test_evaluate_stratifies_and_scores(tmp_path):
    fix = Path(__file__).parent / "fixtures"
    cid = "a"
    _make_folds(tmp_path, cid, "GILGFVFTL", "NLVPMVATV", fix)
    manifest = {cid: {"cognate": "GILGFVFTL", "hla": "HLA-A*02:01",
                      "decoys": ["NLVPMVATV"],
                      "epitopes": {"GILGFVFTL": "", "NLVPMVATV": ""},
                      "novel": True, "tcrdist": None, "samples": 1}}
    from rep2struct.schema import Annotation
    anns = [Annotation(cid, False, "unannotatable", tcrdist=None)]
    out = bm2.evaluate(manifest, tmp_path / "folds", anns)
    assert out["novel"]["n"] == 1
    assert out["novel"]["contact"]["top1"] in (0.0, 1.0)
    assert out["novel"]["seq"]["top1"] == 0.0        # novel TCR unannotated
    assert "ci" in out["novel"]["contact"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./.venv/bin/python -m pytest tests/test_benchmark_driver.py -k evaluate -v`
Expected: FAIL with `AttributeError: ... 'evaluate'`.

- [ ] **Step 3: Write minimal implementation**

```python
# append to src/rep2struct/benchmark.py
from pathlib import Path as _Path

def _paths_by_epitope(manifest_entry, folds_root, cid):
    out = {}
    for ep in manifest_entry["epitopes"]:
        d = _Path(folds_root) / f"{cid}__{ep}"
        out[ep] = sorted(str(p) for p in d.glob("*.cif")) if d.exists() else []
    return out

def _strata_stats(rows, chance):
    def agg(key, use_auroc):
        hits = [r[key] for r in rows]
        pt, lo, hi = bootstrap_ci(hits)
        d = {"top1": pt, "ci": [lo, hi], "p": permutation_p(hits, chance)}
        if use_auroc:
            d["auroc"] = auroc([(r["cognate_contact"], r["decoy_contacts"])
                                for r in rows if r["cognate_contact"] is not None])
        return d
    return {"n": len(rows), "chance": chance,
            "contact": agg("contact_top1", True),
            "plddt": agg("plddt_top1", False),
            "seq": {"top1": (sum(r["seq_top1"] for r in rows) / len(rows)) if rows else 0.0}}

def evaluate(manifest, folds_root, annotations):
    seq_ep = {a.clonotype_id: getattr(a, "epitope", None) for a in annotations}
    per = {}
    for cid, ent in manifest.items():
        pbe = _paths_by_epitope(ent, folds_root, cid)
        contacts = contact_by_epitope(pbe)
        plddts = plddt_by_epitope(pbe)
        cog = ent["cognate"]
        decoy_c = [contacts[e] for e in ent["decoys"] if contacts.get(e) is not None]
        per[cid] = {
            "novel": ent["novel"],
            "contact_top1": 1.0 if retrieval_result(contacts, cog)["top1"] else 0.0,
            "plddt_top1": 1.0 if retrieval_result(plddts, cog)["top1"] else 0.0,
            "seq_top1": 1.0 if sequence_baseline_top1(seq_ep.get(cid), cog) else 0.0,
            "cognate_contact": contacts.get(cog),
            "decoy_contacts": decoy_c,
        }
    rows = list(per.values())
    def chance(rs):
        ks = [1 + len(manifest[c]["decoys"]) for c in manifest]
        return sum(1.0 / k for k in ks) / len(ks) if ks else 0.0
    out = {"per_tcr": per,
           "overall": _strata_stats(rows, chance(rows)),
           "novel": _strata_stats([r for r in rows if r["novel"]], chance(rows)),
           "leaked": _strata_stats([r for r in rows if not r["novel"]], chance(rows))}
    return out
```

```python
# append to scripts/run_benchmark_arm.py (new subcommand)
def score_manifest(out_dir, dextramer_dir):
    out_dir = Path(out_dir)
    manifest = json.loads((out_dir / "manifest.json").read_text())
    clons, labels, hlas = labeled_clonotypes(
        f"{dextramer_dir}/donor1_all_contig_annotations.csv",
        f"{dextramer_dir}/donor1_binarized_matrix.csv")
    sel = [c for c in clons if c.id in manifest]
    anns = annotations_from_cache(sel, nearest_cache(sel))   # B1 sequence baseline
    result = bm.evaluate(manifest, out_dir / "folds", anns)
    lines = ["# Structure-vs-sequence benchmark\n"]
    for strat in ("overall", "novel", "leaked"):
        s = result[strat]
        lines.append(f"## {strat} (n={s['n']}, chance={s['chance']:.3f})")
        c = s["contact"]
        lines.append(f"- contact Top-1 {c['top1']:.2f} CI[{c['ci'][0]:.2f},{c['ci'][1]:.2f}] "
                     f"p={c['p']:.4f} AUROC={c.get('auroc')}")
        lines.append(f"- pLDDT Top-1 {s['plddt']['top1']:.2f} p={s['plddt']['p']:.4f}")
        lines.append(f"- sequence Top-1 {s['seq']['top1']:.2f}\n")
    (out_dir / "benchmark_report.md").write_text("\n".join(lines))
    return result

def _score_cmd(args):
    print(json.dumps(score_manifest(args.out_dir, args.dextramer_dir)["novel"], indent=2))
```

Wire the `score` subcommand in `main()`:

```python
    s = sub.add_parser("score")
    s.add_argument("out_dir"); s.add_argument("dextramer_dir")
    s.set_defaults(func=_score_cmd)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `./.venv/bin/python -m pytest tests/test_benchmark_driver.py -v && ./.venv/bin/python -m pytest -q`
Expected: PASS (all benchmark tests + full suite green).

- [ ] **Step 5: Commit**

```bash
git add scripts/run_benchmark_arm.py src/rep2struct/benchmark.py tests/test_benchmark_driver.py
git commit -m "feat(benchmark): fold scoring, stratified stats, honest report

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Self-Review

**Spec coverage:**
- Claim / retrieval framing -> Tasks 3, 7. ✓
- Free real negatives (decoys from panel) -> Task 2. ✓
- Novelty stratification + leakage guard (TCRdist > 1) -> Tasks 1, 7 (`is_novel`). ✓
- Metric reused unchanged (`qc.ensemble_contact`) -> Task 3. ✓
- Baselines B0 chance, B1 tcrdist, B2 pLDDT -> Tasks 4, 7. ✓
- Reuse released R2S modules so paper cites the tool -> every task imports from `rep2struct.*`. ✓
- Bootstrap CI + permutation test -> Task 5. ✓
- Seed scope one HLA / ~10-15 TCRs / k=3 / 3 samples -> Task 6 driver args. ✓
- Threats: leakage (Task 1), HLA/peptide confound (Task 2 `decoys_for` same-HLA first), small N (Task 5 stats), metric frozen (Global Constraints). ✓
- User-driven Colab fold path -> MANUAL FOLD GATE. ✓
- Class I only / Protenix only / no new model -> Global Constraints, no task adds a model. ✓

**Placeholder scan:** No TBD/TODO; every code step has complete code. ✓

**Type consistency:** `truth` is `dict[id -> (epitope, hla)]` throughout; `tcrdist` read via `getattr(a, "tcrdist", None)` consistently; `contact_by_epitope`/`plddt_by_epitope`/`retrieval_result`/`auroc`/`bootstrap_ci`/`permutation_p` signatures match between definition and callers in `evaluate`/`score_manifest`. ✓

**Known follow-ups (preprint, out of seed scope):** multi-HLA and multi-donor; contact+pLDDT combined ranker; cross-reactivity screening of decoys; epitope-matched strata. Flagged, not silently dropped.

---

## Adversarial-review revisions (BINDING — supersede the tasks above where they conflict)

Two adversarial reviews (science + implementation) ran on the plan. The changes below are mandatory and contain complete code. Where a revision changes a function defined above, use the version here.

### R0. Global additions

- Constructs per TCR are now **cognate + k same-HLA decoys + 1 scramble** (a seeded composition-preserving shuffle of the cognate). The scramble enables the within-pair contrast.
- The null for retrieval is the **empirical TCR-blind predictor** (rank by per-peptide mean contact), tested with a **label-permutation** test — NOT `1/(k+1)` alone.
- `retrieval_result` must treat a `None` cognate as a non-call (not a win) and a tie as a non-win.
- B2 baseline is **true CDR3β-residue pLDDT**, not peptide-chain pLDDT.
- Freeze the raw contact metric definition, but the *headline* is decided only after the TCR-blind and scramble controls; do not read Top-1 vs `1/(k+1)` as the result.

### R1. Task 2 revision — same-HLA-only decoys + scramble builder (replaces `decoys_for`, adds `scramble_peptide`)

```python
# src/rep2struct/benchmark.py  (replace decoys_for; add scramble_peptide)
import random as _random

def decoys_for(cognate, hla, panel, k):
    # same-HLA ONLY: cross-HLA decoys reintroduce the HLA-geometry confound (review M1)
    same = sorted(ep for (ep, h) in panel if h == hla and ep != cognate)
    return same[:k]

def scramble_peptide(cognate, seed=0):
    # composition-preserving shuffle; deterministic; retry so it differs from cognate
    chars = list(cognate)
    rng = _random.Random(f"{cognate}:{seed}")
    for _ in range(20):
        rng.shuffle(chars)
        s = "".join(chars)
        if s != cognate:
            return s
    return "".join(chars)
```

Add tests:

```python
# tests/test_benchmark.py
def test_decoys_same_hla_only_no_cross_fill():
    panel = [("GILGFVFTL","HLA-A*02:01"), ("KLGGALQAK","HLA-A*03:01")]
    assert bm.decoys_for("GILGFVFTL", "HLA-A*02:01", panel, k=3) == []  # no cross-HLA

def test_scramble_preserves_composition_and_differs():
    s = bm.scramble_peptide("GILGFVFTL", seed=1)
    assert sorted(s) == sorted("GILGFVFTL") and s != "GILGFVFTL"
```

Update `test_decoys_same_hla_first` expectation: with the cross-HLA fill removed, `test_decoys_fill_from_other_hla_when_short` is deleted (its behaviour is gone); `test_decoys_same_hla_first` still passes (same-HLA case unchanged).

`build_panel_constructs` also emits the scramble under key `"__scramble__"`:

```python
# src/rep2struct/benchmark.py  (extend build_panel_constructs)
def build_panel_constructs(clonotype, cognate, hla, decoys, tcr_seqs, mhc_seqs, scramble_seed=0):
    jobs = {}
    peptides = {ep: ep for ep in [cognate, *decoys]}
    peptides["__scramble__"] = scramble_peptide(cognate, scramble_seed)
    for key, pep in peptides.items():
        ann = Annotation(clonotype_id=clonotype.id, annotatable=True,
                         confidence_tier="benchmark", epitope=pep, hla=hla)
        jobs[key] = build_construct(clonotype, ann, tcr_seqs, mhc_seqs)
    return jobs
```

Adjust `test_build_panel_constructs_keys_and_peptides` to expect keys `{"GILGFVFTL","NLVPMVATV","__scramble__"}` and that `jobs["__scramble__"].construct_fasta` contains a permutation of the cognate.

### R2. Task 3 revision — None/tie-safe `retrieval_result` (replaces the version above)

```python
# src/rep2struct/benchmark.py  (replace retrieval_result)
def retrieval_result(contacts, cognate):
    # exclude the scramble key from retrieval; it is a separate contrast
    scored = {e: v for e, v in contacts.items() if e != "__scramble__"}
    cval = scored.get(cognate)
    valid = [v for v in scored.values() if v is not None]
    ranked = sorted(scored, key=lambda e: (-1.0 if scored[e] is None else scored[e]),
                    reverse=True)
    if cval is None or not valid:
        top1 = False                      # no-call, never a win (review code-F1)
    else:
        best = max(valid)
        n_at_best = sum(1 for v in valid if v == best)
        top1 = (cval == best) and (n_at_best == 1)   # ties are NOT wins
    return {"ranked": ranked, "top1": top1, "cognate_contact": cval}
```

Add tests:

```python
# tests/test_benchmark.py
def test_retrieval_none_cognate_is_not_a_win():
    assert bm.retrieval_result({"COG": None, "D": 5.0}, "COG")["top1"] is False

def test_retrieval_all_none_is_not_a_win():
    assert bm.retrieval_result({"COG": None, "D": None}, "COG")["top1"] is False

def test_retrieval_tie_is_not_a_win():
    assert bm.retrieval_result({"COG": 5.0, "D": 5.0}, "COG")["top1"] is False
```

### R3. Task 4 revision — CDR3β-residue pLDDT (replaces `model_plddt`/`plddt_by_epitope`)

```python
# src/rep2struct/benchmark.py  (replace model_plddt + plddt_by_epitope)
def _residue_bfactors(cif_path, chain_id):
    from Bio.PDB import MMCIFParser
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = next(MMCIFParser(QUIET=True).get_structure("x", str(cif_path)).get_models())
    for ch in m:
        if ch.id == chain_id:
            # per-residue mean B-factor (pLDDT column), residue order = input seq order
            return [float(np.mean([a.get_bfactor() for a in r])) for r in ch]
    return None

def model_cdr3b_plddt(cif_path, chain_b_seq, cdr3b):
    """Mean pLDDT over the CDR3beta residues, located as a substring of chain B."""
    if not chain_b_seq or not cdr3b:
        return None
    start = chain_b_seq.find(cdr3b)
    if start < 0:
        return None
    bfs = _residue_bfactors(cif_path, "B")
    if bfs is None or start + len(cdr3b) > len(bfs):
        return None
    return mean_confidence(bfs[start:start + len(cdr3b)])

def cdr3b_plddt_by_epitope(paths_by_epitope, chain_b_seq, cdr3b):
    out = {}
    for ep, paths in paths_by_epitope.items():
        vals = [v for v in (model_cdr3b_plddt(p, chain_b_seq, cdr3b) for p in paths)
                if v is not None]
        out[ep] = float(np.median(vals)) if vals else None
    return out
```

Replace the Task 4 pLDDT tests with:

```python
# tests/test_benchmark.py
def test_cdr3b_plddt_none_when_substring_absent():
    assert bm.model_cdr3b_plddt(str(FIX/"cognate_min.cif"), "AAAA", "ZZZZ") is None

def test_cdr3b_plddt_returns_float_when_located():
    # chain B of cognate_min.cif is a poly sequence in the fixture; if its residue
    # count >= len(query) and the query matches, we get a float, else None
    v = bm.model_cdr3b_plddt(str(FIX/"cognate_min.cif"), "BBBB", "BB")
    assert v is None or isinstance(v, float)
```

Note (review code-F4): whether Protenix writes pLDDT into the `.cif` B-factor column is **unverified**. Before trusting B2, parse one real Protenix output at the fold gate and confirm the B-factor column varies per residue (fixtures carry placeholder 1.0). If it does not, B2 is dropped, not silently reported.

### R4. Task 5 addition — TCR-blind null, label-permutation, paired contrast

```python
# src/rep2struct/benchmark.py  (append)
def tcr_blind_prediction(per_tcr_contacts):
    """Constant predictor: the peptide with the highest mean contact across TCRs.
    per_tcr_contacts: list of dict epitope->contact|None (scramble key excluded)."""
    sums, counts = {}, {}
    for d in per_tcr_contacts:
        for ep, v in d.items():
            if v is None:
                continue
            sums[ep] = sums.get(ep, 0.0) + v
            counts[ep] = counts.get(ep, 0) + 1
    if not sums:
        return None
    means = {ep: sums[ep] / counts[ep] for ep in sums}
    return max(means, key=means.get)

def tcr_blind_accuracy(per_tcr_contacts, cognates):
    pred = tcr_blind_prediction(per_tcr_contacts)
    if pred is None:
        return 0.0
    return sum(1 for cog in cognates if cog == pred) / len(cognates) if cognates else 0.0

def label_permutation_p(observed_top1_mean, per_tcr_contacts, cognates,
                        n_perm=10000, seed=0):
    """p that observed Top-1 > a predictor under random reassignment of cognate labels."""
    rng = _random.Random(seed)
    cogs = list(cognates)
    ge = 0
    for _ in range(n_perm):
        perm = cogs[:]
        rng.shuffle(perm)
        hits = 0
        for d, cog in zip(per_tcr_contacts, perm):
            r = retrieval_result({**d, "__scramble__": None}, cog)
            hits += 1 if r["top1"] else 0
        if hits / len(cogs) >= observed_top1_mean:
            ge += 1
    return (ge + 1) / (n_perm + 1)

def paired_contrast(pairs, seed=0, n_boot=2000):
    """pairs: list of (cognate_contact, scramble_contact), Nones dropped.
    Returns dict{n, frac_cognate_higher, mean_delta, ci_delta}."""
    vals = [(c, s) for (c, s) in pairs if c is not None and s is not None]
    if not vals:
        return {"n": 0, "frac_cognate_higher": None, "mean_delta": None, "ci_delta": [None, None]}
    deltas = [c - s for (c, s) in vals]
    frac = sum(1 for d in deltas if d > 0) / len(deltas)
    rng = _random.Random(seed)
    boots = []
    for _ in range(n_boot):
        boots.append(sum(deltas[rng.randrange(len(deltas))] for _ in range(len(deltas))) / len(deltas))
    boots.sort()
    return {"n": len(deltas), "frac_cognate_higher": frac,
            "mean_delta": sum(deltas) / len(deltas),
            "ci_delta": [boots[int(0.025 * n_boot)], boots[min(int(0.975 * n_boot), n_boot - 1)]]}
```

Tests:

```python
# tests/test_benchmark.py
def test_tcr_blind_accuracy_constant_predictor():
    contacts = [{"A": 9.0, "B": 1.0}, {"A": 8.0, "B": 2.0}]  # A always highest
    assert bm.tcr_blind_prediction(contacts) == "A"
    assert bm.tcr_blind_accuracy(contacts, ["A", "B"]) == 0.5  # blind picks A for both

def test_paired_contrast_cognate_higher():
    r = bm.paired_contrast([(10.0, 2.0), (8.0, 3.0)])
    assert r["frac_cognate_higher"] == 1.0 and r["mean_delta"] > 0

def test_paired_contrast_empty():
    assert bm.paired_contrast([(None, 1.0)])["n"] == 0
```

### R5. Task 6 revision — `select_seed_tcrs` seeded-random within stratum; driver `standardize_alleles`; manifest stores cdr3b/chain_b_seq/scramble; drop dead code

```python
# scripts/run_benchmark_arm.py  (replace select_seed_tcrs)
import random as _random
def select_seed_tcrs(clonotypes, truth, annotations, hla, n, prefer_novel=True, seed=0):
    dist = {a.clonotype_id: getattr(a, "tcrdist", None) for a in annotations}
    cands = [c for c in clonotypes if truth.get(c.id, (None, None))[1] == hla]
    novel = [c for c in cands if bm.is_novel(dist.get(c.id))]
    leaked = [c for c in cands if not bm.is_novel(dist.get(c.id))]
    rng = _random.Random(seed)
    rng.shuffle(novel); rng.shuffle(leaked)          # random within stratum (review M3)
    ordered = (novel + leaked) if prefer_novel else (leaked + novel)
    return [c.id for c in ordered[:n]]
```

`_load_truth_and_anns` must call `standardize_alleles` on the selected clonotypes before `build_tcr_seqs` (review code-F3), matching `run_validation_arm.main`:

```python
# scripts/run_benchmark_arm.py  (in _emit_cmd, before build_tcr_seqs)
    from rep2struct.ingest import standardize_alleles
    sel_clonos = [c for c in clonotypes if c.id in set(selected)]
    standardize_alleles(sel_clonos)                  # sets trav_allele/trbv_allele
    tcr_seqs = build_tcr_seqs(sel_clonos)
```

`emit_manifest` stores `cdr3b` and `chain_b_seq` (needed for CDR3β pLDDT), emits the scramble construct, and drops the dead `by_id`/`clono_by_id` lines. Replace its body loop:

```python
# scripts/run_benchmark_arm.py  (emit_manifest, revised loop; also accept clono_by_id map)
def emit_manifest(out_dir, selected, truth, annotations, panel,
                  tcr_seqs, mhc_seqs, k, samples, clono_by_id):
    out_dir = Path(out_dir)
    (out_dir / "constructs").mkdir(parents=True, exist_ok=True)
    dist = {a.clonotype_id: getattr(a, "tcrdist", None) for a in annotations}
    manifest = {}
    for cid in selected:
        cognate, hla = truth[cid]
        decoys = bm.decoys_for(cognate, hla, panel, k)
        clono = clono_by_id[cid]                      # real clonotype (has cdr3b)
        jobs = bm.build_panel_constructs(clono, cognate, hla, decoys,
                                         tcr_seqs, mhc_seqs)
        eps = {}
        for key, job in jobs.items():
            fp = out_dir / "constructs" / f"{cid}__{key}.fasta"
            fp.write_text(job.construct_fasta)
            eps[key] = str(fp)
        manifest[cid] = {"cognate": cognate, "hla": hla, "decoys": decoys,
                         "epitopes": eps, "novel": bm.is_novel(dist.get(cid)),
                         "tcrdist": dist.get(cid), "samples": samples,
                         "cdr3b": clono.cdr3b, "chain_b_seq": tcr_seqs[cid]["B"]}
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest
```

`_emit_cmd` passes `clono_by_id = {c.id: c for c in sel_clonos}` to `emit_manifest`, and the driver test updates its `emit_manifest(...)` call to pass a `clono_by_id={"a": clonos[0]}` mapping and to assert the manifest entry carries `cdr3b` and `chain_b_seq`, and that `epitopes` includes `"__scramble__"`.

### R6. Task 7 revision — per-stratum chance, TCR-blind null, paired contrast, CDR3β pLDDT, drop missing-cognate TCRs (replaces `evaluate` + report)

```python
# src/rep2struct/benchmark.py  (replace evaluate)
def evaluate(manifest, folds_root, annotations):
    seq_ep = {a.clonotype_id: getattr(a, "epitope", None) for a in annotations}
    per = {}
    for cid, ent in manifest.items():
        pbe = _paths_by_epitope(ent, folds_root, cid)
        contacts = contact_by_epitope(pbe)             # includes "__scramble__"
        cog = ent["cognate"]
        if contacts.get(cog) is None:
            continue                                   # no cognate fold -> drop, do not score as win
        panel_contacts = {e: v for e, v in contacts.items() if e != "__scramble__"}
        plddts = cdr3b_plddt_by_epitope(
            {e: p for e, p in pbe.items() if e != "__scramble__"},
            ent.get("chain_b_seq"), ent.get("cdr3b"))
        decoy_c = [contacts[e] for e in ent["decoys"] if contacts.get(e) is not None]
        per[cid] = {
            "novel": ent["novel"],
            "panel_contacts": panel_contacts,
            "contact_top1": 1.0 if retrieval_result(contacts, cog)["top1"] else 0.0,
            "plddt_top1": 1.0 if retrieval_result(plddts, cog)["top1"] else 0.0,
            "seq_top1": 1.0 if sequence_baseline_top1(seq_ep.get(cid), cog) else 0.0,
            "cognate_contact": contacts.get(cog),
            "scramble_contact": contacts.get("__scramble__"),
            "decoy_contacts": decoy_c,
            "cognate": cog,
            "n_panel": 1 + len(ent["decoys"]),
        }

    def strata(rows):
        if not rows:
            return {"n": 0}
        chance = sum(1.0 / r["n_panel"] for r in rows) / len(rows)
        contact_hits = [r["contact_top1"] for r in rows]
        pt, lo, hi = bootstrap_ci(contact_hits)
        blind_acc = tcr_blind_accuracy([r["panel_contacts"] for r in rows],
                                       [r["cognate"] for r in rows])
        perm_p = label_permutation_p(pt, [r["panel_contacts"] for r in rows],
                                     [r["cognate"] for r in rows])
        pc = paired_contrast([(r["cognate_contact"], r["scramble_contact"]) for r in rows])
        return {
            "n": len(rows), "chance": chance, "tcr_blind_acc": blind_acc,
            "contact": {"top1": pt, "ci": [lo, hi],
                        "p_vs_chance": permutation_p(contact_hits, chance),
                        "p_vs_blind": perm_p,
                        "auroc": auroc([(r["cognate_contact"], r["decoy_contacts"])
                                        for r in rows if r["cognate_contact"] is not None])},
            "plddt": {"top1": sum(r["plddt_top1"] for r in rows) / len(rows)},
            "seq": {"top1": sum(r["seq_top1"] for r in rows) / len(rows)},
            "scramble_contrast": pc,
        }

    rows = list(per.values())
    return {"per_tcr": per,
            "overall": strata(rows),
            "novel": strata([r for r in rows if r["novel"]]),
            "leaked": strata([r for r in rows if not r["novel"]])}
```

`score_manifest` report reads the new fields; replace its report block:

```python
# scripts/run_benchmark_arm.py  (score_manifest, report section)
    result = bm.evaluate(manifest, out_dir / "folds", anns)
    lines = ["# Structure-vs-sequence benchmark (seed = calibration pilot)\n"]
    for strat in ("overall", "novel", "leaked"):
        s = result[strat]
        if s.get("n", 0) == 0:
            lines.append(f"## {strat}: no TCRs\n"); continue
        c = s["contact"]; pc = s["scramble_contrast"]
        lines.append(f"## {strat} (n={s['n']}, naive chance={s['chance']:.3f}, "
                     f"TCR-blind null={s['tcr_blind_acc']:.3f})")
        lines.append(f"- contact Top-1 {c['top1']:.2f} CI[{c['ci'][0]:.2f},{c['ci'][1]:.2f}] "
                     f"p_vs_chance={c['p_vs_chance']:.4f} p_vs_blind={c['p_vs_blind']:.4f} "
                     f"AUROC={c['auroc']}")
        lines.append(f"- CO-PRIMARY within-pair cognate>scramble: "
                     f"frac={pc['frac_cognate_higher']} mean_delta={pc['mean_delta']} "
                     f"CI{pc['ci_delta']} (n_pairs={pc['n']})")
        lines.append(f"- CDR3b pLDDT Top-1 {s['plddt']['top1']:.2f} | "
                     f"sequence Top-1 {s['seq']['top1']:.2f}\n")
    (out_dir / "benchmark_report.md").write_text("\n".join(lines))
    return result
```

Update `test_evaluate_stratifies_and_scores`: build the synthetic folds with a `__scramble__` construct dir too (`_make_folds` writes `{cid}__GILGFVFTL`, `{cid}__NLVPMVATV`, and `{cid}____scramble__`), add `"cdr3b"` and `"chain_b_seq"` to the manifest entry, and assert `out["novel"]["scramble_contrast"]["n"] >= 0` and `out["novel"]["tcr_blind_acc"] in (0.0, 1.0)` alongside the existing checks. The headline assertion changes from reading `contact.top1` to also requiring `contact.p_vs_blind` to be present.

### R7. MANUAL FOLD GATE addition

At the gate, before scoring: parse ONE real Protenix `.cif` and confirm (a) five chains A-E, (b) the B-factor column varies per residue (pLDDT present) — else drop B2; (c) chain B sequence contains the clonotype's `cdr3b` substring — else CDR3β pLDDT returns None for that TCR (acceptable, logged). Fold the `__scramble__` construct alongside cognate+decoys for every TCR.
