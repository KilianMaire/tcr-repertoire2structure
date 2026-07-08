# Wiring tcrdock, mhcfine, affinetune: builders, two-layer QC, notebook scaffolds

Date: 2026-07-08
Status: design approved, pending implementation plan

## Problem

The structure-strategist routing feature (merged, main d408e91) knows five tools but
only Protenix has a working Colab notebook. The other executors honestly report
"not-run" when their adapter is absent. This design wires three of them: tcrdock (TCR
interface), mhcfine (precise class I pMHC pose), affinetune (is-this-peptide-presented
classifier). af3 is deferred: AF3 weights are under a non-redistributable license, so
the user cannot supply theirs for now.

The build splits into an offline-testable code layer (input builders, two-layer QC,
registry metadata, report) delivered and merged this round, and a live layer (the real
notebook cells, the Playwright drive, real folds and per-tool scramble calibration)
scaffolded now and validated together one pipeline at a time afterward.

Baseline reference: the merged routing feature
(docs/superpowers/specs/2026-07-08-structure-strategist-design.md). This design fills in
adapters against that architecture; it does not change the routing.

## Decisions taken during brainstorming

1. Class scope: CLASS I FIRST. All three tools wired for class I (aligned with
   Repertoire2Structure and TABLO, both human class I). affinetune runs class I this
   round. MHC II construct geometry and Der p 1 are deferred to a second round.
2. Notebook representation: ONE BUILDER PER TOOL, mirroring Protenix
   (build_protenix_inputs.py + build_colab_notebook.py). Each tool gets a
   `<tool>_inputs.py` (construct to the tool's input format) and a
   `<tool>_notebook.py` (self-contained notebook). Offline-testable builders;
   notebooks validated live.
3. mhcfine QC: BOTH an independent peptide-in-groove-vs-scramble control (primary
   skeptical verdict) AND the model confidence (pLDDT/pae) reported alongside, never as
   the sole basis. mhcfine produces a pMHC pose with no TCR, so the CDR3-peptide metric
   does not apply.
4. A COMMON generalist QC layer runs for every tool in addition to the tool-specific
   metric (see Two-layer QC).
5. Offline/live split as above; live pipelines validated one by one, together, after
   all code is written.

## Architecture

Per-tool file structure, mirroring Protenix, under a new `src/rep2struct/tools/`:

```
tcrdock_inputs.py      construct (cognate + scramble) -> TCRdock input format
tcrdock_notebook.py    -> self-contained Colab notebook (AF2 weights)   [LIVE cells scaffolded]
mhcfine_inputs.py      pMHC only (no TCR chains) -> MHC-Fine input
mhcfine_notebook.py    -> notebook (adapts their PyTorch Colab)          [LIVE cells scaffolded]
affinetune_inputs.py   class I pMHC -> alphafold_finetune input
affinetune_notebook.py -> notebook (adapts their JAX Colab)              [LIVE cells scaffolded]
```

Each `<tool>_inputs.py` reuses `seqs.py` (reconstructed TCR V domains, cached class I
HLA ectodomains) and produces BOTH the cognate and the scramble construct, exactly as
`build_protenix_inputs.py` already does, so each tool carries its own scramble null.

Canonical chain IDs: every notebook emits its CIF with canonical chain IDs so QC stays
uniform. A = TCR alpha, B = TCR beta, C = MHC heavy, D = beta2-microglobulin,
E = peptide. mhcfine has no TCR, so it emits C, D, E only.

## Component: registry qc_metric field

Add a `qc_metric` field to `StructureTool` and set it per tool:

- protenix: `cdr3_peptide`
- tcrdock: `cdr3_peptide`
- mhcfine: `peptide_groove`
- affinetune: `binding_score`
- af3: `cdr3_peptide` (unused until wired)

The QC selects its skeptical control from `qc_metric`, so the metric is chosen by the
tool, not guessed. `as_dicts()` includes `qc_metric` so it is visible to the strategist
and the report.

## Component: two-layer QC

**Layer 1, common generalist QC (`qc.common_checks`), runs for every tool.** A
tool-agnostic validity gate that returns a pass/fail plus a small set of comparable
descriptors:

- structure tools (protenix, tcrdock, mhcfine): expected chains present, coordinates
  finite (no NaN, no exploded geometry), chain lengths consistent with the input, no
  severe steric clash (minimum interatomic distance plausible), peptide present.
- binding tool (affinetune): score is finite and within the expected range ([0, 1]).

If the gate fails, the result is `qc_failed` with the reasons, and the layer-2 metric
does NOT run (measuring a contact on a broken structure is meaningless).

**Layer 2, tool-specific skeptical QC**, runs only if layer 1 passes, dispatched on
`qc_metric`:

- `cdr3_peptide` (protenix, tcrdock): CDR3-Vbeta vs peptide contact against that tool's
  own scramble null. Existing metric.
- `peptide_groove` (mhcfine): peptide vs MHC-groove-residue contact against a scramble
  null (primary verdict), plus the model confidence read and reported alongside, never
  the sole basis.
- `binding_score` (affinetune): `verdict_binding` vs a score null. Already coded.

`qc_structure` runs `common_checks` first (gate plus comparable descriptors reported),
then the `qc_metric` control.

## Component: report

- A new "validity (common)" column shows the layer-1 descriptors. These are
  tool-agnostic, so they MAY be shown side by side across tools: a common denominator.
- The layer-2 skeptical metrics (contact distances vs scramble) remain NOT comparable
  across tools: no cross-tool ranking, each verdict against its own calibration. The
  existing honesty rule holds.
- Honesty labels by tool: tcrdock is a structure (like Protenix, calibration basis
  noted); mhcfine is a POSE, never a proof of TCR recognition (it has no TCR);
  affinetune is a predicted presentation. Encoded via `qc_metric`/`output_type`, not
  just prose.

## Component: per-tool scramble and score calibration

The executor folds or scores cognate plus scramble with its own tool; the qc-agent
passes that group's own `scramble_threshold` per call (the per-call mechanism already
exists). A tcrdock threshold, an mhcfine threshold, and a Protenix threshold are
distinct and never shared.

## Offline vs live split

Delivered and merged this round (offline, TDD):

- `qc_metric` registry field.
- The three `<tool>_inputs.py` builders (cognate plus scramble, reusing seqs.py),
  tested on fixtures.
- `qc.common_checks` (layer 1), tested on the existing CIF fixtures (cognate_min.cif,
  scramble_min.cif, threechain_min.cif).
- The `peptide_groove` metric (layer 2, mhcfine) in qc.py.
- `qc_structure` wiring: common_checks gate, then qc_metric dispatch.
- Report: the common-validity column and the per-tool honesty labels.

Scaffolded now, validated live together one pipeline at a time:

- The real cells inside each `<tool>_notebook.py` that actually invoke the tool
  (their repo, weights, Colab). Scaffold plus explicit TODO markers.
- The Playwright drive sequence per executor against its notebook.
- Real folds and the per-tool scramble calibration numbers (need a GPU).
- The mhcfine model-confidence extraction format (pLDDT/pae depends on their output):
  the reader is coded, validated live.

Guard rail already in place: an executor reports "not-run" when its notebook or adapter
is unavailable and never fabricates output. Merging the offline code breaks nothing:
until a notebook is validated, that tool stays "not-run" and Protenix keeps running.

## Out of scope

- af3 (non-redistributable weights; deferred).
- Class II construct geometry and Der p 1 (deferred to a second round).
- Real MSA runners (still MSA-free, flagged reduced confidence).
- The Protenix fold code path (unchanged).

## Backward compatibility

Protenix keeps `qc_metric="cdr3_peptide"` and its existing path. The common_checks gate
runs before the existing cdr3_peptide metric, so a Protenix fold that passes validity
behaves as today. New builders and notebooks are additive; unwired tools stay "not-run".
