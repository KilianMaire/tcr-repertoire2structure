# Structure Strategist: question-aware fold routing for Repertoire2Structure

Date: 2026-07-08
Status: design approved, pending implementation plan

## Problem

Folding in R2S is hardwired. `agents.py` exposes a single `fold-agent` whose only
instruction is to drive the Protenix Colab notebook; `prep_and_select` builds class I
constructs only. `fold.py` already accepts an injectable `fold_fn`, so the seam for
multiple backends exists but nothing uses it.

We want a specialist that reasons about the biological question and the construct
metadata, then picks and launches the right structure tool per group of clonotypes.
Protenix stays the default workhorse. This design covers the organization of that
specialist. It does not touch the Protenix fold code path itself (a separate session
currently runs it).

Baseline reference: the completed R2S run (orchestrator + fold/qc/report specialists
over 7 in-process MCP tools + Playwright) documented in the session artifact
claude.ai/code/artifact/0b154ba8. This design evolves that baseline; it does not
replace what works.

## User rules (binding)

- The strategist REASONS on the question and PROPOSES a tool. It is not a
  deterministic priority probe.
- Protenix is the default. On any question where no known tool is clearly better,
  stay on Protenix.
- Switch to a specialized tool only when the question justifies it: AF3 when the user
  has the weights and it helps; alphafold_finetune for "is this peptide presented?"
  (class I and II, only option for Der p 1 class II); MHC-Fine for the most precise
  class I peptide pose; TCRdock for the TCR:pMHC interface and V-domain anchoring.
- NEVER Boltz.
- The tool catalog is incomplete and the agent must know it. If a group seems to fall
  outside the validity domain of every wired tool, fall back to Protenix AND flag the
  group in the report with an honest reservation ("this case would be better served by
  an un-wired tool X").

## Decisions taken during brainstorming

1. Question source: BOTH. Default is data-driven inference from construct tags
   (current R2S entry preserved). An optional free-text question, if supplied, takes
   precedence and steers the choice. Backward compatible.
2. Granularity: PER HOMOGENEOUS GROUP. The strategist partitions constructs into
   homogeneous groups (by class / has-TCR / species / output need) and chooses one
   tool per group. Each group gets its own QC calibration and a traceable report row.
3. Out-of-catalog handling: wire the maximum number of tools now, so the case is rare.
   Residual guard = fallback Protenix + explicit report flag.
4. Agent topology: STRATEGIST PLANS, PER-TOOL EXECUTORS. One `structure-strategist`
   agent decides (judgement); it delegates each group to a dedicated executor agent
   per tool. Decision and execution are separated.
5. Tool-capability representation: DATA REGISTRY (approach A). A pure-data
   `StructureTool` registry is the single source of truth; a `list_structure_tools`
   tool feeds it to the strategist, which reasons over facts rather than guessing.
6. Compute target: COLAB-DRIVEN-VIA-BROWSER for all tools (one pattern). There is no
   local GPU in R2S; the GPU is remote (Colab PRO A100 / rented H100), piloted by the
   executor through Playwright MCP.
7. Monitor: DEFERRED. v1 keeps the executor-waits pattern. A shared lightweight
   monitor role (poll resumable markers, re-invoke on completion) is added only if
   runs become too long.
8. MSA: RESTORED as a pre-fold artifact, not a runtime call. Computed local-on-Mac if
   the DB is present, else a Colab CPU step, else MSA-free flagged fallback. This kills
   the remote-MSA-throttle failure mode. Applies to structure tools only.

## Architecture

```
orchestrator  (csv, run_dir, top_n, question?)
   |
   v
prep_and_select  -> tagged constructs {mhc_class, has_tcr, species, allele, output_needed}
   |
   v
structure-strategist            REASONS, does not fold
   |  1. partition tagged constructs into homogeneous groups
   |  2. per group: read list_structure_tools, pick one tool, justify in one sentence,
   |     or flag out-of-domain (fallback Protenix + reservation)
   |  3. delegate each group to the chosen tool's executor
   v
per-tool executor agents (one each)
   protenix-agent . af3-agent . tcrdock-agent . mhcfine-agent . affinetune-agent
   |  each pilots its own Colab notebook via Playwright, writes a normalized result
   |  via record_fold_result; relies on fold.py resumable done-markers
   v
qc-agent  (per-group calibration, output-type-aware) -> report-agent
```

Key organizational properties:

- The strategist knows no launch mechanics. It manipulates tool names and descriptors.
  All browser/Colab complexity lives in the executors.
- Each executor is a small single-responsibility agent: take a list of jobs, pilot its
  tool, call `record_fold_result` in a uniform format. Testable in isolation.
- `fold.py` does not change. Its `run_folds(jobs, fold_fn, run_state)` stays the
  resumable engine; each executor supplies its own `fold_fn`.

## Roles (who does what)

- Scientific analysis = deterministic CODE, not an LLM. Scoring (CDR3-Vbeta vs peptide
  contact, scramble null as ensemble mean, geometry/DockQ) lives in qc.py / qc_folds.py.
  The qc-agent calls the tool, the code computes, the agent only judges reliable/suspect.
- GPU piloting = the executor agent, via Playwright, on REMOTE compute. No local GPU.
  "Wiring a tool" therefore means writing its remote-drive adapter (a Colab notebook),
  not just a registry entry. Per-tool cost is uneven: Protenix notebook is ready;
  MHC-Fine and alphafold_finetune ship Colabs to adapt; TCRdock has no Colab (wrapper
  notebook needed); AF3 depends on gated weights.
- Monitoring = today the executor holds the browser session and polls until cells
  finish, plus fold.py `fold_{id}.done.txt` resumable markers. No separate watchdog in
  v1.

## Component: StructureTool registry

Single source of truth, pure data, in `src/rep2struct/structure_tools.py`.

```python
StructureTool(
    name="protenix",
    validity={"mhc_class": {1, 2}, "needs_tcr": None, "species": "any"},
    output_type="structure",          # structure | binding_score
    strengths="full 3-chain TCR-pMHC fold, default workhorse",
    limits="imposes canonical geometry even on non-binders (basis of skeptical QC)",
    colab_adapter="protenix_colab",   # notebook the executor pilots
    is_default=True,
)
```

v1 entries:

- `protenix`: default, class I+II, output structure.
- `af3`: structure, gated on user having weights.
- `mhcfine`: class I, precise peptide pose, output structure.
- `tcrdock`: TCR:pMHC, interface + V-domain, output structure.
- `affinetune`: class I+II, output binding_score.

Three fields consumed elsewhere:

- `validity`: what the strategist reads (via `list_structure_tools`) to reason. It
  does not guess capabilities, it reads them.
- `output_type`: honesty guard. `affinetune` is `binding_score`, so the report never
  renders it as a "fold". Encoded in the schema, not just prose.
- `colab_adapter`: the notebook name the matching executor pilots.

Adding a future tool = one entry + one notebook + one small executor agent. The
strategist prompt does not change.

## Component: QC calibration per group and report honesty

This is where a bad design would manufacture false science. Three rules.

1. Each group carries its OWN scramble calibration, computed with the SAME tool. The
   contact-distance null is only valid for the tool that produced it (a Protenix
   contact distance and an MHC-Fine one are not on the same scale). So per homogeneous
   group the executor folds cognate + scramble with the same tool, and the qc-agent
   calibrates the threshold on that null. No global threshold shared across tools.

2. `binding_score` tools do NOT go through geometric QC. `affinetune` returns no
   structure; there is nothing to measure in CDR3-peptide contact. Its group follows a
   separate QC path: the presentation score is calibrated against known
   binders/non-binders (or a shuffled-sequence score null), not against geometry. The
   qc-agent branches on the registry's `output_type`.

3. The report traces tool, evidence type, and calibration basis per group, and NEVER
   cross-compares raw numbers between groups. Two guards encoded in the schema (like
   the two existing honesty rules):
   - a `binding_score` result is labelled "predicted presentation", never "fold" or
     "structure";
   - no ranking/sorting of groups on raw distances (not comparable across tools); a
     verdict is compared only to its own calibration;
   - a group routed to fallback Protenix "out of domain" carries a visible reservation.

R2S's core rule is unchanged: a clean fold never confirms specificity, whatever the
tool. The skepticism is structural, not cosmetic.

## Component: MSA strategy

The baseline run dropped the MSA (`--use_msa false`) because the public remote MSA
server throttled and wedged the Colab kernel. That was a robustness workaround, not a
scientific preference. This design restores the MSA while removing the failure mode.

Principle: the failure was the remote MSA server AT RUNTIME, not the MSA itself. So the
MSA is computed OUTSIDE the fold runtime and embedded into the construct inputs, the
same way inputs are already embedded in the self-contained notebook. The piloted fold
notebook makes no throttlable network call for the MSA.

Execution order (first feasible wins):

1. Local on the user's machine: mmseqs2 against a reduced ColabFold-style DB
   (e.g. UniRef30) hosted on the 4TB drive. Produces an a3m per chain, cached in the
   run directory.
2. Colab CPU step: a separate pre-fold notebook runs mmseqs2 on a reduced DB on CPU,
   returns the a3m. Used when the local DB is absent.
3. MSA-free fallback: the current robust path, used only if 1 and 2 both fail, and
   flagged in the report as reduced-confidence.

Scope: MSA is a cross-cutting concern of `output_type = structure` tools (protenix,
af3, mhcfine, tcrdock). It does not apply to `affinetune` (binding_score).

Honesty note: the MSA mainly deepens framework and MHC/B2M chains (well conserved). It
does almost nothing for the hypervariable CDR3 that carries specificity. So it improves
structural quality but does NOT change the core caveat that a clean fold never confirms
specificity. The report must not present the MSA as a specificity gain.

Prerequisite to verify in the plan: whether mmseqs2 and a reduced DB are actually
present/installable on the Mac and the 4TB drive; if not, path 2 is the v1 default.

## Out of scope

- The Protenix fold code path itself (running in a separate session).
- Actually integrating gated AF3 weights or standing up a real H100 job runner. v1
  targets the routing organization + Colab adapters; a genuine remote H100 adapter is
  a later option.
- The shared monitor role (deferred until runs demand it).

## Backward compatibility

No question supplied + a single homogeneous class I TCR-pMHC group = the strategist
picks Protenix for the whole run, i.e. today's behaviour. The entry point gains an
optional `question` argument; everything else is preserved.
