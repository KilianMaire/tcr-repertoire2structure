# Design: Structure recovers specificity where sequence cannot (R2S benchmark)

Date: 2026-07-09
Status: approved design, pre-plan

## Goal and audience

Turn R2S from an honest integration demo into a genuine, falsifiable scientific
contribution. Two tiers, in sequence:

1. **Seed (by 2026-07-13, hackathon):** a small, honest pilot with a stated N and
   confidence intervals. A null result is acceptable and still honest.
2. **Preprint (after):** the same code scaled to a full panel, multi-donor,
   pLDDT-combination analysis, and a proper held-out novel set.

The seed's code IS the preprint's code. Nothing is thrown away.

## The central claim

> On sequence-novel CD8 TCRs (no close TCRdist neighbour), the scramble-calibrated
> CDR3-peptide contact recovers the correct epitope from a candidate panel above
> chance, precisely the population where sequence annotation fails (de-leaked
> recall 0.08 on the same dextramer set).

This is the value R2S adds over `tcrdist`: it targets the unannotatable majority,
not the easy leaked minority.

## Why the negatives are free and real

The 10x 4-donor CD8 dextramer panel is a **closed set of epitopes**; each labelled
clonotype binds exactly one. So for a TCR whose cognate is epitope `e`, every other
panel epitope is a **verified non-binder** (a real decoy, not a shuffled peptide).
The binder/non-binder benchmark needs no wet-lab and no synthetic negatives.

## Task framing: retrieval

For each TCR:
1. Build cognate construct + `k` decoy constructs (decoys share the TCR's
   restricting HLA where possible, so we test peptide discrimination, not HLA
   geometry).
2. Fold each with Protenix (headless, scriptable; see Risks), several samples.
3. Score CDR3(Vbeta)-peptide heavy-atom contact, sample-averaged.
4. Rank the epitopes by contact.

**Primary readouts** (each with bootstrap CI + permutation test):
- Top-1 accuracy: is the cognate ranked #1?
- AUROC of cognate-vs-decoy contact.

**Decisive stratification:** de-leaked novel TCRs (TCRdist > 1 to any reference)
vs TCRs with a close neighbour. The headline is the novel stratum.

## Metric and the baselines it must beat

The contribution is the **metric**, not the folder. On the same folds:

- **B0 chance:** 1/(k+1).
- **B1 sequence (tcrdist):** the existing `annotate.annotate` output. On novel TCRs
  this is ~chance by construction; that is the gap being filled.
- **B2 CDR3 pLDDT reranking:** the 2025 prior-art signal (CDR3 pLDDT correlates with
  docking quality). The contact metric must beat OR complement it, else the result
  is not novel. Report contact, pLDDT, and their combination.
- **Scramble null:** kept as an internal per-pair calibration (cognate vs its own
  shuffled peptide), secondary to decoy retrieval.

Success condition: contact beats chance AND adds over pLDDT on the novel stratum.

## Reuse of released R2S modules (so the paper cites the tool)

The benchmark is implemented on top of the published package, not throwaway scripts.
Concretely it imports:

- `rep2struct.annotate.annotate` — B1 sequence baseline and the novelty strata
  (distance + tier are already returned).
- `rep2struct.foldprep.build_construct` — cognate and decoy construct assembly.
- `rep2struct.tools.protenix_inputs.build` — Protenix fold inputs.
- `rep2struct.qc.ensemble_contact` / `qc.score_model` — the contact metric (the
  contribution), reused verbatim from the QC step.
- `rep2struct.qc.mean_confidence` — the pLDDT-style confidence for the B2 baseline.

New code is additive and lives in the repo: a benchmark driver
(`scripts/run_benchmark_arm.py`) and, if shared logic warrants it, a thin
`rep2struct.benchmark` module (retrieval scoring, strata, bootstrap). Methods
section references the released package version. Code is public before submission.

## Scope

**Seed (2026-07-13):** one HLA (best-populated, likely A\*02:01), ~10-15 TCRs
balanced across ~3-4 epitopes, cognate + 3 decoys, 3 Protenix samples.
~120-240 folds. Report Top-1 / AUROC with CIs, stratified, framed as a pilot with
N stated. Class I only.

**Preprint (after):** full multi-HLA panel, all 4 donors, larger `k`, 5 samples,
pLDDT-combination analysis, held-out novel set. Class II deferred (out of scope for
both tiers unless the seed motivates it).

## Threats to validity (stated up front)

1. **Leakage:** reuse the existing TCRdist <= 1 guard to define the novel stratum;
   report raw and de-leaked.
2. **HLA/peptide confound:** decoys share the restricting HLA where possible.
3. **Small N:** bootstrap CIs + permutation test; no claims past the pilot.
4. **Metric arbitrariness / tuning-on-truth:** freeze the contact definition
   (heavy-atom cutoff, CDR3beta, sample-averaged) BEFORE looking at retrieval
   results.

## Key risks / open dependencies

- **Protenix headless at scale (the main feasibility risk):** the current path is
  Playwright-driven Colab, shown to be fragile (MSA server throttling wedged a
  batch). The benchmark needs Protenix driven programmatically (CLI/API on a
  rented or scripted-Colab A100). Confirming this drive path is the first
  implementation task; if it cannot be made reliable, fall back to a smaller N
  rather than switching engine (literature: AF3-class >= TCRdock for TCR-pMHC
  class I).
- **Enough novel TCRs in one HLA:** the pilot must contain enough de-leaked novel
  clonotypes in the chosen HLA to make the headline claim; verify before folding.

## Non-goals

- No new folding model. No class II. No inter-model benchmark (Protenix only).
- No wet-lab validation. No claim of binding affinity beyond contact-based ranking.
