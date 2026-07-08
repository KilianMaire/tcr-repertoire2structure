# Repertoire2Structure design

Working name, changeable. An agentic pipeline, orchestrated by Claude, that takes a
raw 10x single cell TCR repertoire and produces QC'd predicted TCR pMHC structures for
its top clonotypes, with honest specificity annotation and hallucination flagging.

Built for the "Built with Claude: Life Sciences" hackathon (Researcher track, deadline
2026 July 13). The demonstration is Claude acting as an autonomous bioinformatics
operator end to end: routing between MCP tools, driving a browser for GPU folding,
making skeptical QC judgments, and writing the interpretive report.

## Goal

A researcher hands the pipeline a 10x repertoire CSV. Claude agents return, fully
automatically, a report that links each top clonotype to a candidate epitope
specificity (with confidence), a predicted TCR pMHC structure, and a skeptical QC
verdict on whether that structure is trustworthy or a likely geometry hallucination.

## Two honesty rules, baked into the architecture

These are non negotiable and shape the output schemas, not just the prose.

1. **Specificity is annotation by similarity, never prediction.** A clonotype is
   annotated with a candidate epitope only when a TCRdist neighbor in a labeled
   reference is close enough. The output always carries the distance and a confidence
   tier, and clonotypes with no close neighbor are flagged `unannotatable`. We never
   force a label.

2. **A structure does not confirm specificity.** Protenix (AlphaFold3 reproduction)
   imposes canonical TCR pMHC docking geometry even on non binding sequences (session
   finding, reproduced across ColabFold and Protenix, scramble not discriminable from
   cognate). So the QC step is a skeptical judge: it flags a fold as `suspect` when its
   CDR3 to peptide contact discrimination is no better than the scramble control
   calibration. A clean fold is never treated as evidence that the annotated epitope is
   correct.

## Datasets

Both human, both CD8 oriented, both 10x format.

- **Validation arm (ground truth).** 10x Genomics 4 donor CD8 dextramer set (44 plex
  pMHC dextramers, roughly 150k cells, one known epitope per cell where the dextramer
  count binarizes positive). Format: `filtered_contig_annotations.csv` plus binarized
  dextramer counts. Used to MEASURE the annotation step: precision, recall, and
  unannotatable rate of TCRdist vs VDJdb annotation against the dextramer label.

- **Application arm (scale).** TABLO (Sureshchandra, Zenodo 10.5281/zenodo.13119615,
  CC BY 4.0): multimodal atlas of 5.7M human T cells, blood and tonsil, 10 donors,
  GEX plus TCR, no epitope labels. Used to RUN the full pipeline end to end on a real
  unlabeled repertoire: annotate what is annotatable, flag the rest, fold the top
  clonotypes.

## Architecture

A chain of Claude agents. Each stage has a narrow tool surface and a typed output that
the next stage consumes. State lives on disk (a run directory) so the chain is
resumable, which matters because the fold step depends on a remote GPU session that can
wedge.

```
10x CSV
  -> [0] Ingest and clonotype curation
  -> [1] Specificity annotation (honest)
  -> [2] Fold prep and MSA
  -> [3] Fold (full auto GPU, Protenix via Colab and Playwright)
  -> [4] Skeptical QC (hallucination flag)
  -> [5] Report
```

### Stage 0. Ingest and clonotype curation

- Input: 10x paired alpha beta contig annotations.
- Clonotype definition, fixed upfront: identical (TRAV gene, CDR3 alpha aa, TRBV gene,
  CDR3 beta aa). Cells collapsing to the same tuple are one clonotype; the cell count
  is the clonal size.
- Standardize V and J alleles through TCR Explorer `assign_tcr_alleles`.
- Output: a clonotype table with allele honored V and J, both CDR3s, and clonal size.
- Reuses: `stage55_contig_annotation.py` (10x parsing), TCR Explorer `assign_tcr_alleles`.

### Stage 1. Specificity annotation (honest)

- For each clonotype, run TCRdist against a labeled reference (VDJdb, and IEDB where
  paired) through TCR Explorer `find_similar_paired_tcrs` and `tcrdist_engine`.
- Emit the nearest labeled neighbor, its TCRdist distance, the candidate epitope and
  restricting HLA, and a confidence tier derived from the distance (for example: high,
  medium, low, unannotatable). Thresholds calibrated on the validation arm.
- Output per clonotype: `{epitope, hla, tcrdist, confidence_tier, annotatable: bool}`.
- Reuses: TCR Explorer MCP tools (live), `dossier_epitopes`.

### Stage 2. Fold prep and MSA

- Rank clonotypes by (confidence tier) times (clonal size); select top N. N small for
  the hackathon (target 8 to 12).
- Build the TCR pMHC construct, chains A to E: TCR alpha V domain, TCR beta V domain,
  MHC class I heavy chain, beta 2 microglobulin, peptide. (Class I here, since the demo
  is human CD8. This differs from the existing class II I A b construct and is a build
  item.)
- Precompute the MSA once per shared component and per TCR.
- Output: one Protenix input bundle per selected clonotype.
- Reuses: `stage52` and `stage53` input and precomputed MSA path, `msa.py`. Build item:
  class I construct assembly.

### Stage 3. Fold (full auto GPU)

- A Claude agent drives Playwright to a Colab notebook running Protenix
  (`protenix_base_default_v1.0.0`, 5 seeds), following the established resumable
  procedure: background execution survives disconnect, jobs skip when a `.done.txt`
  exists, downloads land in `~/.playwright-mcp/` and are renamed.
- Output: predicted CIF models per clonotype (25 models each at 5 seeds).
- Reuses: the existing Colab and Playwright folding procedure, `stage45` Protenix path.
  Build item: wrap the procedure as an autonomous agent step with retry and resume.

### Stage 4. Skeptical QC (hallucination flag)

- Score each model: DockQ calibrated geometry, CDR3 to peptide atom contacts, crossing
  and docking angle, peptide vs MHC contact balance.
- Compare CDR3 to peptide contact discrimination against the scramble control
  calibration. Verdict per clonotype: `reliable` when the cognate profile beats the
  scramble null, `suspect` otherwise.
- Output per clonotype: metrics plus `{qc_verdict, reason}`.
- Reuses: `stage41` DockQ calibration, `stage43` and `stage45` scoring, the session
  scramble control result. Build item: fold the scramble calibration into a per model
  verdict function.

### Stage 5. Report

- Auto generate figures (structure panels, geometry, confidence) and a report that ties
  clonotype to candidate specificity to structure to QC verdict, for both arms.
- For the validation arm, also emit the annotation performance table (precision, recall,
  unannotatable rate vs dextramer ground truth).
- Reuses: `build_report_figures.py`, `render_pymol.py`, `render_structures.py`.

## What Claude does (the hackathon hook)

Claude is not a wrapper around one model call. It orchestrates the chain: it routes
clonotypes through the MCP tools, decides which clonotypes are worth folding, drives the
browser to run the GPU fold, applies the skeptical QC judgment (including reading the
scramble calibration and deciding reliable vs suspect), and writes the interpretive
report. The demo shows this running from a raw CSV to a finished report with minimal
human input.

## What exists vs what we build

Roughly 70 percent reuse, 30 percent glue.

Exists: 10x contig parsing, TCR Explorer MCP (clonotype allele assignment, TCRdist,
paired similarity, epitope dossier, MSA), Protenix Colab and Playwright folding
procedure, DockQ calibration and geometry scoring, figure and report builders, the
scramble control result.

Build: the orchestration and run state layer; the class I construct assembly; the
autonomous fold agent wrapper with retry and resume; the honest annotation output schema
with confidence tiers and unannotatable flagging; the scramble calibrated QC verdict
function; the final cross arm report generator.

## Scope and YAGNI

In scope for the hackathon: one engine (Protenix), one clonotype definition, top N small,
the two named datasets, resumable fold, honest annotation and QC, a single report.

Out of scope (stretch, noted not built): multiple folding engines, class II support in
the same run, de novo epitope prediction, a hosted web UI, batch scale beyond top N.

## Error handling

- Ingest: malformed or unpaired 10x rows are dropped with a counted reason, never
  silently.
- Annotation: no close neighbor yields `unannotatable`, not a guessed label.
- Fold: a wedged Colab session is recoverable because the loop is resumable; the agent
  retries, and the demo keeps a replayable fold so it is not hostage to a live session.
- QC: a model that fails to parse or lacks the five chains is reported as `qc_failed`,
  distinct from `suspect`.

## Testing

- Stage 0: unit test clonotype collapsing on a small fixture, including allele
  standardization.
- Stage 1: on the validation arm, assert precision, recall, and unannotatable rate
  against the dextramer labels meet a documented floor.
- Stage 4: assert the scramble calibrated verdict returns `suspect` on a known scramble
  input and `reliable` on a known cognate crystal.
- End to end: a tiny fixture repertoire runs stages 0 to 2 and 4 to 5 offline (fold
  mocked), producing a report.

## Open questions for the plan

- Exact TCRdist confidence tier thresholds (to be calibrated on the validation arm).
- Class I MHC allele handling when the dextramer or annotation HLA is known vs inferred.
- Report format (Markdown, HTML, or a claude.ai Artifact for the demo).
