# Repertoire2Structure

A multi agent pipeline, orchestrated by Claude, that turns a raw 10x single cell TCR repertoire into QC'd predicted TCR pMHC structures for its top clonotypes, with honest specificity annotation and skeptical structure QC.

Built for the Built with Claude: Life Sciences hackathon (Researcher track).

![architecture](docs/architecture.png)

## What it does

A researcher hands the pipeline a 10x contig CSV. Claude agents return a report that links each top clonotype to a candidate epitope specificity (with a confidence tier), a predicted TCR pMHC structure, and a skeptical QC verdict on whether that structure is trustworthy or a likely geometry hallucination.

The chain: ingest and clonotype curation, honest specificity annotation (TCRdist against labeled references), fold prep and MSA, structure folding (Protenix on Colab, driven through the Playwright MCP), skeptical QC, and a self contained HTML report.

## Two honesty rules, enforced in the output schema

1. Specificity is annotation by similarity, never prediction. A clonotype is annotated only when a TCRdist neighbor is close enough, always with the distance and a confidence tier. Clonotypes with no close neighbor are flagged unannotatable. No label is ever forced.
2. A predicted structure does not confirm specificity. Protenix imposes canonical TCR pMHC docking geometry even on non binding sequences, so the QC step is a skeptical judge that flags a fold as suspect when its CDR3 to peptide contact does not beat the scramble control calibration.

## Layers

A deterministic stage layer (pure Python, fully tested offline) carries reliability. On top, a genuine multi agent layer built on the Claude Agent SDK exposes the stages as in process tools and delegates from an orchestrator to specialist agents (a fold agent that drives the browser, a skeptical QC agent, a report agent).

## Datasets

- Validation arm (ground truth): a 10x 4 donor CD8 dextramer set, used to measure precision, recall, and unannotatable rate of the annotation step against the dextramer label.
- Application arm (scale): TABLO (Zenodo 10.5281/zenodo.13119615), a large unlabeled human repertoire, run end to end.

## Development

```
python3.11 -m venv .venv
./.venv/bin/pip install -e . -e ~/imgt-api
./.venv/bin/python -m pytest -q
```

Design and plan live in `docs/superpowers/`. The live integration checklist (real datasets, real V domain reconstruction, real folds) is documented at the end of the plan.
