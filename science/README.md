# The science behind R2S

This directory holds the research that validated the [Repertoire2Structure tool](../README.md), kept separate from the tool's own code (`src/`, `tests/`). Nothing here is needed to install or run R2S; it is the manuscript, the analyses, and the scripts that produced them.

## Layout

```
paper/       the manuscript (main.md, methods.md, supplementary/), figures, tables, and derived data
analysis/    the validation and analysis write ups the paper draws on:
               validation_results.md         annotation precision/recall/unannotatable on the dextramer arm
               fold_qc_results.md            skeptical structure QC outcomes
               confidence_variance_analysis.md   fold confidence variance across seeds
               mhc_peptide_presentation.md   the pose only (non discriminating) MHC Fine finding
               benchmark_preregistration.md  the pre registered structure vs sequence retrieval design
scripts/     figure, analysis, and study arm scripts (see below)
data/        committed inputs for the study arms
```

## Study arms

- Validation arm (ground truth): a 10x 4 donor CD8 dextramer set, used to measure precision, recall, and unannotatable rate of the annotation step against the dextramer label. Driver: `scripts/run_validation_arm.py`.
- Benchmark arm: the pre registered structure vs sequence retrieval test (does structural contact recover the epitope for TCRs that sequence similarity cannot annotate). Driver: `scripts/run_benchmark_arm.py`.
- Application arm (scale): TABLO, a large unlabeled human repertoire, run end to end. Driver: `scripts/run_tablo_arm.py`.

## Reproducing figures

The plotting scripts (`scripts/plot_*.py`, `scripts/render_*.py`) read the committed CSVs in `paper/data/` and write into `paper/figures/`, so the figures rebuild without a GPU. Scripts that read fold outputs expect the (git ignored) `runs/` directory at the repository root, produced by the study arms on a GPU.

Run scripts from this `science/` directory, for example:

```
python scripts/plot_validation.py
```

## Relationship to prior work

The finding (structural confidence reads TCR over pMHC geometry, not the chemistry of recognition, so a clean fold is not evidence of specificity) is convergent with McMaster et al. (bioRxiv 2026.07.08.737208) on AlphaFold 3. R2S corroborates it from a different engine (Protenix) on a repertoire, under a pre registered honest negative. See `paper/methods.md` and `paper/sources.md`.
