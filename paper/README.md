# Paper package

Reproducible results, figures, tables, and text for the Repertoire2Structure
paper. Self-contained: everything here regenerates from committed data without a
GPU.

## Contents

- `outline.md` : the paper plan (claims, figure panels, supplementary, tables).
- `methods.md` : Materials and Methods, including the pre-registration and
  anti-HARKing methodology section.
- `sources.md` : datasets, tools, accessions, and the reproducibility note.
- `data/` : derived tidy CSVs that every figure and table reads. Committed.
- `tables/` : Table 1 (datasets) and Table 2 (pre-registration outcomes).
- `figures/` : generated main figures and a manifest mapping each figure to its
  script and data.
- `supplementary/` : supplementary figures and tables.
- `make_paper_data.py` : regenerates `data/*.csv` from the raw folds in `runs/`.

## Reproduce

Derived data and figures (no GPU):

```
python paper/make_paper_data.py                    # data/*.csv from runs/
python scripts/plot_confidence_variance.py         # Fig 4
python scripts/plot_mhc_scramble.py                # Fig 5
```

Validation arm from scratch (CPU, needs the dextramer dataset under data/dataset):

```
python scripts/run_validation_arm.py data/dataset out.json 1
```

Raw folds (need a GPU, not committed) are regenerated through the Colab notebooks
the pipeline emits; see `docs/fold_procedure.md`.

## Data provenance

Raw Protenix folds live in `runs/` (gitignored: GPU-only, and TABLO is not
redistributable). They are reduced to the committed `data/*.csv` by
`make_paper_data.py`. The validation arm is fully committed in
`docs/validation_donor1_metrics.json`. This two layer split (raw local, derived
committed) is what makes the figures reproducible from a clean clone.

## Status of the findings

- Honest annotation and leakage guard: reproduced from committed metrics.
- Structure vs sequence (TCR recognition): pre-registered confirmation not
  licensed (primary p = 0.34; reconstructed-only reanalysis p = 0.09). Reported
  as a bounded negative.
- MHC-peptide presentation: groove confidence separates a ligand from its
  scramble, AUROC up to 0.99; robust to stub exclusion.
- Data-quality caveat: 25 of 72 folded clonotypes used a poly-G TCR stub; the
  TCR-recognition analysis excludes them (`data/stub_contamination.csv`).
