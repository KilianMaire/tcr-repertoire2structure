# Paper plan

Working title: **Structural confidence in TCR-pMHC prediction reads presentation, not recognition.**

Alternative: An honest TCR repertoire to structure pipeline shows that Protenix
confidence separates TCRs and ligands, but not the epitopes a TCR recognises.

Venue target: a methods or immunoinformatics venue (Bioinformatics, or a short
report). The result is a bounded, pre-registered negative with a positive
counterpart, so honesty of framing is the point, not a leaderboard number.

## One sentence

We built a multi-agent pipeline that turns a raw 10x TCR repertoire into QC'd
TCR-pMHC structures with honest specificity annotation, then used it to ask a
sharp question and answer it against ground truth: an AlphaFold-family confidence
readout is informative about MHC-peptide presentation but not about TCR-peptide
recognition, and its apparent recognition signal fails a pre-registered held-out
test.

## Claims, in the order the paper makes them

1. **Honest annotation.** Similarity annotation (TCRdist to labelled references)
   is precise where it fires but abstains on most novel TCRs, and a leakage guard
   is required or precision is overstated. Backs: `Fig 2`, `Table 1`,
   `paper/data/validation_*.csv`.

2. **A pre-registered test of structure vs sequence.** On sequence-novel TCRs a
   structural confidence readout retrieves the cognate epitope in discovery
   (Top-1 0.52, p=1e-4) but does not survive the pre-registered held-out
   confirmation (Top-1 0.58, p=0.34). Backs: `Fig 3`, `Table 2`,
   `paper/data/tcr_retrieval_top1.csv`.

3. **Mechanism: it separates TCRs, not epitopes.** About 76% of the readout
   variance is a between-TCR docking property; the peptide-specificity component
   is at most 2.3% and does not replicate. A groove negative control shows the
   mirror-image variance structure, validating the decomposition. Backs: `Fig 4`,
   `paper/data/confidence_variance.csv`.

4. **The positive counterpart: it does read presentation.** Using the
   composition-scramble as a matched non-binder, the groove interface confidence
   separates a genuine ligand from its scramble (AUROC up to 0.99), and the
   allele gap (A*02:01 vs A*11:01) is permissiveness at the scoring level, not
   anchor loss. Backs: `Fig 5`, `paper/data/mhc_presentation.csv`,
   `paper/data/scramble_anchor_permissiveness.csv`.

5. **The tool encodes the boundary.** Two honesty rules are enforced in the
   output schema: presentation may be scored, recognition may not; a predicted
   structure never upgrades a specificity call. Backs: `Fig 1`, `Fig 6`.

## Figure panels (main, 4 to 6)

- **Fig 1. Pipeline and honesty rules.** Schematic: 10x contigs, clonotype
  curation, honest annotation, construct build, Protenix fold, skeptical QC,
  report. Callouts for the two honesty rules. Source: `docs/architecture.png`
  redrawn; conceptual.
- **Fig 2. Honest annotation and the leakage guard.** (a) precision, recall,
  unannotatable rate, raw vs de-leaked; (b) precision and recall across the
  TCRdist cut, raw vs de-leaked; (c) distance distribution of correct calls
  (median 0, the leakage signature). Data: `validation_annotation.csv`,
  `validation_threshold_sweep.csv`, `validation_leakage.csv`.
- **Fig 3. Structure vs sequence, discovery and held-out.** (a) Top-1 per
  readout on the A*02:01 discovery panel vs the sequence baseline (0.0);
  (b) the pre-registered A*11:01 held-out confirmation for the single primary
  metric with its permutation p, showing it lands at chance. Data:
  `tcr_retrieval_top1.csv`, `Table 2`.
- **Fig 4. Confidence separates TCRs, not epitopes.** (a) stacked variance
  (TCR, peptide, cognate, residual) per condition; (b) cognate effect size with
  bootstrap CI, significant on discovery, null on held-out. Built:
  `scripts/plot_confidence_variance.py` -> `docs/confidence_variance.png`.
- **Fig 5. Confidence reads MHC-peptide presentation.** binder vs scramble AUROC
  per metric per HLA, with the permissiveness inset (absolute groove medians).
  Built: `scripts/plot_mhc_scramble.py` -> `docs/mhc_peptide_presentation.png`.
- **Fig 6 (optional). The two-axis map.** A schematic placing each structural
  readout on presentation (informative) vs recognition (blind), summarising the
  paper. Conceptual.

## Supplementary

- **S1.** Full confidence battery Top-1 across strata (all nine readouts),
  discovery and held-out. Data: `tcr_retrieval_top1.csv`.
- **S2.** Scramble control calibration on the first live TABLO folds (cognate vs
  scramble contact, the weak-reliable case). Source: `docs/fold_qc_results.md`.
- **S3.** MSA effect on the fold (pLDDT 46 to 95, iPTM 0.17 to 0.915 once the
  precomputed MSA is consumed; the `use_msa false` pitfall). Source:
  `docs/fold_qc_results.md`.
- **S4.** Scramble anchor retention and absolute groove confidence, the
  permissiveness mechanism behind the A*02:01 / A*11:01 gap. Data:
  `scramble_anchor_permissiveness.csv`.
- **S5.** Tool registry and validity domains (Protenix default, and the
  specialised tools with their honest scope). Source: `src/rep2struct/structure_tools.py`.

## Tables

- **Table 1.** Datasets: validation dextramer set (4 donor, per-cell binarized
  specificity) and the TABLO application repertoire, with accessions.
- **Table 2.** Pre-registration: the three frozen predictions and their held-out
  outcomes, with the pass/fail verdict by the pre-committed bar.

## Narrative arc

Build an honest tool, ask whether structure adds what sequence cannot, get an
exciting discovery number, refuse to trust it, pre-register a held-out test, watch
it fail, then dissect why (it reads presentation, not recognition). The negative
is the contribution, and the pipeline is the instrument that made an honest
negative possible.
