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
   (reconstructed Top-1 0.59 to 0.66, exact binomial p below 1e-4) but does not
   survive the pre-registered held-out confirmation (Top-1 0.61, 11/18, p=0.24).
   Backs: `Fig 3`, `Table 2`, `paper/data/tcr_retrieval_top1.csv`.

3. **Mechanism: it separates TCRs, not epitopes.** About 52% of the reconstructed
   readout variance is a between-TCR docking property and 23% is generic peptide
   identity; the cognate-specific component is 7.4% in discovery, falls to 3.5% in
   held-out, and does not replicate. A groove negative control shows the
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

## Figure panels (main, all built as dense multi-panels; see figures/README.md)

- **Fig 1. Pipeline and an annotated TCR-pMHC complex.** (a) the pipeline; (b) the
  whole complex, chains labelled; (c) the groove top-down; (d) the recognition
  interface, TCR CDR loops on the peptide; (e) the two specificity axes.
  `fig1_structure`.
- **Fig 2. Annotation and the leakage guard.** (a) precision, recall, unannotatable
  raw vs de-leaked; (b) precision and recall across the TCRdist cut; (c) distance
  percentiles of correct calls (median 0); (d) abstention waterfall.
  `fig2_validation`.
- **Fig 3. Structure vs sequence.** (a) discovery battery vs sequence 0.00 and
  chance 0.25 with exact-binomial stars; (b) pre-registered held-out primary vs
  chance 0.5 (0.61, not licensed); (c) the two negative controls; (d) discovery vs
  held-out for the primary metric; (e) the held-out binomial null. `fig3_retrieval`.
- **Fig 4. Confidence separates TCRs, not epitopes.** (a) stacked variance; (b)
  cognate effect with bootstrap CI; (c) ICC; (d) cognate-status variance fraction,
  small and not replicated. Reconstructed TCRs only. `fig4_confidence_variance`.
- **Fig 5. Confidence reads MHC-peptide presentation.** (a) binder vs scramble AUROC
  per metric per HLA; (b, c) cognate and scramble groove coloured by per-residue
  pLDDT; (d) per-residue pLDDT along the peptide; (e) anchor retention.
  Reconstructed TCRs only. `fig5_mhc_presentation`.
- **Fig 6. The two-axis map.** Presentation (informative) vs recognition (blind),
  anchored to the Fig 3/5 numbers. `fig6_two_axis_map`.

## Supplementary figures (all built; see supplementary/supplementary.md)

- **S1.** Full retrieval battery, all folds vs reconstructed. `figS1_retrieval_strata`.
- **S2.** Per-residue groove confidence, cognate vs scramble. `figS2_groove_confidence`.
- **S3.** Chain-pair interface confidence matrix. `figS3_chain_pair_iptm`.
- **S4.** Per-sample reproducibility of the readouts. `figS4_reproducibility`.
- **S5.** Gallery of predicted complexes across epitopes. `figS5_complex_gallery`.

Supplementary tables S1 to S4 (tool registry, MSA effect, TABLO scramble QC, anchor
retention) live in `supplementary/supplementary.md`.

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
