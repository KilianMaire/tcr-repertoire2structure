# Structural confidence judges MHC-peptide presentation, not TCR recognition

This is the mirror image of the TCR result in `confidence_variance_analysis.md`,
and together the two draw the honest boundary of what a Protenix or AlphaFold
confidence readout can and cannot say about a TCR-pMHC fold.

## Question

The TCR side showed that the interface confidence separates TCRs, not epitopes.
The natural follow up: does the same confidence carry the *other* axis of
specificity, whether the peptide is a genuine ligand for the HLA (presentation)?
If it does, the model is informative about pMHC binding while being blind to TCR
recognition, which is a clean and mechanistically sensible split.

## Method

Read only over the folded panels (no refold). The scramble control is a
composition preserving shuffle of the cognate: it holds amino-acid content and
length fixed and only breaks the anchor order, so it is a matched non-binder. The
cognate and the same-HLA decoys are all genuine binders. A presentation metric
should therefore separate binders (cognate + decoys) from the scramble,
independent of TCR recognition.

Per candidate metric we take the per (TCR, epitope) median over the 5 Protenix
samples, then report the within-TCR paired fraction (cognate > scramble, and
decoy > scramble) and the pooled binder-vs-scramble AUROC with a TCR-bootstrap CI.

Reproduce: `python scripts/analyze_mhc_scramble.py runs/panel1 runs/hla_a1101`
and `python scripts/plot_mhc_scramble.py`.

## Results

AUROC for separating a genuine binder from its composition-scramble, reconstructed
TCRs only (`docs/mhc_peptide_presentation.png`; groove metrics do not touch the TCR
chains, so excluding the poly-G stubs barely moves them):

| metric | what it measures | A*02:01 AUROC | A*11:01 AUROC |
| --- | --- | --- | --- |
| iptm_groove | MHC-peptide interface ipTM | 0.77 [0.69, 0.86] | 0.99 [0.97, 1.00] |
| neg_gpde_groove | MHC-peptide interface PAE (negated) | 0.82 [0.74, 0.89] | 1.00 [0.98, 1.00] |
| iptm_b2m_pep | b2m-peptide ipTM | 0.79 [0.70, 0.88] | 0.99 [0.96, 1.00] |
| pep_plddt | peptide chain pLDDT | 0.72 [0.64, 0.80] | 0.88 [0.78, 0.98] |
| pep_ptm | peptide chain pTM | 0.75 [0.67, 0.83] | 0.88 [0.80, 0.94] |
| pep_iptm | peptide interface pTM (whole complex) | 0.52 [0.47, 0.58] | 0.64 [0.52, 0.75] |
| ranking_score | Protenix global ranking | 0.47 [0.43, 0.51] | 0.54 [0.41, 0.68] |

## What the numbers say

1. The confidence does judge presentation. The groove interface readouts
   (`iptm_groove`, `neg_gpde_groove`, `iptm_b2m_pep`) separate a real ligand from
   its scramble, cleanly on A*11:01 (AUROC 0.97 to 0.99) and moderately on
   A*02:01 (0.71 to 0.77). The decoys beat the scramble too, which confirms this
   is a binding signal and not a cognate or TCR signal.

2. The cleanest single metric is the groove interface confidence: `iptm_groove`
   (the MHC-peptide ipTM) or `neg_gpde_groove` (its PAE). Peptide chain pLDDT and
   pTM work but are weaker. `ranking_score` sits at chance (a good global ranking
   score does not imply a presented peptide), and `pep_iptm` is at chance on
   A*02:01 (0.52) and only weakly informative on A*11:01 (0.64, CI touching 0.5),
   so neither is a reliable presentation filter.

3. The HLA gap is permissiveness at the scoring level, not anchor loss in the
   scramble. We first guessed the scramble simply destroys the anchor more often
   on A*11:01, and tested it: the fraction of scrambles that still satisfy the
   allele anchor motif is similar for the two alleles (60% for A*02:01, 50% for
   A*11:01), so anchor destruction does not account for the AUROC gap. The real
   mechanism is in the absolute groove confidence. A*02:01 is a permissive
   groove: Protenix seats an A*02:01 scramble almost as confidently as a genuine
   ligand (median iptm_groove 0.959 vs 0.976, a 0.018 gap with heavily
   overlapping spread), which compresses the separation. On A*11:01 the scramble
   confidence drops and spreads (0.940 vs 0.972, scramble lower quartile 0.88),
   giving clean separation. The readout inherits the allele's binding
   permissiveness: it separates ligand from non-ligand well only where the groove
   itself is discriminating.
   Reproduce: `python scripts/check_scramble_anchors.py runs/panel1 runs/hla_a1101`

## Conclusion

Taken with the TCR result, the picture is surgical rather than merely negative.
A Protenix or AlphaFold confidence readout is informative about MHC-peptide
presentation (does this peptide bind this HLA) and uninformative about TCR-peptide
recognition (does this TCR read this peptide). The groove interface confidence is
a usable presentation filter, strongest on discriminating (less permissive)
alleles. It is not, and
the variance analysis shows it cannot be, a TCR-specificity oracle. This is the
evidence base for keeping the two honesty rules separate in the output schema:
presentation may be scored, recognition may not.
