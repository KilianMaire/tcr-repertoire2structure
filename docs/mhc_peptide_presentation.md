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

AUROC for separating a genuine binder from its composition-scramble
(`docs/mhc_peptide_presentation.png`):

| metric | what it measures | A*02:01 AUROC | A*11:01 AUROC |
| --- | --- | --- | --- |
| neg_gpde_groove | MHC-peptide interface PAE (negated) | 0.77 [0.72, 0.82] | 0.98 [0.95, 1.00] |
| iptm_b2m_pep | b2m-peptide ipTM | 0.74 [0.66, 0.81] | 0.97 [0.93, 1.00] |
| iptm_groove | MHC-peptide interface ipTM | 0.71 [0.64, 0.78] | 0.99 [0.97, 1.00] |
| pep_plddt | peptide chain pLDDT | 0.69 [0.62, 0.76] | 0.88 [0.79, 0.96] |
| pep_ptm | peptide chain pTM | 0.69 [0.64, 0.75] | 0.88 [0.79, 0.95] |
| pep_iptm | peptide interface pTM (whole complex) | 0.50 [0.46, 0.54] | 0.62 [0.55, 0.71] |
| ranking_score | Protenix global ranking | 0.46 [0.42, 0.50] | 0.51 [0.43, 0.60] |

## What the numbers say

1. The confidence does judge presentation. The groove interface readouts
   (`iptm_groove`, `neg_gpde_groove`, `iptm_b2m_pep`) separate a real ligand from
   its scramble, cleanly on A*11:01 (AUROC 0.97 to 0.99) and moderately on
   A*02:01 (0.71 to 0.77). The decoys beat the scramble too, which confirms this
   is a binding signal and not a cognate or TCR signal.

2. The cleanest single metric is the groove interface confidence: `iptm_groove`
   (the MHC-peptide ipTM) or `neg_gpde_groove` (its PAE). Peptide chain pLDDT and
   pTM work but are weaker. Two metrics are useless here: `pep_iptm` and
   `ranking_score` sit at chance, so a good global ranking score does not imply a
   presented peptide.

3. The HLA gap is biologically sensible. A*02:01 is a permissive allele with
   loose anchors (P2 and the C-terminus tolerate several residues), so a shuffle
   of an A*02:01 nonamer often remains a passable binder and the separation
   shrinks. A*11:01 has a strict C-terminal lysine or arginine anchor that a
   shuffle usually destroys, giving near perfect separation. This is offered as
   interpretation, not proof; confirming it would need an explicit anchor
   retention check on the scrambles.

## Conclusion

Taken with the TCR result, the picture is surgical rather than merely negative.
A Protenix or AlphaFold confidence readout is informative about MHC-peptide
presentation (does this peptide bind this HLA) and uninformative about TCR-peptide
recognition (does this TCR read this peptide). The groove interface confidence is
a usable presentation filter, strongest on anchor-strict alleles. It is not, and
the variance analysis shows it cannot be, a TCR-specificity oracle. This is the
evidence base for keeping the two honesty rules separate in the output schema:
presentation may be scored, recognition may not.
