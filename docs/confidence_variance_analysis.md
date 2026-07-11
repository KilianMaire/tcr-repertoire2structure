# Structural confidence separates TCRs, not epitopes

This is the mechanistic result behind the failed pre-registered confirmation. It
explains why a structural confidence readout appeared to retrieve the cognate
epitope on the discovery panel yet did not survive held-out testing, and it does
so with an effect size decomposition rather than a single binary Top-1, so it is
not vulnerable to the "underpowered null" objection.

## Question

An AlphaFold or Protenix interface confidence (`iptm_TCRpep`, the ipTM between a
TCR chain and the peptide) was proposed as a way to place sequence novel TCRs on
an epitope. Does that number carry information about which peptide is the
cognate, or only about how well the TCR docks onto the shared HLA?

## Method

Read only over the already folded panels (no refold). For each TCR the panel is
its cognate plus same HLA decoys. We take the per (TCR, epitope) median readout
across the 5 Protenix samples, then partition its variance: an ICC across TCRs
(between TCR variance over total; near 1 means the readout is a property of the
TCR, roughly constant over which peptide is shown), then a sequential eta-squared
that removes the TCR mean first and measures how much of the total sum of squares
cognate status and peptide identity explain on the within TCR residual, then the
within panel cognate effect (cognate value minus panel mean, in ipTM units, with a
bootstrap CI over TCRs and a within TCR permutation p where the pseudo cognate is
a random member of that same panel).

**Reconstructed TCRs only.** An audit found that a substantial minority of the
folded clonotypes used a poly-G stub TCR rather than a reconstructed V domain
(19 of 48 on A*02:01, 6 of 24 on A*11:01); a stub is a glycine backbone with a
floating CDR3, so its TCR-peptide readouts are uninformative. This analysis is run
on reconstructed clonotypes only. The stub contamination and its effect are in
`paper/data/stub_contamination.csv`.

Reproduce (reconstructed-only is the plot default; the analysis script exposes
both): `python scripts/analyze_stub_contamination.py runs/panel1 runs/hla_a1101`.

## Results (reconstructed TCRs only)

| panel | readout | ICC | var TCR | var cognate | var peptide | cognate minus panel mean | perm p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| panel1 (A*02:01, n=29) | iptm_TCRpep_max | 0.367 | 51.7% | 7.4% | 23.4% | +0.094 [+0.051, +0.135] | 0.0001 |
| panel1 (A*02:01, n=29) | iptm_beta_pep | 0.351 | 50.4% | 7.7% | 23.2% | +0.093 [+0.052, +0.132] | 0.0001 |
| a1101 (A*11:01, n=18) | iptm_TCRpep_max | 0.350 | 66.3% | 3.5% | 7.2% | +0.030 [-0.011, +0.072] | 0.091 |
| a1101 (A*11:01, n=18) | iptm_beta_pep | 0.349 | 66.2% | 4.2% | 6.7% | +0.035 [-0.008, +0.079] | 0.070 |
| control: groove (A*02:01) | iptm_groove_ctrl | -0.21 | 8.7% | 1.2% | 83.8% | +0.002 [-0.001, +0.005] | 0.14 |
| control: groove (A*11:01) | iptm_groove_ctrl | -0.71 | 13.6% | 10.1% | 68.0% | +0.002 [-0.001, +0.004] | 0.08 |

## What the numbers say

1. The TCR interface confidence is dominated by the TCR. Half to two thirds of its
   variance (52% and 66%) is between TCR. Show a given TCR any same HLA peptide and
   the confidence lands in roughly the same neighbourhood. That is the single HLA
   docking confound, measured.

2. The peptide specificity component is real but weak, and it does not confirm on
   the held out set. On the discovery panel the cognate effect is 7.4% of variance,
   +0.094 ipTM, p=0.0001. On the held out A*11:01 panel it is 3.5% of variance,
   +0.030 ipTM, p=0.09: a positive trend that does not clear significance at this n.
   So the signal is small but not absent (an earlier full-panel analysis that
   included the poly-G stubs put it near zero; the stubs were diluting it). The
   honest statement is that the peptide contributes a few percent of the variance,
   swamped by the TCR docking term, and too weak to confirm on a held out HLA.

3. The negative control validates the decomposition. The MHC peptide groove ipTM
   shows the mirror image structure: its variance is driven by peptide identity
   (84% and 68%), its ICC across TCRs is negative (it is not a TCR property), and
   its cognate effect in ipTM units is null (+0.002, not significant). Its cognate
   variance FRACTION reads a noisy 1% to 10% only because the groove is saturated
   near 0.97 for every peptide, so the total variance denominator is tiny; the
   effect size in ipTM units is the honest readout and it is flat. This is why the
   companion figure leads with the effect size panel, not the variance fractions.

## A note on significance testing

The retrieval significance in this project should be read as an exact binomial
test of Top-1 against naive per-panel chance (0.25 on the 4-way discovery panel,
0.5 on the binary held-out panel), not as the label-permutation p or the TCR-blind
null printed in the raw `runs/*/benchmark_report.md`. An audit showed those two are
miscalibrated on the discovery panel: because a permuted cognate often is not even
a member of that TCR's panel, the label-permutation null almost never reaches
chance, so it assigns p near 1e-4 even to the negative control. By the correct
binomial test, discovery `iptm_TCRpep_max` (17 of 29 reconstructed, vs chance 0.25)
is genuinely significant and the groove control is correctly not. The within panel
cognate-effect permutation used above does not have this defect, because its
pseudo cognate is always drawn from the TCR's own panel.

## Conclusion

The structural confidence readout separates TCRs, not epitopes. It is mostly a
measure of how confidently a given TCR docks onto the shared HLA. A genuine
peptide specificity signal exists but is only a few percent of the variance, is
swamped by the TCR docking term, and does not confirm on a held out HLA. It cannot
be used as a TCR specificity oracle. This is the honest boundary of the tool, and
it is the reason Honesty Rule 2 (a predicted structure does not confirm
specificity) is enforced in the annotation and QC logic rather than left to the
reader.
