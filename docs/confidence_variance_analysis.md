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
across the 5 Protenix samples, the value the benchmark ranks on, then partition
its variance:

1. ICC across TCRs (between TCR variance over total). A value near 1 means the
   readout is essentially a property of the TCR, roughly constant over which
   peptide is shown.
2. Sequential variance decomposition. Remove the TCR mean first, then measure how
   much of the total sum of squares cognate status and peptide identity explain
   on the within TCR residual.
3. Within panel cognate effect. The cognate value minus the panel mean, in ipTM
   units, with a bootstrap CI over TCRs and a within TCR permutation p (the
   pseudo cognate is a random panel member).

Reproduce: `python scripts/analyze_confidence_variance.py runs/panel1 runs/hla_a1101`

## Results

| panel | readout | ICC | var TCR | var cognate | var peptide | cognate minus panel mean | perm p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| panel1 (A*02:01, n=48) | iptm_TCRpep_max | 0.685 | 76.0% | 2.3% | 6.1% | +0.062 [+0.032, +0.092] | 0.0001 |
| panel1 (A*02:01, n=48) | iptm_beta_pep | 0.693 | 76.6% | 1.9% | 6.0% | +0.056 [+0.027, +0.085] | 0.0002 |
| a1101 (A*11:01, n=24) | iptm_TCRpep_max | 0.443 | 71.3% | 0.3% | 0.6% | +0.010 [-0.031, +0.051] | 0.32 |
| a1101 (A*11:01, n=24) | iptm_beta_pep | 0.504 | 74.4% | 1.1% | 0.5% | +0.020 [-0.018, +0.060] | 0.17 |
| control: groove (A*02:01) | iptm_groove_ctrl | -0.06 | 19.9% | 0.2% | 25.8% | +0.002 | 0.28 |
| control: groove (A*11:01) | iptm_groove_ctrl | -0.62 | 18.5% | 0.3% | 63.8% | +0.000 | 0.39 |

## What the numbers say

1. The TCR interface confidence is dominated by the TCR. Across both panels
   about three quarters of its variance (76% and 71%) is between TCR, with an ICC
   of 0.4 to 0.7. Show a given TCR any same HLA peptide and the confidence lands
   in roughly the same place. That is the single HLA docking confound, measured.

2. The peptide specificity component is small and does not replicate. On the
   discovery panel there is a genuine but tiny cognate effect (2.3% of variance,
   +0.062 ipTM, p=0.0001). On the held out A*11:01 panel it collapses to 0.3% of
   variance and +0.010 ipTM with p=0.32. The apparent discovery signal was a two
   percent sliver riding on top of the docking effect, and it did not carry to a
   new HLA and peptide set.

3. The negative control validates the decomposition. The MHC peptide groove ipTM
   shows the mirror image structure: its variance is driven by peptide identity
   (26% and 64%), its ICC across TCRs is at or below zero (it is not a TCR
   property), and its cognate effect is null. So the method attributes variance
   to the correct factor for a readout whose nature is known, which is what
   licenses reading the TCR interface result the same way.

## Conclusion

The structural confidence readout separates TCRs, not epitopes. It is mostly a
measure of how confidently a given TCR docks onto the shared HLA. Any peptide
specificity signal is at most a couple of percent of its variance, is swamped by
the TCR docking term, and does not survive a held out test. It cannot be used as
a TCR specificity oracle. This is the honest boundary of the tool, and it is the
reason Honesty Rule 2 (a predicted structure does not confirm specificity) is
enforced in the output schema rather than left to the reader.
