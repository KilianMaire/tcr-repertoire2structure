# Pre-registration: held-out confirmation of the structural-confidence readout

Committed **2026-07-11, before any confirmation fold was run.** Its purpose is
to close the tuning-on-truth loophole in the discovery result: the winning
readout below was chosen post-hoc from a battery of 12, so it must be
re-tested on data that played no part in that choice, with the metric fixed
in advance and no re-selection at scoring time.

## Discovery result (not confirmatory)

Donor1, HLA-A*02:01, n=48 unannotatable (tcrdist-novel) CD8 TCRs. Panel per
TCR: cognate + 3 same-HLA decoys + composition-scramble, 5 Protenix samples.

- Sequence baseline (tcrdist): Top-1 0.000 (all panel TCRs are unannotatable).
- CDR3b-peptide **contact: refuted** (Top-1 0.19, below chance 0.25; within-pair
  cognate>scramble mean_delta -26.9, negative).
- Structural **confidence** recovers the epitope. Best of a 12-readout battery:
  `iptm_TCRpep_max` Top-1 0.521, CI[0.375,0.667], TCR-blind null 0.021,
  label-permutation p=1e-4. CDR3b pLDDT 0.417. Negative control
  `iptm_groove_ctrl` (MHC-peptide interface) 0.188, below chance.

Because `iptm_TCRpep_max` was selected as the best of 12, its 0.521 is an
optimistic point estimate. The test below is what licenses claiming it.

## Confirmation design (frozen here)

- **Held-out set:** donor1, HLA-A*11:01, unannotatable-only, cognates balanced
  across the two available A*11:01 epitopes (AVFDRKSDAK, IVTDFSVIK). Same-HLA
  decoys only, so the panel is a balanced binary retrieval: naive chance 0.5,
  and the TCR-blind (rank-by-peptide-mean) null is also ~0.5. This is a
  different HLA and a different peptide set from the discovery panel; it tests
  whether the readout generalizes beyond A*02:01. It does **not** close the
  single-donor axis (still donor1).
- **Primary metric (pre-committed, single):** `iptm_TCRpep_max`
  = max(chain_pair_iptm[TCRalpha][peptide], chain_pair_iptm[TCRbeta][peptide]),
  median over the 5 samples. No re-selection from the battery.
- **Secondary (pre-committed):** `iptm_TCRpep_mean`, `neg_gpde_beta_pep`.
- **Negative control:** `iptm_groove_ctrl` (MHC-peptide interface ipTM) must
  NOT exceed chance; if it does, the signal is peptide-in-groove, not TCR
  recognition, and the confirmation fails.

## Pre-registered predictions

1. **Primary:** `iptm_TCRpep_max` Top-1 exceeds both chance (0.5) and the
   TCR-blind null, with label-permutation p < 0.05.
2. **Contact** stays at or below chance (the refutation replicates).
3. **Negative control** `iptm_groove_ctrl` does not exceed chance.

Scoring uses the frozen `render_report` in `scripts/run_benchmark_arm.py`; it
prints every readout, but the confirmatory claim reads only the three
pre-committed metrics and the control named above. A 2-of-3 (primary +
control) pass is the bar; secondaries are supporting, not gating.
