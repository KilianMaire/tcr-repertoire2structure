# Table 2. Pre-registration: predictions and held-out outcomes

Frozen 2026-07-11 before any confirmation fold (`docs/benchmark_preregistration.md`).
Held-out set: donor 1, HLA-A*11:01, unannotatable only, epitopes balanced so naive
chance and the TCR-blind null are both 0.5.

| # | prediction (pre-committed) | bar | held-out result (full panel, frozen) | verdict |
| --- | --- | --- | --- | --- |
| 1 (primary) | `iptm_TCRpep_max` Top-1 exceeds chance and the TCR-blind null | label-permutation p < 0.05 | Top-1 0.583, p = 0.34 | **FAIL** |
| 2 | contact stays at or below chance (refutation replicates) | Top-1 <= 0.5 | Top-1 0.375, delta -17.6 | pass |
| 3 (control) | `iptm_groove_ctrl` does not exceed chance | Top-1 <= 0.5 | Top-1 0.500, p = 1.0 | pass |

**Overall verdict: the confirmation is not licensed.** The single prediction that
would support the positive claim (the primary) fails. The control and the contact
refutation pass, but they do not license the positive on their own. We therefore
do not claim structural confidence retrieves the epitope. The held-out best of the
battery (`neg_gpde_beta_pep`, 0.667, p = 0.10) is not claimed either, since
re-selecting from the battery is the HARKing the pre-registration prevents.

**Post-hoc data-quality reanalysis (not part of the frozen test).** Excluding the
poly-G stub constructs (Methods) raises the held-out primary to +0.030 ipTM at
p = 0.09, still short of the pre-committed 0.05. The conclusion is unchanged; the
mechanism reading is revised from "no peptide signal" to "a weak peptide signal,
underpowered to confirm." See `paper/data/stub_contamination.csv`.
