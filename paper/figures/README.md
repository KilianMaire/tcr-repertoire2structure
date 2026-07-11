# Figure manifest

Each main figure, its generating script, its data source, and its status.

| figure | message | script | data | status |
| --- | --- | --- | --- | --- |
| Fig 1 | pipeline and the two honesty rules | (to draw) | `docs/architecture.png` | to build |
| Fig 2 | honest annotation and the leakage guard | (to build) | `data/validation_annotation.csv`, `validation_threshold_sweep.csv`, `validation_leakage.csv` | to build |
| Fig 3 | structure vs sequence, discovery and held-out, with the pre-registration verdict | (to build) | `data/tcr_retrieval_top1.csv`, `tables/table2_preregistration.md` | to build |
| Fig 4 | confidence separates TCRs, not epitopes | `scripts/plot_confidence_variance.py` | raw folds via `analyze_confidence_variance.py` | **needs recon-only regen** |
| Fig 5 | confidence reads MHC-peptide presentation | `scripts/plot_mhc_scramble.py` | raw folds via `analyze_mhc_scramble.py` | final (robust to stubs) |
| Fig 6 | the two-axis map: presentation yes, recognition no | (to draw) | conceptual | to build |

## Known correction pending

`fig4_confidence_variance.png` currently shows the FULL-panel numbers (var TCR
76%, cognate 2.3%). An audit found 25 of 72 folded clonotypes used a poly-G TCR
stub; the canonical figure must be regenerated on reconstructed TCRs only, where
the cognate variance is larger (7.4% discovery, 3.5% held-out). See
`data/stub_contamination.csv` and `scripts/analyze_stub_contamination.py`. Fig 5
(presentation, groove-based) is unaffected.

## Significance reporting note

Report retrieval significance as an exact binomial test of Top-1 against naive
per-panel chance (0.25 on the discovery 4-way panel, 0.5 on the held-out binary),
not the label-permutation p or the TCR-blind null in the raw benchmark reports:
an audit showed those two are miscalibrated on the discovery panel (they assign
p near 1e-4 even to the negative control). By the correct binomial test, discovery
`iptm_TCRpep_max` (25/48 vs 0.25) is genuinely significant (p about 5e-5) and the
groove control (9/48) is correctly not (p about 0.88); the held-out primary
(14/24 vs 0.5) is not significant (p about 0.34), consistent with the failed
pre-registration.
