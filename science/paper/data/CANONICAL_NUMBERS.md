# Canonical numbers (single source of truth for all figures and text)

Every figure and every sentence must use these values. They are computed from the
committed CSVs in `paper/data/`. Convention: **reconstructed TCRs only** (poly-G
stub folds excluded) and **exact one-sided binomial** tests against per-panel chance.
Do not use the stale pre-stub numbers (76 percent between-TCR, 2.3 percent cognate)
that appear in older drafts.

## Annotation validation (Fig 2) — source: validation_*.csv, docs/validation_donor1_metrics.json
- Labeled clonotypes: 3,325 (donor 1).
- Raw: precision 0.89, recall 0.45, unannotatable 0.50 (1,665 predicted, 1,485 correct).
- De-leaked: precision 0.78, recall 0.08, unannotatable 0.90 (193 predicted, 151 correct).
- Leakage: 1,472 of 3,325 (44 percent) within TCRdist 1 of a reference; median correct-call distance 0.
- Threshold sweep (de-leaked precision): 1.00 at cut 6, ~0.79 to 0.81 through cut 24,
  0.78 at 48, 0.74 at 60, 0.58 at 90.

## Retrieval (Fig 3) — source: tcr_retrieval_top1.csv, panel="reconstructed"
Discovery A*02:01, n=29, chance 0.25, sequence baseline 0.00 (all novel):
- iptm_TCRpep_mean 0.66 (19/29, p=4.9e-6)
- iptm_alpha_pep 0.62 (18/29, p=2.6e-5)
- iptm_beta_pep 0.59 (17/29, p=1.2e-4); iptm_TCRpep_max 0.59 (17/29, p=1.2e-4)
- neg_gpde_beta_pep / iptm_global / ptm_global / ranking_score 0.55 (16/29, p=5.0e-4)
- iptm_groove_ctrl (negative control) 0.24 (7/29, p=0.61, below chance)
- CDR3b-peptide contact: refuted, 0.19, below chance (source: docs/benchmark_preregistration.md)

Held-out A*11:01, n=18, chance 0.5 (pre-registered primary = iptm_TCRpep_max):
- iptm_TCRpep_max 0.61 (11/18), exact binomial p=0.24 -> NOT licensed (bar was p<0.05)
- Robustness (all fail): frozen all-folds n=24 gives 0.583 (14/24, p=0.27);
  reconstructed permutation p=0.09. Every variant fails.

## Variance decomposition (Fig 4) — source: confidence_variance.csv (reconstructed)
Discovery A*02:01 iptm_TCRpep_max: TCR identity 51.7%, peptide identity 23.4%,
  cognate status 7.4%, residual 17.6%; ICC 0.367; cognate delta +0.094 ipTM
  [95% CI 0.051, 0.135] -> significant.
Held-out A*11:01 iptm_TCRpep_max: TCR 66.2%, peptide 7.2%, cognate 3.5%, residual ~23%;
  ICC 0.350; cognate delta +0.030 [-0.011, 0.072] -> ns, does not replicate.
Groove control A*02:01: TCR 8.7%, peptide 83.8%, cognate 1.2% (mirror image).
Groove control A*11:01: TCR 13.6%, peptide 68.0%, cognate 10.1%.

## Presentation (Fig 5) — source: mhc_presentation.csv (reconstructed)
Binder-vs-scramble AUROC (A*02:01 n=29 / A*11:01 n=18):
- neg_gpde_groove 0.82 / 0.995
- iptm_b2m_pep 0.79 / 0.99
- iptm_groove 0.77 / 0.99
- pep_ptm 0.75 / 0.88
- pep_plddt 0.72 / 0.88
- pep_iptm 0.52 / 0.64
- ranking_score 0.47 / 0.54 (at chance)

Per-residue pLDDT for one A*11:01 clonotype (source: peptide_plddt.csv):
- cognate AVFDRKSDAK: 75 to 97 (high)
- scramble (same 10 residues shuffled): 59 to 65 (low)

Anchor permissiveness (source: scramble_anchor_permissiveness.csv):
- A*02:01: cognate anchor frac 0.958, scramble 0.604; groove median 0.976 vs 0.959
- A*11:01: cognate 1.0, scramble 0.5; groove median 0.972 vs 0.940
- Reading: the allele gap is scoring-level permissiveness, not anchor loss.

## Data quality (stub contamination) — source: stub_contamination.csv
- 25 of 72 folded clonotypes used a poly-G TCR stub (panel1 19/48, a1101 6/24).
- All TCR-recognition analyses exclude them: discovery n drops 48 -> 29, held-out 24 -> 18.

## Structures (renders in paper/figures/, produced by scripts/render_structure_pymol.py)
- _struct_view1.png : whole complex, TCR-up canonical side view (flu GILGFVFTL, A*02:01)
- _struct_view2.png : groove looking down the cleft (same complex)
- _interface_flu.png : recognition interface, TCR CDR loops (blue Va, green Vb) contacting peptide (orange)
- _groove_conf_cognate.png / _groove_conf_scramble.png : A*11:01 groove, peptide colored by per-residue pLDDT, residues labelled

## Figure conventions (enforced by scripts/figstyle.py)
- No figure titles. No result text baked into the figure (it goes in the caption).
- Bold "(a)", "(b)" panel tags only, via figstyle.panel_label.
- Okabe-Ito palette from figstyle.PALETTE; SURFACE background; despine.
- Save with figstyle.save(fig, stem) -> writes both PDF (manuscript) and PNG.
- No em-dashes or hyphen-as-punctuation in any caption text (house style).
