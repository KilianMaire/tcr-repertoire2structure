# Supplementary material

All values follow `paper/data/CANONICAL_NUMBERS.md` (reconstructed TCRs only,
exact binomial tests). House style: no em-dashes, no hyphen as punctuation.

## Supplementary figures

### Figure S1. Full retrieval battery, all folds versus reconstructed only.
`figS1_retrieval_strata.pdf`. For each HLA panel (A*02:01 discovery, chance 0.25;
A*11:01 held-out, chance 0.5), every confidence readout is shown twice, once on all
folds and once on the reconstructed-only subset. Excluding the poly-G stub folds
lifts most readouts (for example the discovery primary iptm_TCRpep_max rises from
0.52 on all 48 folds to 0.59 on 29 reconstructed), which is why every
TCR-recognition analysis in the main text uses the reconstructed subset. The
held-out primary still fails to confirm on either stratum. Source:
`paper/data/tcr_retrieval_top1.csv`.

### Figure S2. Per-residue confidence in the groove, cognate versus scramble.
`figS2_groove_confidence.pdf`. One A*11:01 clonotype, its cognate peptide
AVFDRKSDAK and a composition-scramble of the same ten residues, folded identically.
The peptide sticks are coloured by per-residue pLDDT on a shared scale. The cognate
is placed at high confidence (pLDDT 75 to 97) in an extended groove conformation,
the scramble at low confidence (59 to 65) in a non-canonical path. The same amino
acids are placed confidently when they are the true epitope and uncertainly when
shuffled, which is the presentation signal the groove metrics read. Source:
`paper/data/peptide_plddt.csv`; renders from `scripts/render_structure_pymol.py`.

### Figure S3. Chain-pair interface confidence for a confident cognate complex.
`figS3_chain_pair_iptm.pdf`. The mean five-sample chain-pair ipTM matrix for one
A*11:01 cognate complex. The two TCR-to-peptide cells (the basis of iptm_TCRpep)
and the MHC-to-peptide cell (iptm_groove) are outlined. The groove cell is high
(0.96) while the TCR-to-peptide cells are moderate (0.43), which is the
presentation-strong, recognition-weak pattern in matrix form. Source:
`paper/data/chain_pair_iptm_example.csv`.

### Figure S4. Per-sample reproducibility of the confidence readouts.
`figS4_reproducibility.pdf`. The five Protenix samples of one A*11:01 clonotype,
cognate versus its composition-scramble, for three readouts. iptm_groove separates
the two by a wide, stable margin across all five samples (about 0.96 versus 0.77),
iptm_TCRpep separates them by a small margin, and ranking_score does not separate
them. The retrieval signal is not sample noise. Source:
`paper/data/per_sample_readouts.csv`.

### Figure S5. Gallery of predicted TCR-pMHC complexes across epitopes.
`figS5_complex_gallery.pdf`. Four predicted complexes in the same canonical TCR-up
orientation, three A*02:01 epitopes (GILGFVFTL, ELAGIGILTV, FLYALALLL) and one
A*11:01 epitope (IVTDFSVIK), coloured by chain. Each panel carries an inset that
zooms on the peptide (orange sticks) in the MHC groove, shown as a recessive grey
cartoon. The pipeline produces consistent well-folded complexes across epitopes and
both alleles, not a single cherry-picked fold. Renders and per-epitope selection
from `scripts/render_gallery.py`.

## Supplementary tables

### Table S1. Structure-tool registry and validity domains.
Source: `src/rep2struct/structure_tools.py`. The pipeline selects a tool by validity
and honours each tool's scope; only the default folds the panels in this paper.

| tool | output | MHC class | needs TCR | honest scope |
| --- | --- | --- | --- | --- |
| protenix (default) | structure | I and II | either | full three-chain TCR-pMHC fold; imposes canonical geometry even on non-binders, which is the basis of the scramble QC |
| af3 | structure | I and II | either | AlphaFold3-class accuracy when weights are available; gated, only if the user supplies weights |
| mhcfine | structure | I | no | most precise class I peptide pose (0.54A backbone RMSD vs 6VRN in live calibration); placement only, does not separate binder from scramble, never evidence of recognition |
| tcrdock | structure | I and II | yes | TCR:pMHC interface and V-domain anchoring; recognition judged by interface PAE against the tool's own folded scramble null |
| affinetune | binding score | I and II | no | is-this-peptide-presented classifier; returns a presentation score, not a structure |

### Table S2. MSA effect on one flu M1 fold (A*11:01 example clonotype 9ab6b3bfa998, five samples each).
Source: `docs/fold_qc_results.md`. Consuming a precomputed MSA transforms the fold;
the `use_msa false` flag silently suppresses a provided MSA.

| configuration | median CDR3-beta to peptide contacts | pLDDT | ipTM |
| --- | --- | --- | --- |
| MSA consumed (no use-msa-false flag) | 39.0 | 95.3 | 0.915 |
| MSA in JSON but use-msa-false set | 0.0 | 46.2 | 0.17 |
| MSA free (single sequence) | 0.0 | 46.2 | 0.17 |

### Table S3. Scramble-control calibration on the first live TABLO folds.
Source: `docs/fold_qc_results.md`. Three public CD8 clonotypes, each folded cognate
and scramble; QC scores the CDR3-beta to peptide contact and reads the margin over
the scramble null. Scramble contacts are not zero, which is the structural
hallucination the QC exists to catch. All three cognates beat their scramble, one
only weakly.

| clonotype | epitope (HLA) | cognate contact | scramble contact | verdict |
| --- | --- | --- | --- | --- |
| 63ed9ff4e42b | KLGGALQAK, CMV IE1 (A*03:01) | 68 | 29 | reliable |
| 98b9ddcabb19 | NQKLIANQF, EBV (B*15:01) | 48 | 14 | reliable |
| d952b775645a | CTELKLSDY, EBV (A*01:01) | 145 | 109 | reliable, thin margin |

### Table S4. Anchor retention and absolute groove confidence, cognate versus scramble.
Source: `paper/data/scramble_anchor_permissiveness.csv`. The scramble keeps its
predicted anchor residues in most cases and its absolute groove confidence is close
to the binder's, so the cleaner binder-versus-scramble separation on A*11:01 is a
permissiveness effect at the scoring level, not a loss of anchors.

| HLA | cognate anchor fraction | scramble anchor fraction | binder groove median | scramble groove median |
| --- | --- | --- | --- | --- |
| A*02:01 | 0.958 | 0.604 | 0.976 | 0.959 |
| A*11:01 | 1.000 | 0.500 | 0.972 | 0.940 |

### Table S5. Companion engine benchmark on class II TCR-pMHC crystals (DockQ).

An independent benchmark on six solved TCR-pMHC-II crystals (murine I-A^b and I-A^u,
human HLA-DR) scored Protenix against Boltz and AF2-multimer by DockQ, on the same five
chain construct order. Protenix reproduces the peptide to TCR-beta interface, the
contact that governs antigen specificity, at high accuracy and far above AF2-multimer.
Boltz scores marginally higher overall, but its calibration cases include long-deposited
entries almost certainly inside training, so its edge is partly recall rather than
generalisation. This motivates Protenix as the folding engine. The benchmark is a
companion study in a different system (mouse class II, a house dust mite allergen
peptide) and is external to the human class I analyses of this paper.

| engine | DockQ total (mean) | peptide to TCR-beta DockQ (mean) |
| --- | --- | --- |
| Boltz | 0.906 | 0.914 |
| Protenix | 0.875 | 0.887 |
| AF2-multimer (ColabFold) | 0.759 | 0.655 |

## Supplementary methods notes

- **Stub contamination.** 25 of 72 folded clonotypes used a poly-G placeholder in
  place of a reconstructed variable domain (panel1 19 of 48, a1101 6 of 24). These
  carry no TCR-interface signal and are excluded from every TCR-recognition analysis
  (`paper/data/stub_contamination.csv`).
- **Significance testing.** The label-permutation p and the TCR-blind null used in
  early analysis are miscalibrated on the discovery panel (they assign p near 1e-4
  even to the negative control). The correct test is an exact one-sided binomial
  against naive per-panel chance, used throughout the main text.
- **MSA cost.** With full MSA depth (about 150 sequences per chain) a fold slows from
  about 4 minutes to tens of minutes on an A100; capping depth to 32 to 64 is the
  practical lever for many clonotypes (`docs/fold_qc_results.md`).
