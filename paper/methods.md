# Materials and Methods

## Datasets

**Validation set (ground truth for annotation).** The 10x Genomics 4 donor CD8
pMHC dextramer dataset. Each single cell carries a binarized specificity call
against a panel of pMHC dextramers. Cells were collapsed to clonotypes and each
clonotype took the dextramer epitope its cells agreed on (dominant call). Donor 1
contributed 3,325 clonotypes with a ground-truth label. HLA restriction for each
clonotype was read from the dextramer barcode (for example A0201_GILGFVFTL maps
to HLA-A*02:01). See Table 1 and `paper/sources.md` for accessions.

**Application set (scale).** TABLO, a large unlabelled human repertoire
(Zenodo 10.5281/zenodo.13119615), used to run the full pipeline end to end.

The raw contig CSVs are large and, for TABLO, not redistributable; they are not
committed. `paper/sources.md` lists exactly where to obtain each file.

## Pipeline

The pipeline is two layers. A deterministic stage layer (pure Python, tested
offline) carries reliability: ingest and clonotype curation, honest annotation,
construct build, fold orchestration, skeptical QC, and report. On top, a
multi-agent layer built on the Claude Agent SDK exposes the stages as in-process
tools and delegates from an orchestrator to specialist agents (a fold agent that
drives the browser to run Protenix on Colab, a skeptical QC agent, a report
agent). Structure folding is run on a GPU through Colab and is the one step that
is not fully autonomous; it is human-supervised through the browser.

## Clonotype curation and honest annotation

Contigs were parsed to paired alpha and beta clonotypes and V and J alleles were
standardised. Specificity annotation is by similarity, never prediction: for each
clonotype the paired TCRdist to a labelled reference set is computed (through the
TCR Explorer package), and the nearest neighbour's epitope is attached only when
the distance falls under a tier threshold (high at or below 12, medium at or
below 24, low at or below 48). Clonotypes with no close neighbour are left
unannotatable. This is Honesty Rule 1, and it is enforced in the annotation logic
(`src/rep2struct/annotate.py`), which never attaches an epitope to an
unannotatable clonotype and always records the distance and tier for an annotated
one. It is not enforced by the dataclass in `schema.py`, which is a passive
container.

**Leakage guard.** The dextramer study is itself a reference source for TCRdist,
so many of its exact TCRs sit in the reference. A clonotype whose nearest labelled
neighbour is at TCRdist at most 1 is flagged leakage-suspected and excluded from
the de-leaked metrics, which measure generalisation to genuinely novel TCRs
rather than recall of the reference.

## Construct building

Each fold construct is a five chain FASTA in a fixed order: chain A the TCR alpha
variable domain, chain B the TCR beta variable domain, chain C the MHC class I
heavy chain ectodomain, chain D beta 2 microglobulin, chain E the peptide. The V
domains are reconstructed from V gene, J gene and CDR3 through TCR Explorer's
germline reconstruction. The heavy chain ectodomain is fetched once per allele
from the EBI IPD/IMGT-HLA API and cached; the mature beta 2 microglobulin is an
invariant constant. A clonotype whose V domain cannot be reconstructed falls back
to a clearly marked poly-G stub (ten glycines plus CDR3) for both TCR chains. The
fixed chain order is the index convention every downstream confidence readout
depends on (chain pair index 0,1,2,3,4 for A,B,C,D,E).

**Stub contamination and its handling.** An audit found that a substantial
minority of the folded panels used a poly-G stub rather than a reconstructed V
domain (19 of 48 on A*02:01, 6 of 24 on A*11:01), because germline
reconstruction failed for those clonotypes and the stubbed constructs were still
folded. A stubbed TCR is a glycine backbone with a floating CDR3, so any
TCR-peptide interface readout for it is uninformative. The TCR-recognition
analysis is therefore run on reconstructed clonotypes only (an objective,
outcome-independent filter: chain_b_seq does not start with ten glycines), and
the effect of the contamination is documented in
`paper/data/stub_contamination.csv` and `scripts/analyze_stub_contamination.py`.
The MHC-peptide presentation analysis scores the groove (chains C, D, E) and does
not involve the TCR chains, so it is unaffected by the stubs (AUROC changes by at
most 0.06 when they are excluded).

## Structure prediction

Constructs were folded with Protenix (`protenix_base_default_v1.0.0`) at 5 seeds
per construct on an A100 or H100 GPU. When a precomputed MSA was used, one
unpaired a3m per unique protein chain was computed through the ColabFold MMseqs2
API and injected; the `use_msa false` flag fully suppresses a provided MSA and
was therefore dropped. Consuming the MSA transforms the fold (peptide chain pLDDT
from 46 to 95, complex ipTM from 0.17 to 0.915, and the TCR docking in all five
poses rather than one). Protenix writes one summary confidence JSON per sample
next to the CIFs.

## Structural readouts and skeptical QC

For the retrieval and presentation analyses, each candidate epitope was folded as
its own construct and every readout was summarised as the median over the 5
samples, which is insensitive to a lone degenerate pose. Readouts are named by the
chain pair they score, for example `iptm_TCRpep_max` is the maximum of the TCR
alpha to peptide and TCR beta to peptide interface ipTM, and `iptm_groove` is the
MHC heavy to peptide interface ipTM.

**Scramble control.** For every cognate a composition-preserving shuffle of the
peptide was folded as a matched non-binder (same amino acids, same length, broken
order). This is Honesty Rule 2, enforced in the QC logic (`src/rep2struct/qc.py`):
a fold is called reliable only when its CDR3 to peptide heavy-atom contact beats
the scramble null, so a clean-looking structure never on its own confirms
specificity. The ensemble contact used for the verdict is the MEDIAN across
samples, not the mean; a single scramble sample once produced 591 spurious
contacts and inverted a mean-based verdict, which is why the median is used.

## Statistics

Top-1 retrieval scores the cognate as a hit only when it is the strict unique
maximum of its panel. Confidence intervals are 2000-sample bootstraps over TCRs.
The label-permutation p reshuffles the cognate assignment across TCRs and
recomputes Top-1. A TCR-blind null ranks epitopes by their panel-mean readout,
ignoring TCR identity, and is reported alongside naive chance (the mean inverse
panel size). Variance decomposition partitions the per (TCR, epitope) median
readout: an intraclass correlation across TCRs, then a sequential eta-squared for
TCR identity, then, on the within-TCR residual, for cognate status and peptide
identity. Presentation is scored as the AUROC separating genuine binders (cognate
and same-HLA decoys) from the scramble, with a TCR-bootstrap CI.

Because the retrieval panels are built from unannotatable (sequence-novel) TCRs,
the sequence baseline is 0.000 by construction. Any comparison to sequence is
therefore reported with that caveat and never as a bare claim that structure beats
sequence.

## Pre-registration and avoiding HARKing

The structural-confidence retrieval readout was chosen post hoc as the best of a
battery of twelve on the discovery panel. That is a tuning-on-truth loophole:
selecting the winning statistic after seeing the outcome is HARKing (hypothesising
after the results are known) and inflates the point estimate. Rather than report
the selected 0.52 as if it were a test, we pre-registered a held-out confirmation
before running a single confirmation fold (`docs/benchmark_preregistration.md`,
committed 2026-07-11).

The pre-registration fixed, in advance and in writing: the held-out set (donor 1,
HLA-A*11:01, unannotatable only, epitopes balanced so naive chance and the
TCR-blind null are both 0.5); a single primary metric (`iptm_TCRpep_max`, no
re-selection from the battery); the secondaries; the negative control
(`iptm_groove_ctrl` must not exceed chance); and three predictions with a decision
rule. Scoring used the frozen `render_report` path.

**Outcome, reported whatever it was.** The primary metric reached Top-1 0.583 with
a label-permutation p of 0.34, which does not clear the pre-committed p below 0.05.
The negative control stayed at chance and the contact refutation replicated, but
the single prediction that would license the positive claim failed. We therefore
do not claim structural confidence retrieves the epitope. Choosing a different
readout now (the held-out best was `neg_gpde_beta_pep` at 0.667, p 0.10, still not
significant) would repeat the exact HARKing the pre-registration exists to
prevent, so we do not. The negative is the result, and the mechanism analysis
(confidence separates TCRs, not epitopes) explains it.

The frozen pre-registered analysis ran on the full panel, so the full-panel
number (Top-1 0.583, p 0.34) is the official outcome. A post-hoc data-quality
reanalysis on reconstructed TCRs only (excluding the poly-G stubs above) raises
the held-out cognate effect to +0.030 ipTM at p 0.09, and the discovery cognate
variance from 2.3% to 7.4%: the stubs were diluting a weak but real peptide
signal. This reanalysis is post-hoc and does not license the positive claim (it
still does not clear p below 0.05), but it revises the mechanism reading: the
confidence carries a small genuine peptide signal that is underpowered to
confirm at this n, rather than none at all.

## Reproducibility

Raw folds live in `runs/` and are not committed (they need a GPU to regenerate).
The derived tables every figure and table reads are committed under `paper/data/`
and are regenerated from the raw folds by `paper/make_paper_data.py`. Figures are
produced by the scripts in `scripts/` from those tables. The validation arm is
fully reproducible from the committed `docs/validation_donor1_metrics.json` and
`scripts/run_validation_arm.py`.
