# Structural confidence in TCR-pMHC prediction reads presentation, not recognition

## Abstract

A T cell receptor (TCR) repertoire records which T cells a person carries, but not
what those T cells see. Turning a repertoire into the structures of its TCR-peptide-MHC
complexes, and reading specificity off those structures, would close that gap. We
built Repertoire2Structure, a multi-agent pipeline that takes a raw 10x single-cell
repertoire to quality-controlled TCR-pMHC structures with an explicitly honest
specificity annotation, and we used it to ask one sharp question against ground
truth: does an AlphaFold-family structural confidence readout tell us which peptide
a TCR recognises? On sequence-novel TCRs the readout retrieves the cognate epitope
well above chance in a discovery panel (Top-1 0.66, exact binomial p below 1e-5),
but the effect does not survive a pre-registered held-out test on a second HLA
allele (Top-1 0.61, p 0.24). A variance decomposition explains why: about half of
the readout is a between-TCR docking property and only 7 percent is cognate-specific,
and that cognate component does not replicate. The same confidence signal is strongly
informative about the other specificity axis, MHC-peptide presentation, where it
separates a genuine ligand from a composition-scramble of itself with AUROC up to
0.99. The negative is the contribution: structural confidence is a presentation
detector wearing the costume of a recognition detector, and a pipeline that states
this boundary in its output is more useful than one that hides it. Every number here
reproduces from committed data without a GPU.

## Introduction

Antigen specificity is the central variable of adaptive immunity and the hardest to
measure at scale. High-throughput single-cell platforms now read paired TCR alpha
and beta chains from thousands of cells at once, so the repertoire (the catalogue of
receptors) is cheap. What the repertoire does not contain is the ligand: the peptide,
presented on a particular MHC allele, that each receptor was selected to see. Sequence
based methods try to fill this gap by matching a query TCR to receptors of known
specificity, but they are precise only where a close labelled neighbour exists and
abstain on the large novel fraction that matters most.

Structure is the appealing alternative. If a model could fold a TCR onto a candidate
peptide-MHC and score the result, specificity would follow from physics rather than
from a lookup table, and novelty would stop being a barrier. AlphaFold-family models,
including the Protenix implementation used here, emit interface confidence readouts
(the interface predicted TM-score, ipTM, and related quantities) that are widely
treated as proxies for whether a modelled interface is real. The open question is
whether those readouts carry recognition information (which peptide a given TCR reads)
or only the easier presentation information (whether a peptide binds an MHC groove at
all), and the two are constantly conflated because a confident-looking complex satisfies
both intuitions at once.

We take the conflation seriously and separate the two axes by construction. We fold
matched panels in which the MHC and the decoy peptides are held fixed and only the
cognate assignment changes, we pre-register a held-out confirmation before running it,
and we build a composition-scramble control that keeps a peptide's amino acid content
while destroying its identity. The pipeline that produces these folds is itself part
of the claim: its annotation abstains rather than guess, and its output schema forbids
a predicted structure from upgrading a specificity call. Figure 1 shows one predicted
complex to orient the reader, with the T cell receptor docked over the peptide held
in the MHC groove.

## Results

### An honest annotation with a mandatory leakage guard

Before any structure is folded, the pipeline annotates each clonotype by TCRdist
similarity to a labelled reference and reports precision and recall against 3,325
dextramer-labelled clonotypes from a public 10x CD8 set (donor 1). The raw numbers
look strong (precision 0.89, recall 0.45), but the reference and the test set share
receptors: 1,472 of 3,325 labelled clonotypes (44 percent) sit within TCRdist 1 of a
reference neighbour, and the median distance of a correct call is 0 (Figure 2a, 2c).
Removing those near-identical matches gives the honest figure, precision 0.78 at
recall 0.08 on genuinely novel receptors (Figure 2a). Precision stays near 0.8 as the
distance cut is relaxed to 48 and only then degrades, reaching 0.58 at 90 (Figure 2b),
which is why the pipeline keeps its confidence tiers tight and leaves 90 percent of
novel clonotypes unannotatable rather than force a label the reference cannot support.
The lesson carried forward is that the sequence route is precise but silent on novelty,
so the novel fraction is exactly where a structural route would have to earn its keep.

### A pre-registered test of structure versus sequence

We assembled a discovery panel of 48 A*02:01-restricted CD8 clonotypes that the
annotation leaves unannotatable, so the sequence baseline retrieves the cognate
epitope at Top-1 0.00 by construction. Each clonotype was folded against its cognate
peptide, three same-HLA decoys, and its composition-scramble, five Protenix seeds
each, and every candidate ranked by a battery of confidence readouts. A data-quality
audit found that 25 of 72 folded clonotypes across both panels had used a poly-G
placeholder in place of a reconstructed variable domain, and all TCR-recognition
analyses below exclude those stub folds (29 reconstructed clonotypes remain in the
discovery panel, 18 in the held-out panel; Figure S1 shows the effect of the exclusion
across readouts).

On the discovery panel the structural confidence readouts retrieve the cognate epitope
well above the naive per-panel chance of 0.25: the best, the mean of the two TCR-to-peptide
interface ipTMs, reaches Top-1 0.66 (19 of 29, exact binomial p below 1e-5), and the
whole family of TCR-to-peptide readouts clears chance with p below 1e-3 (Figure 3a).
The negative control, the MHC-to-peptide groove interface, sits at 0.24 (below chance),
and the CDR3-beta-to-peptide heavy-atom contact is refuted (0.19, below chance), so
the signal is specific to the confidence readouts rather than to folding a plausible
complex.

A discovery number chosen as the best of a battery is an optimistic point estimate, so
we pre-registered a single primary metric (the maximum of the two TCR-to-peptide ipTMs),
a held-out set (a second allele, A*11:01, with a different peptide pair), and the
pass bar, all before folding the held-out panel. The primary did not confirm: Top-1
0.61 (11 of 18), exact binomial p 0.24 against the balanced chance of 0.5, which does
not clear the pre-committed 0.05 (Figure 3b). The verdict is robust to how the panel
is cut. The frozen analysis on all 24 held-out folds gives 0.583 and p 0.27, and the
reconstructed-only permutation test gives p 0.09; every variant fails, so we do not
claim that structural confidence retrieves the epitope. Reporting the exciting
discovery number without this test is precisely the harking the pre-registration
exists to prevent.

### The mechanism: confidence separates TCRs, not epitopes

Why does a signal that is real in discovery fail to generalise? A sequential variance
decomposition of the readout answers it (Figure 4a). On the discovery panel, TCR
identity accounts for 52 percent of the readout variance and generic peptide identity
for 23 percent, while the cognate-status component (the part that would carry
recognition, the confidence going up specifically for the true peptide over the decoys)
is 7 percent, and on the held-out panel that cognate component falls to 3.5 percent.
The cognate effect on ipTM is a significant plus 0.094 in discovery (bootstrap 95 percent
interval 0.051 to 0.135) but an indistinguishable-from-zero plus 0.030 in held-out
(interval minus 0.011 to 0.072; Figure 4b). The groove negative control shows the
mirror image, with generic peptide identity dominating and the between-TCR term small,
which validates the decomposition rather than reflecting an artefact of it. The readout
is largely a per-TCR docking property: some receptors fold confidently against almost
anything, and that between-TCR spread, not per-peptide recognition, is most of what
the confidence measures.

### The positive counterpart: confidence reads presentation

The same readouts that are blind to recognition are strongly informative about the
other specificity axis. Treating the composition-scramble as a matched non-binder that
holds amino acid content fixed, the groove interface confidence separates a genuine
ligand from its scramble with AUROC up to 0.82 on A*02:01 and up to 0.99 on A*11:01
(Figure 5a), while the raw model ranking score stays at chance. The interface the
metric scores is the peptide sitting in the MHC groove (Figure 5b). Per-residue
confidence makes the effect concrete: for one A*11:01 clonotype the cognate peptide
AVFDRKSDAK is placed at pLDDT 75 to 97, whereas a scramble of the same ten residues is
placed at 59 to 65 and in a non-canonical conformation (Figure S2). The model is
confident about where each residue sits when the peptide is a true binder and uncertain
when it is not, which is a presentation signal, not a recognition one. The allele gap
between A*02:01 and A*11:01 is a permissiveness effect at the scoring level rather than
a loss of anchors: the scramble retains its predicted anchor residues in most cases
(Figure S4), and the absolute groove confidence of scramble and binder are close, so
the cleaner separation on A*11:01 reflects its stricter groove, not a different
mechanism.

### The tool encodes the boundary

The pipeline is built so that this boundary cannot be crossed by accident. Presentation
may be scored, because the scramble control gives a matched null for it, and the QC
step reports the margin of a fold over its own scramble rather than the raw confidence
of the fold. Recognition may not be upgraded by a structure: a predicted complex never
overrides the sequence annotation, and a clonotype the annotation left unannotatable
stays unannotatable no matter how confident its fold looks. Figure 6 places the two
routes on the two specificity axes. Sequence abstains on novel receptors, structural
confidence occupies the presentation-strong, recognition-blind corner, and the
recognition-strong corner where a true specificity oracle would sit remains unoccupied.

## Discussion

The result is a bounded, pre-registered negative with a positive counterpart. An
AlphaFold-family confidence readout is a good presentation detector and a poor
recognition detector, and the apparent recognition signal that shows up in a discovery
panel is mostly a between-TCR docking property that does not survive a held-out test.
This matters because the two axes are routinely conflated: a confident TCR-pMHC fold is
read as evidence that the TCR recognises the peptide, when the confidence is largely
reporting that the peptide can be presented and that this particular receptor docks
readily. A method that scores folds without separating the axes will report recognition
it has not measured.

Our claim is deliberately narrow. The held-out test moves to a second HLA allele and a
second peptide set but stays within one donor, so it does not close the single-donor
axis, and it uses one model family at a fixed setting. A stronger recognition signal
may exist for other models, deeper multiple-sequence alignments, or explicit interface
energetics, and nothing here rules that out. What the work does establish is a
measurement discipline: fold matched panels that vary only the cognate assignment,
pre-register the confirmation, keep a composition-scramble as the presentation null,
and decompose the variance before believing a headline number. Under that discipline
the honest reading of current structural confidence is that it presents, and does not
yet recognise.

The instrument is the second contribution. Repertoire2Structure runs the whole path
from a raw repertoire to QC'd structures as cooperating agents, and it was the honesty
of its annotation and the presence of its scramble control that made an honest negative
possible rather than an accidental positive. A pipeline that abstains, that scores
presentation against a matched null, and that refuses to let a structure upgrade a
specificity call is the right shape for this problem, because it fails loudly where a
naive pipeline would fail silently.

## Figures

- **Figure 1.** Pipeline and an annotated predicted TCR-pMHC complex (flu GILGFVFTL on
  HLA-A*02:01). (a) the pipeline, repertoire to QC'd structure; (b) the whole complex
  with chains labelled; (c) looking down onto the peptide groove; (d) the recognition
  interface, TCR CDR loops contacting the peptide; (e) the two specificity axes the
  paper measures.
- **Figure 2.** Annotation validation and the leakage guard. (a) precision, recall and
  unannotatable rate, raw versus de-leaked; (b) precision and recall across the TCRdist
  cut; (c) TCRdist percentiles of correct calls (median 0, the leakage signature);
  (d) the abstention waterfall, labelled to de-leaked scored to predicted to correct.
- **Figure 3.** Structure versus sequence retrieval. (a) discovery A*02:01 confidence
  battery against the sequence baseline (0.00) and chance (0.25), stars mark exact
  binomial p below 0.05; (b) the pre-registered held-out A*11:01 primary against chance
  (0.5), which does not confirm (0.61); (c) the two negative controls below chance;
  (d) the primary metric, discovery versus held-out, with exact binomial p; (e) the
  held-out binomial null with the observed 11 of 18 marked.
- **Figure 4.** Where the confidence variance goes. (a) sequential decomposition into
  TCR identity, peptide identity, cognate status and residual; (b) the cognate effect
  on ipTM with bootstrap 95 percent interval, significant in discovery, null in
  held-out; (c) the ICC, the between-TCR share of variance; (d) the cognate-status
  variance fraction, small everywhere and not replicated held-out. Reconstructed TCRs
  only.
- **Figure 5.** Confidence and MHC-peptide presentation. (a) binder versus scramble
  AUROC per metric and HLA; (b, c) the A*11:01 groove, cognate and scramble peptide
  coloured by per-residue pLDDT; (d) per-residue pLDDT along the peptide; (e) anchor
  retention, cognate versus scramble. Reconstructed TCRs only.
- **Figure 6.** The two-axis map. Structural readouts are strong on presentation and
  blind on recognition; the recognition corner is unoccupied.

See `paper/figures/README.md` for the script and data behind each panel, and
`paper/supplementary/supplementary.md` for Figures S1 to S5 and Tables S1 to S4.
