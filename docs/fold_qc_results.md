# Fold arm results (real Protenix structures + skeptical QC)

Three TABLO public CD8 clonotypes were folded on Protenix
(`protenix_base_default_v1.0.0`) on an A100, each as a cognate construct and a
scramble control (peptide shuffled, everything else identical). QC scores the
CDR3(Vbeta) to peptide heavy-atom contact and calibrates it against the scramble
null. Contact is averaged across the Protenix sample ensemble, since the docking
pose varies substantially between samples of one fold.

| clonotype | epitope (HLA) | cognate contact | scramble contact | verdict |
| --- | --- | --- | --- | --- |
| 63ed9ff4e42b | KLGGALQAK, CMV IE1 (A\*03:01) | 68 | 29 | reliable |
| 98b9ddcabb19 | NQKLIANQF, EBV (B\*15:01) | 48 | 14 | reliable |
| d952b775645a | CTELKLSDY, EBV (A\*01:01) | 145 | 109 | reliable (weak) |

Mean Vbeta to peptide heavy-atom contacts across 5 samples per construct.

## The scramble control is the point

The scramble contacts are not zero. They are 29, 14 and 109. That is the
structural hallucination this QC exists to catch: Protenix docks the TCR onto
the peptide and forms a plausible interface even when the peptide sequence is a
shuffle. A clean-looking fold therefore does not confirm specificity.

The honest signal is the margin between cognate and its own scramble. All three
cognates exceed their scramble, so they pass, but d952b775645a is marginal: its
scramble reaches 109 against a cognate 145, so the structure is called reliable
only weakly and should not be read as evidence that CTELKLSDY is the true
epitope. This is exactly the skepticism Honesty Rule 2 encodes: the verdict is
about beating the scramble null, never about confirming the annotation.

Note on this run: folded MSA-free (single sequence) after the MSA server
throttled and wedged an earlier MSA-based batch. MSA-free models are lower
confidence, which makes the pose spread across samples wider and the QC margin
the honest thing to report.

Reproduce: build inputs with `scripts/build_protenix_inputs.py`, fold with the
notebook from `scripts/build_colab_notebook.py`, then
`python scripts/qc_folds.py <folds_dir> <run_dir> report.html`.

## MSA pre-fold validation (flu M1, 2026-07-09)

The v1 MSA slice added a Colab CPU cell that computes one unpaired a3m per
unique protein chain via the ColabFold MMseqs2 API and injects it as
`unpairedMsaPath`, before the GPU fold. Validated live on an A100 with the flu
M1 clonotype 9ab6b3bfa998 (epitope GILGFVFTL, HLA-A\*02:01, TRBV19, a known
cognate pair, tcrdist 0 to a reference).

The cognate was folded in three configurations, five samples each:

| configuration | CDR3(Vbeta) to peptide contacts (5 samples) | median | pLDDT | iPTM |
| --- | --- | --- | --- | --- |
| MSA consumed (no use-msa-false flag) | 39, 42, 42, 38, 36 | 39.0 | 95.3 | 0.915 |
| MSA in JSON but use-msa-false set | 17, 0, 0, 0, 0 | 0.0 | 46.2 | 0.17 |
| MSA free (single sequence) | 13, 0, 0, 0, 0 | 0.0 | 46.2 | 0.17 |

Two findings, both load bearing.

First, the `--use_msa false` flag fully suppresses a provided `unpairedMsaPath`.
The middle row (MSA present in the input JSON, flag still set) is identical to
the MSA free row: pLDDT 46, iPTM 0.17, four of five poses undocked. So the flag
had to be dropped for Protenix to consume the precomputed MSA. The write inputs
cell prints `MSA_IN_INPUT chains: 8` as positive proof the paths reached the
folded JSON, so a wrong cell run order can never silently fold MSA free again.

Second, once the MSA is actually consumed the fold is transformed: pLDDT 46 to
95, iPTM 0.17 to 0.915, and the TCR docks in all five poses (median contact 0 to
39) rather than one of five. This is the docking the MSA free run could not
produce.

Cost note: with the full MSA depth (about 150 sequences per chain) the fold
slowed from about 4 minutes to tens of minutes per construct on the same A100.
Capping the MSA depth (for example 32 to 64 sequences) is the v2 lever to make
it practical across many clonotypes. The scramble control with MSA was not
folded in this session (interrupted for the confidence check), so the calibrated
cognate versus scramble margin under MSA is still open.
