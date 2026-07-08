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
