# tcrdock live validation checklist

Status: builder reworked to the gene-level TSV contract (offline, merged). Steps 2
to 4 need a Colab GPU session and the AF2 params, done together with the user on a
real class I 10x CSV.

## Step 1 (DONE, offline)

`tcrdock_inputs.build(clonotype, annotation)` emits cognate + scramble rows for the
10 column tcrdock TSV: organism, mhc_class, mhc, peptide, va, ja, cdr3a, vb, jb,
cdr3b. Guards missing peptide/HLA, enforces the column contract. 95/95 green.

## Validation seed data

`data/validation_tcrdock_classI.csv` (8 clonotypes, 16 contig rows): real paired
class I TCRs pulled from tcr_explorer's records index (human MHCI), 2 clonotypes each
for GILGFVFTL (flu M1, A*02:01), NLVPMVATV (CMV pp65, A*02:01), GLCTLVAML (EBV BMLF1,
A*02:01) and KLGGALQAK (A*03:01). Because they come from the index `annotate()`
searches, they self-match at tcrdist 0 -> high tier, so all 8 build valid tcrdock
rows offline. Annotation is circular (irrelevant to fold/QC validation). The local
dextramer fixture is 2 synthetic clonotypes that do NOT annotate (dead end); TABLO is
3.1M mostly-unannotatable rows.

## Step 2: real Colab cell (`tcrdock_notebook.py`)

Replace the `raise NotImplementedError` scaffold with the real two step invocation:

1. `git clone https://github.com/phbradley/TCRdock`, install its bundled AF2 fork.
2. Write the embedded rows to `targets.tsv`.
3. `python setup_for_alphafold.py --targets_tsvfile targets.tsv --output_dir setup`.
4. `python run_prediction.py --targets setup/targets.tsv --outfile_prefix run
   --model_names model_2_ptm --data_dir $ALPHAFOLD_DATA_DIR`.
5. Emit the CIF with canonical chain IDs A=TCRa, B=TCRb, C=MHC heavy, D=b2m,
   E=peptide (tcrdock's own chain labels must be remapped).

## Open unknowns to resolve on the first real run (do NOT guess-fix in the builder)

Raised by adversarial review of the builder; each needs a real tcrdock run to settle:

1. `mhc_class` is emitted as int `1`. Confirm tcrdock's parser wants `"1"` and not
   the roman `"I"` or an enum. Fix the builder only if the real TSV parser rejects 1.
2. J genes (`ja`/`jb`) are emitted BARE (e.g. `TRAJ33`); there is no `traj_allele`
   field in `Clonotype`. If tcrdock's gene table is keyed by `gene*allele` and does
   not tolerate a bare J gene (or silently assumes `*01`), that is a data model gap:
   add `traj_allele`/`trbj_allele` upstream. Check against tcrdock's genes TSV.
3. HLA formatting: RESOLVED for the multi-value case. Real annotation output returns
   comma-joined values like `HLA-A*02,HLA-A*02:01`; `_mhc_allele` now strips `HLA-`
   and keeps the most specific (colon-bearing) token, giving a clean `A*02:01`. Still
   to confirm live: 3 field (`A*02:01:01`) and lowercase input against tcrdock's lookup.
4. Scramble degeneracy (pre-existing, in `construct_io.scramble_peptide`): for very
   short or homopolymer like peptides the reverse+rotate can equal the original,
   weakening the calibration null. Backlog; verify the real peptides are not degenerate.

## Live run 1 (2026-07-08, Colab CoquetLab, T4 15GB)

Drive proven end-to-end: login (CoquetLab) -> upload ipynb -> GPU runtime -> run cell
-> read output, all via Playwright. Cell 1 OK (Tesla T4, TCRdock cloned). Cell 2:
`download_blast.py` OK (ncbi-blast 2.11.0), but `pip install -r requirements.txt`
FAILED with `ModuleNotFoundError: No module named 'setuptools.extern.six'`.

Root cause: Colab is Python 3.12 + modern setuptools; tcrdock's requirements target
py3.8 and pin old packages whose setup.py imports `setuptools.extern.six` (gone in
setuptools >= 58). Even past that, tcrdock's bundled AF2/JAX fork is built for py3.8;
its jaxlib will not match py3.12 + CUDA 12. Conclusion: tcrdock needs a real py3.8/3.9
environment on Colab (condacolab), not a one-line patch. The notebook builder is at
`scripts/build_tcrdock_notebook.py`; env bring-up is the open task before a fold runs.

## Step 3: Playwright drive

Validate the executor (`tcrdock-agent`) opens the notebook, runs all cells, waits,
downloads the CIF, and calls `record_fold_result(tool="tcrdock")`. It must report
not-run (never fabricate) if the notebook/adapter is absent.

## Step 4: calibration

Run cognate + scramble on the real CSV, read the CDR3-to-peptide contact for both,
and set the tcrdock scramble_threshold from that null. Distinct from Protenix's; never
shared. QC uses `qc_metric="cdr3_peptide"` for tcrdock.
