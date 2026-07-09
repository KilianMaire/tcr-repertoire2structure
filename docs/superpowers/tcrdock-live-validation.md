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

## Live run 2 (VALIDATED 2026-07-09, Colab CoquetLab, A100 83GB, notebook tcrdock_bringup.ipynb)

End to end GREEN on real data (CELL1 from data/validation_tcrdock_classI.csv, a flu M1
TRBV19 TCR: va TRAV8-3*01, ja TRAJ42*01, cdr3a CAVGARGGSQGNLIF, vb TRBV19*01, jb
TRBJ2-7*01, cdr3b CASSTRAGVEQYF), cognate GILGFVFTL vs scramble TFVFGLIGL, HLA A*02:01.

### Env recipe (the whole point; TCRdock declares almost none of it)

Same side conda py3.8 + shell-out pattern as affinetune (miniforge, `mamba create -n
tcrdock python=3.8`, clone TCRdock, `cudatoolkit=11.1 cudnn`, LD_LIBRARY_PATH
`/opt/conda/envs/tcrdock/lib:/usr/lib64-nvidia` on every shell-out). Two things the thin
TCRdock requirements.txt (only biopython/numpy/pandas/scipy/matplotlib) does NOT give
you, both surfaced live and fixed:

1. The ENTIRE AlphaFold stack is missing (README defers it to "the AlphaFold README").
   run_prediction dies at `import tensorflow`. TCRdock's bundled AF fork is the 2.3.x
   line (it annotates `jax.Array`, so jax 0.2.x is too old and affinetune's AF-2.0-era
   stack is the WRONG one). Install DeepMind's AlphaFold **v2.3.2** python stack:
   `pip install -r https://raw.githubusercontent.com/google-deepmind/alphafold/v2.3.2/requirements.txt`
   (tensorflow-cpu 2.11, dm-haiku 0.0.9, chex 0.0.7, numpy 1.21.6, biopython 1.79),
   then jaxlib for CUDA 11:
   `pip install jax==0.3.25 jaxlib==0.3.25+cuda11.cudnn805 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html`.
   This set is internally consistent (jax 0.3.25 has `jax.Array`; tf-cpu 2.11 + numpy
   1.21.6 agree; cudnn805 matches cudatoolkit 11.1).
2. Biopython gotcha: TCRdock imports `Bio.SubsMat` (removed in biopython >= 1.80). The
   AF 2.3.2 requirements already pin biopython==1.79, so installing that set fixes it;
   but if you install anything that upgrades biopython, force it back to 1.79.

### Open unknowns, all SETTLED live

1. `mhc_class=1` (int): ACCEPTED. setup_for_alphafold.py rc 0.
2. Bare J allele `TRAJ42*01`: ACCEPTED. setup rc 0, no gene-table rejection. No
   `traj_allele`/`trbj_allele` data-model change needed.
3. HLA `A*02:01`: resolved to the `A0201` template internally. setup rc 0.
4. Output chain layout (was the big unknown): TCRdock writes the docking region as a
   SINGLE merged PDB chain `A` (~3284 atoms ~410 residues = MHC groove + peptide + Vα +
   Vβ), NOT the A..E multichain we assumed. There is nothing to "remap to A=TCRa
   B=TCRb...". The block structure is instead recoverable from `target_chainseq` in
   `{key}_final.tsv`, four `/`-joined segments in order **0=MHC, 1=peptide, 2=TCRalpha,
   3=TCRbeta**. TCRdock also emits 3 template models per target (indices 0/1/2), each a
   `..._model_1_model_2_ptm.pdb` + `_plddt.npy` + `_predicted_aligned_error.npy` +
   `_ptm.npy`.

### The discrimination metric (this is TCRdock's real output)

`{key}_final.tsv` carries `model_2_ptm_plddt`, `model_2_ptm_pae`, per-block
`model_2_ptm_plddt_{b}` and the inter-block PAE matrix `model_2_ptm_pae_{i}_{j}`. The
recognition signal is the peptide (block 1) to TCR (blocks 2,3) PAE plus the peptide
pLDDT. Calibrated live (model _0):

| metric | cognate GILGFVFTL | scramble TFVFGLIGL |
| --- | --- | --- |
| peptide pLDDT (plddt_1) | 86.24 | 64.55 |
| pae peptide->TCRa (pae_1_2) | 11.95 | 20.86 |
| pae peptide->TCRb (pae_1_3) | 10.49 | 20.77 |
| pae TCRa->peptide (pae_2_1) | 10.55 | 12.49 |
| pae TCRb->peptide (pae_3_1) | 8.10 | 11.44 |
| overall pae (model_2_ptm_pae) | 7.06 | 9.24 |

Lower peptide<->TCR PAE and higher peptide pLDDT for the true cognate TCR. The cleanest
scalar null is the peptide<->TCR interface PAE (e.g. mean of pae_1_2 and pae_1_3):
cognate ~11.2 vs scramble ~20.8.

### Wiring decision still open (surfaced to user, not yet coded)

The registry lists tcrdock as `qc_metric="cdr3_peptide"`, whose qc path
(`qc.load_chains` + `common_checks` expecting chains {A,B,C,D,E}) CANNOT consume
TCRdock's single-chain PDB. The natural wiring is to treat tcrdock like affinetune:
adapter reads the peptide<->TCR interface PAE from `final.tsv`, writes a one-float score
(inverted so higher = more recognized, since lower pae = better), and QC uses
`verdict_binding` against a tcrdock-specific `scramble_threshold` from the null above.
That means changing tcrdock's registered `qc_metric` from `cdr3_peptide` to a
PAE/binding style. Alternative: keep it structural and teach the `cdr3_peptide` path to
read the single-chain-plus-block layout. Pending user's call before wiring.
