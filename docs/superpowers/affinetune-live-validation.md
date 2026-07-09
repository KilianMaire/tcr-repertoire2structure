# affinetune live validation checklist

Status: builder exists (offline, merged); the notebook is still the fail-loud stub. Steps
2 to 4 need a Colab GPU session and the fine-tuned params, done together with the user on
the real class I seed CSV. This tool shares tcrdock's environment blocker (JAX built for
an old Python), so the condacolab bring-up is the gate for BOTH.

## Real tool identity

affinetune maps to **phbradley/alphafold_finetune** (Motmaen et al. 2023 PNAS, "Peptide
binding specificity prediction using fine-tuned protein structure prediction networks").
It is a JAX AlphaFold fine-tune that scores whether a peptide is PRESENTED by a class I
(and class II) MHC. It returns a binding score, not a structure. This is why the registry
gives it `output_type="binding_score"` and `qc_metric="binding_score"`, and why the honesty
label is "predicted presentation", never "fold" or "structure".

## Step 1 (DONE, offline)

`affinetune_inputs.build(construct_fasta)` parses the construct and emits cognate +
scramble records `{mhc, b2m, peptide}`, each carrying its own scramble null (mirroring the
Protenix builder). Merged; covered by the offline suite.

Caveat surfaced while grounding the recipe (see Open unknowns): alphafold_finetune's pMHC
input is MHC-alpha + peptide only. b2m is NOT part of its target_chainseq. The builder keeps
b2m for now (harmless, other tools use it); the notebook adapter is what must drop b2m when
it assembles the chainseq. Do NOT strip b2m from the builder before the live run confirms it.

## Validation seed data

Reuse `data/validation_tcrdock_classI.csv` (8 class I clonotypes: GILGFVFTL/A*02:01,
NLVPMVATV/A*02:01, GLCTLVAML/A*02:01, KLGGALQAK/A*03:01). affinetune needs only the HLA
allele and the peptide, both present. A*02:01 is the best-supported allele for a first run
because alphafold_finetune ships a prebuilt A*02:01 template alignment (see Step 2).

## Step 2: real Colab cell (to replace the `raise NotImplementedError` stub)

Grounded against the current alphafold_finetune README (commands quoted verbatim):

1. `git clone https://github.com/phbradley/alphafold_finetune`.
2. Fetch + unpack the params/datasets bundle:
   `wget https://files.ipd.uw.edu/pub/alphafold_finetune_motmaen_pnas_2023/datasets_alphafold_finetune_v2_2023-02-20.tgz`
   then `tar -xzvf datasets_alphafold_finetune_v2_2023-02-20.tgz`.
   The pMHC params land at
   `datasets_alphafold_finetune/params/mixed_mhc_pae_run6_af_mhc_params_20640.pkl`.
3. Write a `targets.tsv` with the header
   `mhc  start  peptide  targetid  target_chainseq  templates_alignfile`
   where `target_chainseq` is `<MHC-alpha-seq>/<peptide>` (slash-separated, b2m omitted)
   and `templates_alignfile` points at the per-allele alignment TSV (the repo ships
   `examples/pmhc_hcv_polg_10mers/alignments/A0201_10mer_alignments.tsv` for A*02:01
   10-mers).
4. Run:
   `python run_prediction.py --targets targets.tsv --outfile_prefix run --model_names model_2_ptm_ft --model_params_files datasets_alphafold_finetune/params/mixed_mhc_pae_run6_af_mhc_params_20640.pkl --ignore_identities`
5. Read the predicted presentation score out of `run*_final.tsv` for cognate AND scramble.

## Environment bring-up (VALIDATED live 2026-07-09 on Colab T4, driver 580 / CUDA 12.8)

alphafold_finetune is JAX + an AlphaFold fork pinned to old wheels (jax 0.2.22 / jaxlib
0.1.72+cuda111, cp38). Today's Colab kernel is Python 3.12, so those wheels cannot install
in-kernel, exactly the wall tcrdock hit. The working pattern is the maintainers' own one
(their `alphafold_ft_colab_pipeline_v1.ipynb`): DO NOT switch the kernel. Build a conda
python 3.8 env off to the side and SHELL OUT to it. No condacolab kernel-restart is needed.
This same env unblocks tcrdock (same JAX/AF family). Proven recipe, in order:

1. Miniforge + a py3.8 env (no kernel restart):
   `wget -qO /tmp/mf.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh && bash /tmp/mf.sh -b -p /opt/conda`
   then `/opt/conda/bin/mamba create -y -n af python=3.8` (yields Python 3.8.20).
2. `git clone https://github.com/phbradley/alphafold_finetune /content/aff` and
   `/opt/conda/envs/af/bin/pip install -r /content/aff/requirements_colab_python38_v2.txt`
   (this pins the jaxlib-cuda111 wheel; ~1.5 GB incl. torch 1.10, TF 2.5; ~7 min).
3. Base AF params `params_model_2_ptm.npz` (dropbox, in the notebook) into
   `/content/alphafold_params/params/`, plus the dataset bundle (params pkl + example
   alignments) as in Step 2 of the run recipe above.
4. THE CUDA FIX (this is what makes it run, not crash): the py3.8 env has jaxlib-cuda111 but
   NOT the CUDA 11 runtime libs, so run_prediction aborts (RC 134) with
   `Could not load dynamic library 'libcudart.so.11.0' / 'libcublas.so.11'`. Install them
   into the env: `/opt/conda/bin/mamba install -y -n af -c conda-forge cudatoolkit=11.1 cudnn`.
   That provides `libcudart.so.11.0` (-> 11.1.74) and `libcublas.so.11` (-> 11.3). AlphaFold's
   Evoformer has no convolutions, so the exact cudnn version is not load-critical (mamba
   pulled cudnn 9.x and inference still runs).
5. THE LD_LIBRARY_PATH RECIPE (the second gotcha): every shell-out MUST set
   `LD_LIBRARY_PATH=/opt/conda/envs/af/lib:/usr/lib64-nvidia` and PREPEND, never replace.
   The env lib gives libcudart/libcublas; `/usr/lib64-nvidia` gives the driver `libcuda.so.1`.
   Drop the driver dir and you get `cuInit: UNKNOWN ERROR (303)` and a silent CPU fallback
   (very slow), which looks like a hang. With both dirs present, jax reports
   `GpuDevice(id=0)` and `Successfully opened dynamic library libcudart.so.11.0`.
6. Run the prediction as a shell-out to that python (Step 4 of the run recipe), e.g.
   `cd /content/aff && LD_LIBRARY_PATH=/opt/conda/envs/af/lib:/usr/lib64-nvidia /opt/conda/envs/af/bin/python run_prediction.py --targets one.tsv --outfile_prefix one_test --model_names model_2_ptm_ft --model_params_files datasets_alphafold_finetune/params/mixed_mhc_pae_run6_af_mhc_params_20640.pkl --data_dir /content/alphafold_params/ --ignore_identities`.
   Proven: the maintainers' example ran end-to-end on the GPU for ~13 min (no crash) once
   Steps 4-5 were in place.

Playwright reading caveat (learned here): Colab renders shelled-out stdout in a cross-origin
output iframe. Short cell outputs read fine from the accessibility snapshot right after they
finish, but a long AF cell's streamed stderr goes STALE in the snapshot (you re-read an old
cell's traceback). To read a real result, DO NOT trust streamed stdout; either check the
on-disk file, or make the cell's LAST EXPRESSION a DataFrame (`pd.read_csv(out, sep='\t')`)
so Colab renders it as a table in the MAIN document, which the snapshot can read.

## Calibration fold (VALIDATED live 2026-07-09 on Colab T4)

A real cognate + scramble prediction on the class I seed peptide, run end to end
through the rebuilt env. Cognate `GILGFVFTL` (flu M1, A*02:01) vs its deterministic
scramble `TFVFGLIGL` (the pipeline null from `construct_io.scramble_peptide`). Both
targets use the true A*02:01 alpha (175 aa, from the polg example), a b2m-free
chainseq (`alpha/peptide`), and the shipped 9-mer alignment file
`examples/tiny_pmhc_finetune/alignments/A0203_alignments.tsv` with `--ignore_identities`.
Output `cal_final.tsv`, column `model_2_ptm_ft_pae`:

| target | model_2_ptm_ft_pae | peptide pLDDT | pae_1_0 (pep to MHC) |
| --- | --- | --- | --- |
| cognate GILGFVFTL | 1.017 | 91.2 | 4.50 |
| scramble TFVFGLIGL | 2.252 | 36.9 | 24.10 |

Reads live (confirms all three open unknowns):
1. Score column IS `model_2_ptm_ft_pae` and LOWER = presented. Cognate 1.02 beats
   scramble 2.25 (2.2x). Two independent structural corroborations: the peptide
   pLDDT collapses (91.2 to 36.9) and the peptide-to-MHC pae blows up (4.50 to
   24.10) for the scramble. The adapter must invert (score = -pae) before it feeds
   `verdict_binding`, which treats HIGHER as more presented.
2. b2m-free chainseq (alpha + peptide only) folds correctly.
3. templates_alignfile for a 9-mer: RESOLVED. The `tiny_pmhc_finetune` bundle ships
   per-allele 9-mer alignment files (`target_len` 184 = 175-aa alpha + 9-mer), each
   with its own committed template PDBs. Reuse the closest A2-family file
   (`A0203_alignments.tsv`) with the true A*02:01 alpha; `run_prediction.py:124`
   asserts only `target_len == len(query)` (a LENGTH check), and `--ignore_identities`
   neutralizes the donor-allele mismatch. No alignment is hand-faked.

Calibration threshold: cognate pae 1.02 clears the scramble null 2.25 with a
1.24-unit margin. The affinetune `scramble_threshold` (a binding_score null) is set
from this pair, on the inverted score. Distinct from the Protenix/tcrdock/mhcfine
thresholds; never shared.

## Open unknowns to resolve on the first real run (do NOT guess-fix in the builder)

1. Output score column + DIRECTION. RESOLVED from run_prediction.py source (lines 162-193):
   the output file is `{outfile_prefix}_final.tsv` (tab-separated); it writes
   `{model_name}_pae` (mean predicted aligned error over the complex, e.g.
   `model_2_ptm_ft_pae`) and per-chain-pair `{model_name}_pae_{chain1}_{chain2}`. The score
   IS the PAE and LOWER = more likely presented (the bundled BinderClassifier has slope
   -7.9, so higher pae -> lower binding probability). So the result reader must read
   `*_final.tsv`, take `<model>_pae`, and INVERT it (e.g. score = -pae) before it feeds
   `verdict_binding`, which treats HIGHER as more presented. CAPTURED live 2026-07-09:
   cognate pae 1.02 vs scramble pae 2.25 (see the Calibration fold section above).
2. b2m in the chainseq. The A*02:01 example row is MHC-alpha/peptide with no b2m. Confirm
   the pMHC model wants alpha only; if so the notebook assembles chainseq from `mhc` +
   `peptide` and ignores `b2m`.
3. templates_alignfile for arbitrary allele/peptide length. RESOLVED for A*02:01 9-mers
   (see the Calibration fold section: reuse the `tiny_pmhc_finetune` 9-mer alignment
   files, `target_len` 184). Still open for other lengths/alleles (A*03:01 for
   KLGGALQAK, class II): pick the shipped file whose `target_len` matches the query
   length and whose donor allele is closest, with `--ignore_identities`; the bundle
   ships class I files A02xx/A03xx/A11xx/A23xx/A24xx/A26xx/A29xx/A30xx/A31xx/A33xx/A68xx
   and B-alleles, plus DRB/DPB class II template files. Do not hand-fake an alignment.
4. Scramble degeneracy (pre-existing, `construct_io.scramble_peptide`): verify the real
   peptides are not homopolymer-like so the scramble null is a genuine non-binder.

## Step 3: Playwright drive

Validate the executor (`affinetune-agent`) opens the notebook, runs the cells by keyboard
(Ctrl+Enter / Shift+Enter), waits, reads the score, calls `record_fold_result(tool="affinetune")`,
then releases the runtime (Runtime > Disconnect and delete runtime). It must report not-run
and never fabricate a score if the notebook/adapter/env is unavailable.

## Step 4: calibration

Score cognate + scramble with affinetune, read the presentation score for both, and set the
affinetune `scramble_threshold` (a binding_score null) from that pair. Distinct from the
Protenix/tcrdock/mhcfine thresholds; never shared. QC uses `qc_metric="binding_score"` and
`verdict_binding` labels the row "predicted presentation", not geometry.
