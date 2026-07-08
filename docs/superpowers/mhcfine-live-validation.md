# mhcfine live validation checklist

Status: pivoted to mhcfine first because it is PyTorch (OpenFold-based), so it installs
on modern Colab, unlike tcrdock's legacy JAX fork. Builder reworked offline (merged).
Repo: https://bitbucket.org/abc-group/mhc-fine (paper PMC10705405, fine-tuned AlphaFold
in PyTorch).

## Step 1 (DONE, offline)

`mhcfine_inputs.build(construct_fasta)` emits `{protein_sequence, peptide_sequence}`
cognate + scramble (MHC heavy = chain C, peptide = chain E), matching MHC-Fine's real
inputs (unique_id, protein_sequence, peptide_sequence). Multi-value HLA now normalized
upstream in `annotate()` (shared `construct_io.normalize_hla`) so all 8 seed clonotypes
resolve a clean allele. 100/100 green.

## Contract (from their Inference.ipynb)

- Inputs: `unique_id` (str), `protein_sequence` (MHC heavy AA), `peptide_sequence` (8-11).
- Setup: `git clone` the bitbucket repo; `gdown` `mhc_fine_weights.pt` -> `data/model/`;
  `chmod +x a3m_generation/msa_run` (bundled MSA binary, path data/msa/{id}/mmseqs/aggregated.a3m).
- MSA is REQUIRED, built by `preprocess.get_a3m()` (msa_run; no large local DB, internet).
- Run: `model.Model().inference(np_sample, unique_id)` -> `./output/{unique_id}.pdb` + mean_plddt.

## Step 2: drive their Inference.ipynb (do NOT reconstruct blind)

Their notebook holds the real get_a3m + sample-building + inference code we do not have
verbatim. Plan: download their `Inference.ipynb`, upload to Colab, read the actual cells
live, then append ONE batch cell that loops our embedded TARGETS
(`scratchpad/mhcfine_targets.json`, 16 = 8 cognate + 8 scramble across GILGFVFTL,
NLVPMVATV, GLCTLVAML, KLGGALQAK) calling their get_a3m + inference and saving each PDB.

## Step 2 LIVE RUN 1 (2026-07-08, DONE, T4): pipeline proven end to end

Two real folds succeeded on Colab (stock GPU runtime, T4 15GB):
- Their known good example (6VRN, 179 aa alpha1alpha2 groove, peptide HMTEVVRHC):
  mean_plddt 98.18, mean_masked_plddt 96.93, output/EX_6VRN_A.pdb written.
- Our real R2S target clon0_cognate (full 275 aa ectodomain from hla_ectodomains.json,
  HLA-A*02:01, peptide GILGFVFTL), fresh msa_run MSA: mean_plddt 97.77,
  mean_masked_plddt 95.24, output/clon0_cognate.pdb written.

RESOLVED unknowns:
1. Deps: on stock Colab (torch 2.11, py3.12) the ONLY extra deps are
   `gdown biopython ml_collections dm-tree einops`. `from src import preprocess, model`
   then imports clean. (Contrast tcrdock, whose legacy JAX stack does not import.)
2. `msa_run` (bundled 17 MB binary) works on Colab as is after `chmod +x`; it queries
   its own remote MSA endpoint (no local DB), one aggregated.a3m per unique_id, cached
   under data/msa/{id}/mmseqs/. Re-runs with the same id skip the MSA.
3. Full 275 aa ectodomain folds fine (plddt ~98). NO need to trim to the alpha1alpha2
   groove; the R2S mhcfine builder's full-heavy-chain output (chain C) is compatible as is.

## First-try recipe (bake this in; it needs NO runtime restart)

The numpy trap, and the fix (re-validated live 2026-07-09, end-to-end pose). Colab ships
NumPy 2.x, and the AF2-derived code uses removals from that release: `np.string_`,
`np.unicode_`, `np.float_`, `np.complex_`, `np.bool8`, and `np.sum(generator)` (numpy 2
turned the last into a hard `TypeError`). The OBVIOUS fix, `pip install "numpy<2"`, is a
TRAP: on today's Colab image the downgrade poisons numpy's OWN compiled `mtrand.so` with a
dtype-size ABI wall ("numpy.dtype size changed, Expected 96 got 88") that resurfaces at
`model.inference` (lazy `numpy.random` import) and survives `--force-reinstall` and kernel
restarts. Wasted a long session on it.

The ROBUST recipe KEEPS stock numpy 2 and shims the removals in Python BEFORE
`from src import ...`. No downgrade, no ABI wall, no restart. Also the preprocess template
path needs the `kalign` binary (apt), else `kalign.py` builds a command with a None path ->
`TypeError: sequence item 0: expected str instance, NoneType`. Ordered cells for a cold
runtime:

1. `git clone https://bitbucket.org/abc-group/mhc-fine.git` ; `cd mhc-fine`.
2. deps, keeping stock numpy: `pip install -q gdown biopython ml_collections dm-tree einops`
   and `apt-get -qq install -y kalign`. Do NOT touch numpy.
3. numpy-2 shim, THEN import: `np.string_=np.bytes_; np.unicode_=np.str_; np.float_=
   np.float64; np.complex_=np.complex128; np.bool8=np.bool_`; wrap `np.sum` so a generator
   arg is `list()`-ed first; then `import torch; from src import preprocess, model`.
4. `gdown` weights (id 1gz8uF8DKE0CzyX_WeDGOX7xP69LjpaZT, 388 MB) -> data/model/
   mhc_fine_weights.pt ; `chmod +x a3m_generation/msa_run`.
5. Per target: `get_a3m(prot,a3m,uid)` -> `preprocess_for_inference(prot,pep,a3m)` ->
   `model.Model().inference(np_sample, uid)` -> ./output/{uid}.pdb + {mean_plddt}.

Live proof (2026-07-09): the A*02:01 + GILGFVFTL cognate folded under stock numpy 2.0.2 to
`{'mean_plddt': 97.77, 'mean_masked_plddt': 95.24}`, output/c0_cognate.pdb = chains A+B,
2314 atoms. This exact shim+fold sequence is what `notebook.py::_mhcfine_notebook` now
emits.

Playwright measurement trap (also cost time): Colab keeps a SEPARATE cross-origin
googleusercontent output iframe PER cell. A scraper that scans all frames for "Error" will
grab a STALE traceback from an earlier failed cell and declare a still-running fold dead.
Trust the target cell's own inline stdout stream (RESULT_OK) or a fresh on-disk check
(`glob output/*.pdb`), not "any frame contains Error". Also: Monaco auto-closes brackets, so
bracket-heavy one-liners corrupt on `keyboard.type` AND on paste; deliver such payloads as
base64 and `exec(base64.b64decode('...').decode())` (pure alphanumeric, nothing to corrupt).

Colab driving note (Playwright): code-cell output renders in a cross-origin
googleusercontent iframe, so read it via `page.frames()`+`frame.evaluate(innerText)`, not
the main-DOM `.output_subarea` (that only carries plain stdout streams). The reliable run
primitive is: click cell body -> `Escape` -> `b` (new cell) -> `keyboard.type(one_physical_line)`
-> `Control+Enter`. Keep the fold code on ONE physical line (semicolons, ternary for the
`if`) so monaco auto-indent cannot corrupt it. A clean one-shot notebook lives at
`/Users/fzd181/.playwright-mcp/mhcfine_fold.ipynb`.

## Open unknowns (next session)

1. Output PDB chain labeling: mhcfine writes chains A (MHC) and B (peptide) from
   `constants.temp_chains=["A","B"]`. For R2S `qc_structure` peptide_groove, remap their
   A->C (MHC), B->E (peptide). Confirm atom records on a downloaded PDB.
2. Batch the remaining 16 (8 cognate + 8 scramble) and calibrate the scramble threshold.
   WATCH OUT: mean_plddt alone likely will NOT separate cognate from scramble (MHC-Fine
   confidently places almost any 9-mer in the groove). The discriminator must be groove
   geometry / peptide-MHC contact pattern (peptide_groove), not plddt.

## Step 2 LIVE RUN 2 (2026-07-08, DONE, T4): batch 16 + scramble calibration

All 16 folded (8 cognate + 8 scramble, 16/16 ok, 0 fail; clon1/3/5/7 are exact
duplicates of clon0/2/4/6, deterministic seed so identical metrics). Output chains
confirmed A (MHC, 275 res) + B (peptide, 9 res). Metrics + PDBs pulled back
(scratchpad/mhcfine_results.json, scratchpad/mhcfine_pdbs/).

Calibration result (the payoff, and it is a NEGATIVE that matters):

| epitope   | cog masked_plddt | scr masked_plddt | cog res-in-groove | scr res-in-groove | cog contacts | scr contacts |
|-----------|------------------|------------------|-------------------|-------------------|--------------|--------------|
| GILGFVFTL | 95.24            | 90.54            | 9                 | 8                 | 139          | 142          |
| NLVPMVATV | 93.99            | 92.15            | 8                 | 8                 | 123          | 129          |
| GLCTLVAML | 93.00            | 89.89            | 8                 | 9                 | 121          | 129          |
| KLGGALQAK | 91.97            | 90.72            | 8                 | 7                 | 130          | 125          |

Conclusions:
1. mean_plddt (global) does NOT separate cognate from scramble (all ~97.0-97.8).
2. Peptide-MHC heavy-atom contact count (<4A) does NOT separate them either: scramble
   peptides seat just as deep in the groove (often MORE contacts). MHC-Fine imposes the
   canonical groove pose on ANY 9-mer. This is exactly the project thesis: a pose is not
   proof of binding.
3. Only mean_masked_plddt (peptide-region confidence) is directionally correct
   (cognate > scramble on 4/4 pairs), but the margin is small (1.2 to 4.7) and absolute
   scramble values stay high (~90-92) with cross-clonotype overlap. NOT a reliable
   absolute threshold.

CALIBRATION VERDICT for the mhcfine executor: peptide_groove (contact-based) CANNOT gate
specificity, and neither can plddt. The mhcfine QC must stay conservative: report the pose
+ mean_masked_plddt as a weak, within-pair-only relative signal, and NEVER let mhcfine
output read as "confirms recognition". mhcfine is a POSE, full stop. If a discriminator is
wanted later, the candidate is a per-pair cognate-vs-scramble masked_plddt DELTA (not an
absolute cutoff), validated on more pairs.

## Paper replication check (pose accuracy vs crystal, 2026-07-08)

Distinct question from the scramble control. The MHC-Fine paper (PMC10705405) benchmarks
peptide POSE ACCURACY (backbone RMSD vs the experimental structure) on known binders, and
reports a median around 0.66 A. We folded 6VRN's real peptide (HMTEVVRHC) and compared our
pose to the 6VRN crystal (RCSB): superpose on the MHC (182 shared CA, MHC CA RMSD 0.563 A),
then peptide backbone(N,CA,C,O) RMSD = 0.542 A, CA-only 0.507 A. That REPLICATES the
paper's accuracy (slightly better on this case) and confirms our pipeline drives MHC-Fine
correctly.

This does NOT conflict with the scramble result. The paper claims accurate poses for
binders; it never claims to tell binders from non-binders. Our scramble control adds the
orthogonal QC fact: MHC-Fine seats a non-binder in the groove just as confidently, so the
pose alone is not evidence of specificity. Paper = "places a binder accurately"; R2S QC =
"...and places a non-binder just as well, so a pose is not proof". Complementary.

## Step 3: QC made honest (DONE, branch feat/honest-peptide-groove)

LIVE RUN 2 showed groove-contact and plddt both fail to separate cognate from scramble, so
the old `verdict_groove` (pose_reliable if contact beats a "scramble null") was misleading.
Rewritten to be honest:
- `verdict_groove(pose_atoms, clonotype_id, tool, confidence=None)` (threshold dropped):
  returns `pose_only` for any in-groove pose (placement, never a specificity claim) and
  `pose_failed` when no peptide is in the groove. `calibration_basis="pose_quality"`.
- Also fixed a latent gate bug the live run exposed: the common gate expected chains
  {C,D,E} for peptide_groove, but a real mhcfine pose has only MHC heavy (C) + peptide (E),
  no b2m/TCR, so every real pose would have failed the gate. Gate now expects {C,E}.
- report EVIDENCE maps `pose_only` and `pose_failed` to "pose (...)"; qc-agent prompt states
  a groove tool's verdict is pose-only. Suite 100/100.

peptide_groove is distinct from Protenix/tcrdock QC and never shared. If a real specificity
discriminator is wanted later, the candidate is a per-pair cognate-vs-scramble
masked_plddt DELTA (not an absolute threshold), validated on more pairs.

## Step 4: notebook adapter WIRED (DONE, branch feat/wire-mhcfine-notebook)

`tools/notebook.build_notebook("mhcfine", inputs)` now returns the real validated recipe
(no longer a NotImplementedError stub): clone, deps+kalign keeping STOCK numpy 2, a numpy-2
compat shim before `from src import` (see the corrected recipe above; superseded the
original numpy<2 pin, which was proven a dead end on 2026-07-09), weights + chmod msa_run,
then a fold loop over INPUTS ({key: {protein_sequence, peptide_sequence}}, e.g. cognate +
scramble) writing ./output/{key}.pdb. Other tools (tcrdock/affinetune/af3) keep the
fail-loud stub.

## Step 6: adapter recipe CORRECTED to stock-numpy-2 shim (DONE, branch fix/mhcfine-numpy2-shim)

Live end-to-end drive of the mhcfine-agent surfaced that the numpy<2 pin from Steps 1-5 is a
dead end on today's Colab image (mtrand ABI wall at inference). Replaced the pin with the
keep-stock-numpy-2 + Python shim recipe, proven live to fold the A*02:01+GILGFVFTL cognate
(plddt 97.77, chains A+B). Tests locked the lesson: `test_mhcfine_keeps_stock_numpy2_no_downgrade`
(no `numpy<2`/`numpy==1` in the emitted notebook) and `test_mhcfine_shim_precedes_import_so_no_restart`
(shim + kalign before `from src import`). Suite 105/105.

## Step 5: adapter EXPOSED to the executor (DONE, branch feat/expose-mhcfine-notebook)

MCP tool `build_fold_notebook(run_dir, clonotype_id, tool)`: loads the clonotype's FoldJob,
shapes the tool inputs (mhcfine -> mhcfine_inputs.build(construct_fasta), keys PREFIXED by
clonotype id so per-clonotype output/{cid}_cognate.pdb never collide), calls build_notebook,
writes `<run_dir>/notebooks/{cid}_{tool}.ipynb`, returns the path. Wired into the
mhcfine-agent executor (tool added to its allowed set + prompt: list_fold_jobs ->
build_fold_notebook -> Playwright upload/run/download -> record_fold_result; report not-run
on a fail-loud scaffold or a Colab error). Suite 104/104.

STILL PENDING (live drive only, no more offline code): run the mhcfine-agent end to end on
a real Colab session to confirm the Playwright upload/run/download loop. Env caveat: the run
dir must sit under a Playwright-allowed root (e.g. /Users/fzd181) for the browser file
upload, and the Step 2 monaco/clipboard + iframe-read lessons apply to the live drive.
