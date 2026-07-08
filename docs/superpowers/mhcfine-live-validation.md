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

The ONLY ordering trap: Colab ships NumPy 2.x, and the AF2-derived code uses the removed
`np.string_`. Pin `numpy<2` and, critically, do it BEFORE importing numpy/torch/src, so
no kernel restart is needed. Also the preprocess template path needs the `kalign` binary
(apt), else `kalign.py` builds a command with a None path -> `TypeError: sequence item 0:
expected str instance, NoneType`. Ordered cells for a cold runtime:

1. `git clone https://bitbucket.org/abc-group/mhc-fine.git` ; `cd mhc-fine`.
2. FIRST, before any import: `pip install -q "numpy<2" gdown biopython ml_collections
   dm-tree einops` and `apt-get -qq install -y kalign`. (numpy 1.26.4; the opencv/tobler/
   pytensor "requires numpy>=2" pip warnings are irrelevant to mhcfine.)
3. `import torch; from src import preprocess, model` (loads numpy 1.26 fresh -> no restart).
4. `gdown` weights (id 1gz8uF8DKE0CzyX_WeDGOX7xP69LjpaZT, 388 MB) -> data/model/
   mhc_fine_weights.pt ; `chmod +x a3m_generation/msa_run`.
5. Per target: `get_a3m(prot,a3m,uid)` -> `preprocess_for_inference(prot,pep,a3m)` ->
   `model.Model().inference(np_sample, uid)` -> ./output/{uid}.pdb + {mean_plddt}.

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

## Step 3: QC + calibration (superseded by LIVE RUN 2 above)

mhcfine `qc_metric = peptide_groove`. LIVE RUN 2 showed groove-contact and plddt both fail
to separate cognate from scramble, so the peptide_groove verdict must be conservative
(pose only, never proof). Distinct from Protenix/tcrdock QC; never shared.
