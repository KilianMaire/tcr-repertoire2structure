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

## Open unknowns to settle on the first real run

1. Exact deps beyond torch (ml_collections, dm-tree, biopython, gdown): the notebook has
   no explicit pip block, so imports may need pinning live.
2. `msa_run` binary behaviour on Colab (arch, network endpoint) and per-target MSA time.
3. Output PDB chain labeling for the pMHC pose -> remap to canonical C (MHC), E (peptide)
   for `qc_structure` peptide_groove; confirm on the first real output.

## Step 3: QC + calibration

mhcfine `qc_metric = peptide_groove`: peptide-in-groove-vs-scramble as the skeptical
verdict, mean_plddt reported alongside (never the sole basis). Calibrate the mhcfine
scramble threshold from the cognate-vs-scramble groove contact on this seed. Distinct
from Protenix/tcrdock; never shared. mhcfine is a POSE, never a proof of TCR recognition.
