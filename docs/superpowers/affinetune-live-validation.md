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

## Environment blocker (shared with tcrdock, this is the gate)

alphafold_finetune is JAX + an AlphaFold fork pinned to an old Python (3.8/3.9). Today's
Colab is Python 3.12 + CUDA 12; its jaxlib will not match the pinned jax, exactly the wall
tcrdock hit (`tcrdock-live-validation.md`, live run 1). So affinetune's live drive is NOT a
clean pip path like mhcfine was. It needs `condacolab` to bring up a py3.9 env, then
`pip install` alphafold_finetune's requirements and a jaxlib that matches the Colab GPU.
Solve the condacolab bring-up ONCE and both tcrdock and affinetune unblock. Until then the
executor correctly reports not-run; nothing is faked.

## Open unknowns to resolve on the first real run (do NOT guess-fix in the builder)

1. Output score column + DIRECTION. `run*_final.tsv` carries the predicted metric; the
   likely field is the peptide-MHC PAE (Motmaen uses a PAE-derived score) where LOWER means
   more likely presented. Confirm the exact column name and whether to invert it before it
   feeds `verdict_binding` (which treats HIGHER as more presented). The result reader must
   read the real column, never assume one. (Fetch of run_prediction.py was rate-limited;
   settle this by reading the actual output file on the first run.)
2. b2m in the chainseq. The A*02:01 example row is MHC-alpha/peptide with no b2m. Confirm
   the pMHC model wants alpha only; if so the notebook assembles chainseq from `mhc` +
   `peptide` and ignores `b2m`.
3. templates_alignfile for arbitrary allele/peptide length. The repo ships a prebuilt
   A*02:01 10-mer alignment; 9-mers, other alleles (A*03:01 for KLGGALQAK), and class II
   need their own alignment file. Find alphafold_finetune's alignment-generation helper
   before running anything outside A*02:01 10-mers; do not hand-fake an alignment.
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
