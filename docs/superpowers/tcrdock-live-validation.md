# tcrdock live validation checklist

Status: builder reworked to the gene-level TSV contract (offline, merged). Steps 2
to 4 need a Colab GPU session and the AF2 params, done together with the user on a
real class I 10x CSV.

## Step 1 (DONE, offline)

`tcrdock_inputs.build(clonotype, annotation)` emits cognate + scramble rows for the
10 column tcrdock TSV: organism, mhc_class, mhc, peptide, va, ja, cdr3a, vb, jb,
cdr3b. Guards missing peptide/HLA, enforces the column contract. 95/95 green.

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
3. HLA formatting: builder strips only a leading `HLA-`. Confirm tcrdock accepts
   2 field (`A*02:01`); check 3 field (`A*02:01:01`) and any lowercase input from the
   real annotation output. Extend `_mhc_allele` only against tcrdock's allele lookup.
4. Scramble degeneracy (pre-existing, in `construct_io.scramble_peptide`): for very
   short or homopolymer like peptides the reverse+rotate can equal the original,
   weakening the calibration null. Backlog; verify the real peptides are not degenerate.

## Step 3: Playwright drive

Validate the executor (`tcrdock-agent`) opens the notebook, runs all cells, waits,
downloads the CIF, and calls `record_fold_result(tool="tcrdock")`. It must report
not-run (never fabricate) if the notebook/adapter is absent.

## Step 4: calibration

Run cognate + scramble on the real CSV, read the CDR3-to-peptide contact for both,
and set the tcrdock scramble_threshold from that null. Distinct from Protenix's; never
shared. QC uses `qc_metric="cdr3_peptide"` for tcrdock.
