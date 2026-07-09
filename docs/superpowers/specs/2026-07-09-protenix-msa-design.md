# Protenix pre-fold MSA (v1) — design

## Problem

Protenix is the default fold workhorse, but it folds MSA-free (single sequence).
The first live fold proved the cost: on the flu M1 GILGFVFTL TCR (a well
characterized binder), 4 of 5 Protenix samples had zero CDR3(Vbeta) to peptide
contact for both cognate and scramble. The TCR does not dock onto the peptide.
A single sequence gives the model no co-evolution signal to place the interface,
so the fold is close to a guess and every QC verdict is uninformative.

The fix is to give Protenix a real MSA, computed outside the fold runtime so a
slow or throttled MSA server can never wedge the GPU fold (the failure that made
the project switch to MSA-free in the first place: Protenix's own `--use_msa true`
server throttled and wedged the Colab kernel mid batch).

This slice covers Protenix only. tcrdock and affinetune are template based (no
deep MSA), mhcfine builds its own MSA via its bundled `msa_run` binary, and af3
is not wired.

## Verified interface

Protenix consumes a precomputed MSA per protein chain via the input JSON
(`docs/infer_json_format.md`):

- `unpairedMsaPath`: absolute path to a per-chain non-pairing `.a3m`.
- `pairedMsaPath`: absolute path to a pairing `.a3m` (out of scope for v1).

If neither is present the chain folds without MSA. So the mechanism is: write
`unpairedMsaPath` into each `proteinChain` of the input JSON, pointing at an a3m
produced before the fold.

## Decisions

- Compute the MSA in a Colab CPU cell inside `_protenix_notebook`, before the
  fold cell. Not embedded in the notebook, not computed at prep. The cell is
  separate from the `protenix pred` GPU cell, so a slow MSA server delays that
  cell but cannot wedge the fold.
- Unpaired MSA only, one a3m per protein chain. For a TCR-pMHC there is no real
  inter-chain co-evolution signal, so pairing adds little. Paired MSA is deferred
  to v2.
- MSA server: the ColabFold MMseqs2 API (`api.colabfold.com`), the robust public
  server, via the `colabfold` package's `run_mmseqs2` (exact import pinned at
  implementation).
- Report honesty is in scope for v1: the report states the actual per-clonotype
  MSA basis (colab_cpu or none), driven by a manifest the fold produces, not the
  prep-time intent.

## Architecture

`_protenix_notebook` gains one cell and one edit, in this cell order:

1. INPUTS (Protenix JSON per record, no MSA paths). Unchanged.
2. Install Protenix + GPU check. Unchanged.
3. NEW: MSA cell (Colab CPU). Detailed below.
4. Write `inputs/{key}.json`, now with `unpairedMsaPath` injected per protein
   chain. Edit of the existing write cell.
5. Fold: `protenix pred ... --use_msa false`. Protenix consumes the provided a3m
   paths and does not run its own search.
6. Repatriation. Unchanged, plus the MSA manifest travels in the zip.

### The MSA cell

- Collect the unique protein chain sequences across all INPUTS records. Cognate
  and scramble share chains A, B, C, D; only the peptide (E) differs, so this is
  4 sequences, not 8.
- Exclude the peptide by length: skip any sequence with `len(seq) < 20`. This
  covers the 9-mer peptide and any short chain, and needs no chain-id coupling.
- For each remaining unique sequence, query the ColabFold MMseqs2 API (unpaired)
  and write `/content/msa/<seq_hash>.a3m`. Build a `seq -> a3m path` map.
- Walk INPUTS and set `proteinChain["unpairedMsaPath"]` for every chain whose
  sequence has an a3m.
- Graceful degradation per chain: if the API fails for a sequence (timeout,
  throttle, error), that chain gets no path and folds MSA-free; the cell logs it
  and continues. The run never dies on MSA.
- Write the manifest per clonotype, named from the INPUTS key prefix (the
  clonotype id), e.g. `out/{cid}_msa_manifest.json`, so several clonotypes
  unzipping into the same run dir do not collide. It holds, per chain,
  `{got_msa: bool, n_seqs: int}` plus a summary, and sits under `out/` so the
  existing repatriation zip carries it back.

### Fold command

Keep `--use_msa false`. The one residual unknown is whether Protenix honors
`unpairedMsaPath` while `--use_msa false`. Verified live at implementation; if it
does not, drop the flag (the provided paths still take precedence over a fresh
search). This is the only interface risk and it is checked before declaring the
slice done.

### Report honesty

The repatriated per-clonotype `{cid}_msa_manifest.json` lands in the run dir.
`render_report` reads it and shows, per clonotype, the actual MSA basis:
`colab_cpu (k/n chains)` when MSA was used, or `MSA-free` when it was not,
overriding the prep-time `msa_basis`. The report never claims an MSA that a chain
did not receive.

### msa.py

Consequence of computing in the notebook: `msa.py`'s prep-time runner slots
(`local_runner`, `colab_runner`) are not used by this slice. They remain for a
future local mmseqs option. `msa.py` is unchanged; its prep-time `msa_basis`
stamp is superseded for Protenix by the manifest-driven report value.

## Testing

Unit (offline, deterministic):

- The Protenix notebook contains the MSA cell (colabfold / run_mmseqs2,
  `/content/msa`, writes `unpairedMsaPath`), and the write cell injects
  `unpairedMsaPath` into protein chains.
- Sequence dedup: cognate + scramble yield 4 unique chain sequences, not 8.
- Peptide exclusion: a `< 20` residue chain gets no MSA path.
- Graceful fallback marker present (a chain without an a3m folds MSA-free).
- Report renders the manifest-driven per-clonotype MSA basis (colab_cpu vs
  MSA-free), including the mixed and all-failed cases.

Live (the real success criterion, run on the A100 Colab as before):

- Re-fold the flu M1 clonotype WITH MSA. The cognate must now dock: median
  CDR3-peptide contact > 0 and > its scramble, in contrast to today's MSA-free
  run where 4 of 5 poses had zero contact. This is the proof the slice worked.

## Scope and non-goals

- Protenix only. Not tcrdock, affinetune, mhcfine, af3.
- Unpaired MSA only. Paired MSA deferred to v2.
- Peptide chain excluded from MSA.
- No local mmseqs. The ColabFold API is the only MSA source in v1.

## Risks

- `--use_msa false` may not honor `unpairedMsaPath`. Mitigation: verified live;
  drop the flag if needed.
- ColabFold API throttling. Mitigation: the MSA cell is decoupled from the GPU
  fold and degrades per chain to MSA-free; it never wedges the fold.
