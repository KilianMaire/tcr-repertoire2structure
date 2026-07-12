"""Convert construct FASTAs (chains A-E) into Protenix input JSONs.

For each cognate construct we also emit a scramble control: the peptide (chain E)
deterministically shuffled, everything else identical. The scramble is the null
the QC step calibrates against (Honesty Rule 2: a fold does not confirm
specificity; a cognate must beat its scramble on CDR3-peptide contact).

MSA paths are omitted, so Protenix queries its own MSA server at fold time
(`--use_msa true`), which removes the precompute-and-bundle step.

Usage:
  python scripts/build_protenix_inputs.py <run_dir_with_construct_fastas> <out_dir>
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))
from rep2struct.tools.construct_io import parse_fasta, scramble_peptide


def to_protenix(name, chains):
    order = ["A", "B", "C", "D", "E"]
    seqs = [{"proteinChain": {"sequence": chains[c], "count": 1, "id": [c]}}
            for c in order if c in chains]
    return [{"name": name, "sequences": seqs, "covalent_bonds": []}]


def main():
    run_dir = Path(sys.argv[1])
    out_dir = Path(sys.argv[2])
    out_dir.mkdir(parents=True, exist_ok=True)
    fastas = sorted(run_dir.glob("construct_*.fasta"))
    n = 0
    for fp in fastas:
        cid = fp.stem.replace("construct_", "")
        chains = parse_fasta(fp.read_text())
        if set("ABCDE") - set(chains):
            print(f"skip {cid}: missing chains {set('ABCDE') - set(chains)}")
            continue
        # cognate
        (out_dir / f"{cid}_cognate.json").write_text(
            json.dumps(to_protenix(f"{cid}_cognate", chains), indent=1))
        # scramble control (peptide shuffled)
        sc = dict(chains); sc["E"] = scramble_peptide(chains["E"])
        (out_dir / f"{cid}_scramble.json").write_text(
            json.dumps(to_protenix(f"{cid}_scramble", sc), indent=1))
        n += 1
        print(f"{cid}: cognate peptide {chains['E']} | scramble {sc['E']}")
    print(f"wrote {n} cognate + {n} scramble Protenix inputs to {out_dir}")


if __name__ == "__main__":
    main()
