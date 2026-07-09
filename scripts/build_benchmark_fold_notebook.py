"""Build the Protenix fold notebook for a benchmark seed batch.

Reads every construct FASTA emitted by run_benchmark_arm.py (cognate + decoys +
scramble, keyed {clonotype}__{epitope}) and embeds them into the validated
MSA-based Protenix notebook (tools.notebook.build_notebook), which dedups the MSA
by unique protein chain, so the shared TCR/MHC chains are searched once and reused
across a TCR's constructs. The user runs the notebook on Colab/H100 and returns the
out/ zip; scoring reads out/{key}/**/*.cif.

Usage:
  python scripts/build_benchmark_fold_notebook.py <constructs_dir> <out.ipynb>
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from rep2struct.tools.construct_io import parse_fasta
from rep2struct.tools.protenix_inputs import _to_protenix
from rep2struct.tools.notebook import build_notebook


def build_inputs(constructs_dir):
    inputs = {}
    for fasta in sorted(Path(constructs_dir).glob("*.fasta")):
        chains = parse_fasta(fasta.read_text())
        inputs[fasta.stem] = _to_protenix(fasta.stem, chains)
    return inputs


def main():
    constructs_dir, out_ipynb = sys.argv[1], sys.argv[2]
    inputs = build_inputs(constructs_dir)
    nb = build_notebook("protenix", inputs)
    Path(out_ipynb).write_text(json.dumps(nb, indent=1))
    print(f"embedded {len(inputs)} constructs -> {out_ipynb}")
    print("keys:", sorted(inputs))


if __name__ == "__main__":
    main()
