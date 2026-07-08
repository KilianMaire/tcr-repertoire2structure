"""Build a self-contained Protenix fold notebook with the inputs embedded.

No Drive upload, no MSA bundle: the JSONs are written inline in cell 2, Protenix
queries its own MSA server, and results zip for download. Pick one representative
clonotype per distinct epitope (cognate + scramble each).

Usage:
  python scripts/build_colab_notebook.py <protenix_inputs_dir> <out.ipynb> [id1 id2 ...]
"""
from __future__ import annotations
import json
import sys
from pathlib import Path


def code(src):
    return {"cell_type": "code", "metadata": {}, "execution_count": None,
            "outputs": [], "source": src}


def md(src):
    return {"cell_type": "markdown", "metadata": {}, "source": src}


def main():
    in_dir = Path(sys.argv[1])
    out_ipynb = sys.argv[2]
    ids = sys.argv[3:] if len(sys.argv) > 3 else None

    files = sorted(in_dir.glob("*.json"))
    if ids:
        files = [f for f in files if any(f.stem.startswith(i) for i in ids)]
    inputs = {f.name: json.loads(f.read_text()) for f in files}
    seeds = "101"

    nb = {"nbformat": 4, "nbformat_minor": 0,
          "metadata": {"accelerator": "GPU", "colab": {"provenance": []},
                       "kernelspec": {"name": "python3", "display_name": "Python 3"}},
          "cells": [
        md(["# Repertoire2Structure: Protenix fold (cognate + scramble)\n",
            "Self-contained. Runtime > Change runtime type > GPU (A100/L4), then Run all.\n",
            "Folds the embedded TCR-pMHC constructs; the scramble controls calibrate the ",
            "skeptical QC (a fold does not confirm specificity)."]),
        code(["#@title 1. Install Protenix + show GPU\n",
              "import subprocess\n",
              "print(subprocess.run('nvidia-smi --query-gpu=name,memory.total --format=csv,noheader',",
              " shell=True, capture_output=True, text=True).stdout)\n",
              "r = subprocess.run('pip install -q protenix', shell=True, capture_output=True, text=True)\n",
              "print(r.stderr[-800:]); print('install done')"]),
        code(["#@title 2. Write the embedded fold inputs\n",
              "import json, os\n",
              f"INPUTS = {json.dumps(inputs)}\n",
              "os.makedirs('inputs', exist_ok=True)\n",
              "for name, obj in INPUTS.items():\n",
              "    json.dump(obj, open(f'inputs/{name}', 'w'))\n",
              "print('wrote', len(INPUTS), 'inputs:', sorted(INPUTS))"]),
        code(["#@title 3. Fold each input (MSA server, 1 seed), zip results\n",
              "import glob, os, subprocess, shutil\n",
              "os.makedirs('out', exist_ok=True)\n",
              "for f in sorted(glob.glob('inputs/*.json')):\n",
              "    name = os.path.basename(f)[:-5]\n",
              "    if os.path.isdir(f'out/{name}'):\n",
              "        print('skip', name); continue\n",
              f"    cmd = f'protenix pred -i {{f}} -o out/{{name}} -s {seeds} '\\\n",
              "          '-n protenix_base_default_v1.0.0 --use_msa true --use_default_params true'\n",
              "    print('RUN', name, flush=True)\n",
              "    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)\n",
              "    if r.returncode != 0:\n",
              "        print('FAIL', name, (r.stdout+r.stderr)[-1500:]); continue\n",
              "    print('done', name)\n",
              "subprocess.run('zip -qr rep2struct_folds.zip out', shell=True)\n",
              "print('zipped ->', os.path.getsize('rep2struct_folds.zip'), 'bytes')"]),
        code(["#@title 4. Download results\n",
              "from google.colab import files\n",
              "files.download('rep2struct_folds.zip')"]),
      ]}
    Path(out_ipynb).write_text(json.dumps(nb, indent=1))
    print(f"wrote {out_ipynb} with {len(inputs)} folds: {sorted(inputs)}")


if __name__ == "__main__":
    main()
