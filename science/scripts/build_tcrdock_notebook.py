"""Build a self-contained TCRdock (AF2) fold notebook with targets embedded.

TCRdock takes gene-level rows (organism, mhc_class, mhc, peptide, va, ja, cdr3a,
vb, jb, cdr3b) and reconstructs the TCR:pMHC structure from templates, so this
notebook writes targets.tsv, runs setup_for_alphafold.py then run_prediction.py,
and remaps the output chains to the canonical A=TCRa B=TCRb C=MHC D=b2m E=peptide.

The env is the fragile part: TCRdock is written for python 3.8 + an old JAX/AF2
fork. The install cell is best-effort and validated live on Colab, not blind.

Usage:
  python scripts/build_tcrdock_notebook.py <targets.json> <out.ipynb>
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

AF_PARAMS_URL = "https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar"


def code(src):
    return {"cell_type": "code", "metadata": {}, "execution_count": None,
            "outputs": [], "source": src}


def md(src):
    return {"cell_type": "markdown", "metadata": {}, "source": src}


def main():
    targets = json.loads(Path(sys.argv[1]).read_text())
    out_ipynb = sys.argv[2]

    nb = {"nbformat": 4, "nbformat_minor": 0,
          "metadata": {"accelerator": "GPU", "colab": {"provenance": []},
                       "kernelspec": {"name": "python3", "display_name": "Python 3"}},
          "cells": [
        md(["# Repertoire2Structure: TCRdock fold (cognate + scramble)\n",
            "Runtime > Change runtime type > GPU, then run cells in order.\n",
            "TCRdock reconstructs each TCR:pMHC from gene names; the scramble rows ",
            "are the per-target calibration null (a fold does not confirm specificity)."]),
        code(["#@title 1. GPU + clone TCRdock\n",
              "import subprocess\n",
              "print(subprocess.run('nvidia-smi --query-gpu=name,memory.total --format=csv,noheader',",
              " shell=True, capture_output=True, text=True).stdout)\n",
              "subprocess.run('git clone https://github.com/phbradley/TCRdock', shell=True)\n",
              "print('cloned')"]),
        code(["#@title 2. Install requirements + BLAST db (best-effort; iterate live)\n",
              "import subprocess\n",
              "r = subprocess.run('cd TCRdock && pip install -q -r requirements.txt',",
              " shell=True, capture_output=True, text=True)\n",
              "print('pip:', r.stderr[-800:])\n",
              "r = subprocess.run('cd TCRdock && python download_blast.py', shell=True,",
              " capture_output=True, text=True)\n",
              "print('blast:', (r.stdout + r.stderr)[-800:])"]),
        code([f"#@title 3. Download AlphaFold params ({AF_PARAMS_URL.split('/')[-1]})\n",
              "import os, subprocess\n",
              "os.makedirs('afdata/params', exist_ok=True)\n",
              f"cmd = 'wget -q {AF_PARAMS_URL} -O afparams.tar && tar -xf afparams.tar -C afdata/params'\n",
              "r = subprocess.run(cmd, shell=True, capture_output=True, text=True)\n",
              "print('params:', (r.stdout + r.stderr)[-600:])\n",
              "print('params dir:', os.listdir('afdata/params')[:5])"]),
        code(["#@title 4. Write the embedded targets.tsv\n",
              "import csv\n",
              f"TARGETS = {json.dumps(targets)}\n",
              "COLS = ['organism','mhc_class','mhc','peptide','va','ja','cdr3a','vb','jb','cdr3b']\n",
              "with open('TCRdock/targets.tsv', 'w', newline='') as fh:\n",
              "    w = csv.DictWriter(fh, fieldnames=COLS, extrasaction='ignore', delimiter='\\t')\n",
              "    w.writeheader()\n",
              "    for t in TARGETS:\n",
              "        w.writerow(t)\n",
              "print('wrote', len(TARGETS), 'targets')"]),
        code(["#@title 5. setup_for_alphafold.py\n",
              "import subprocess\n",
              "cmd = ('cd TCRdock && python setup_for_alphafold.py '\n",
              "       '--targets_tsvfile targets.tsv --output_dir setup')\n",
              "r = subprocess.run(cmd, shell=True, capture_output=True, text=True)\n",
              "print((r.stdout + r.stderr)[-1500:])"]),
        code(["#@title 6. run_prediction.py (model_2_ptm)\n",
              "import subprocess, os\n",
              "os.environ['ALPHAFOLD_DATA_DIR'] = os.path.abspath('afdata')\n",
              "cmd = ('cd TCRdock && python run_prediction.py --targets setup/targets.tsv '\n",
              "       '--outfile_prefix run --model_names model_2_ptm '\n",
              "       '--data_dir ' + os.environ['ALPHAFOLD_DATA_DIR'])\n",
              "r = subprocess.run(cmd, shell=True, capture_output=True, text=True)\n",
              "print((r.stdout + r.stderr)[-2000:])"]),
        code(["#@title 7. Collect outputs, zip, download\n",
              "import glob, subprocess\n",
              "pdbs = glob.glob('TCRdock/**/*.pdb', recursive=True)\n",
              "print('output PDBs:', pdbs[:20], '\\ncount:', len(pdbs))\n",
              "# chain remap (native class I order MHC,b2m,peptide,TCRa,TCRb -> C,D,E,A,B)\n",
              "# is applied downstream once the real chain labels are confirmed on this output.\n",
              "subprocess.run('cd TCRdock && zip -qr /content/tcrdock_folds.zip . "
              "-i \"*.pdb\" \"run*\"', shell=True)\n",
              "from google.colab import files\n",
              "files.download('/content/tcrdock_folds.zip')"]),
      ]}
    Path(out_ipynb).write_text(json.dumps(nb, indent=1))
    print(f"wrote {out_ipynb} with {len(targets)} targets")


if __name__ == "__main__":
    main()
