from __future__ import annotations

# Google Drive file id for the MHC-Fine weights (data/model/mhc_fine_weights.pt, 388 MB).
_MHCFINE_WEIGHTS_ID = "1gz8uF8DKE0CzyX_WeDGOX7xP69LjpaZT"


def _code(*lines):
    return {"cell_type": "code", "metadata": {}, "execution_count": None,
            "outputs": [], "source": list(lines)}


def _stub_notebook(tool: str, inputs: dict) -> dict:
    """Fail-loud scaffold for a tool whose live Colab cell is not yet validated.
    It embeds the inputs but refuses to fake a result."""
    return {
        "nbformat": 4, "nbformat_minor": 5, "metadata": {},
        "cells": [
            _code(f"# {tool} inputs (embedded, MSA-free at runtime)\n",
                  "INPUTS = ", repr(inputs)),
            _code(f"# TODO(live): invoke {tool} here against INPUTS and write the model/score.\n",
                  "raise NotImplementedError('live cell not yet validated')\n"),
            _code("# TODO(live): save the output (CIF for structure, float for binding) to OUT_PATH\n"),
        ],
    }


def _mhcfine_notebook(inputs: dict) -> dict:
    """Validated MHC-Fine adapter (live 2026-07-08, T4).

    INPUTS is {key: {protein_sequence, peptide_sequence}} (e.g. cognate + scramble).
    Critical ordering learned live: pin numpy<2 and install kalign BEFORE importing
    numpy/torch/src, otherwise the AF2-derived code hits the removed np.string_ and the
    kernel needs a restart. The bundled msa_run binary builds the MSA (no local DB). The
    output pose is chain A (MHC heavy) + chain B (peptide) only, no b2m/TCR."""
    return {
        "nbformat": 4, "nbformat_minor": 5, "metadata": {"accelerator": "GPU"},
        "cells": [
            _code("# mhcfine inputs (embedded): {key: {protein_sequence, peptide_sequence}}\n",
                  "INPUTS = ", repr(inputs)),
            _code("# 1. clone MHC-Fine\n",
                  "import os, subprocess\n",
                  "if not os.path.isdir('/content/mhc-fine'):\n",
                  "    subprocess.run('git clone https://bitbucket.org/abc-group/mhc-fine.git',\n",
                  "                   shell=True, cwd='/content', check=True)\n",
                  "os.chdir('/content/mhc-fine')\n"),
            _code("# 2. deps FIRST, before any numpy/torch import.\n",
                  "# numpy<2: AF2-derived code uses the removed np.string_. kalign: required by\n",
                  "# preprocess (else kalign.py joins a None path -> TypeError). No restart needed.\n",
                  "import subprocess\n",
                  "subprocess.run('pip install -q \"numpy<2\" gdown biopython ml_collections dm-tree einops',\n",
                  "               shell=True, check=True)\n",
                  "subprocess.run('apt-get -qq install -y kalign', shell=True, check=True)\n"),
            _code("# 3. import (numpy 1.26 loads fresh)\n",
                  "import torch\n",
                  "assert torch.cuda.is_available(), 'need a GPU runtime'\n",
                  "from src import preprocess, model\n"),
            _code("# 4. weights (388 MB) + make the msa_run binary executable\n",
                  "import os, gdown\n",
                  "os.makedirs('data/model', exist_ok=True)\n",
                  "mp = 'data/model/mhc_fine_weights.pt'\n",
                  "if not os.path.exists(mp):\n",
                  f"    gdown.download('https://drive.google.com/uc?id={_MHCFINE_WEIGHTS_ID}', mp, quiet=False)\n",
                  "subprocess.run('chmod +x a3m_generation/msa_run', shell=True, check=True)\n"),
            _code("# 5. fold each embedded record; write ./output/{key}.pdb (chain A=MHC, B=peptide)\n",
                  "import os, json\n",
                  "m = model.Model()\n",
                  "results = {}\n",
                  "for key, rec in INPUTS.items():\n",
                  "    prot, pep = rec['protein_sequence'], rec['peptide_sequence']\n",
                  "    a3m = os.path.join(os.getcwd(), 'data', 'msa', key, 'mmseqs', 'aggregated.a3m')\n",
                  "    if not os.path.exists(a3m):\n",
                  "        preprocess.get_a3m(prot, a3m, key)\n",
                  "    ns = preprocess.preprocess_for_inference(prot, pep, a3m)\n",
                  "    met = m.inference(ns, key)\n",
                  "    results[key] = {'metrics': met, 'pdb': f'output/{key}.pdb'}\n",
                  "    print('FOLDED', key, met, flush=True)\n",
                  "json.dump(results, open('/content/mhcfine_result.json', 'w'), indent=2, default=str)\n",
                  "print('DONE', len(results))\n"),
        ],
    }


# Tools whose live Colab cell has been validated and wired to a real recipe.
_WIRED = {"mhcfine": _mhcfine_notebook}


def build_notebook(tool: str, inputs: dict) -> dict:
    """Return a self-contained Colab notebook for a fold tool. Validated tools get their
    real recipe; unwired tools get a fail-loud scaffold (never a faked result)."""
    builder = _WIRED.get(tool)
    if builder:
        return builder(inputs)
    return _stub_notebook(tool, inputs)
