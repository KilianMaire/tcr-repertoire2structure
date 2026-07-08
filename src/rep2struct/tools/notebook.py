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
    """Validated MHC-Fine adapter (live re-validated 2026-07-09, end-to-end pose).

    INPUTS is {key: {protein_sequence, peptide_sequence}} (e.g. cognate + scramble).

    numpy lesson (learned the hard way, live): do NOT downgrade numpy. On today's Colab
    image `pip install "numpy<2"` poisons numpy's OWN compiled mtrand.so (a dtype-size
    ABI wall, 'expected 96 got 88') that survives force-reinstall and restarts. The robust
    recipe KEEPS stock numpy 2 and shims the AF2-era removals the code relies on
    (np.string_ / np.unicode_ / np.float_ / np.complex_ / np.bool8, plus np.sum over a
    generator, which numpy 2 turned into a hard error) BEFORE `from src import ...`. No
    downgrade, no ABI wall, no kernel restart. kalign is still required by preprocess
    (else kalign.py joins a None path -> TypeError). The bundled msa_run binary builds the
    MSA (no local DB). The output pose is chain A (MHC heavy) + chain B (peptide), no
    b2m/TCR; a cognate 9-mer lands at ~97 mean pLDDT."""
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
            _code("# 2. deps: keep STOCK numpy 2 (downgrading poisons numpy's mtrand.so). kalign is\n",
                  "# required by preprocess (else kalign.py joins a None path -> TypeError).\n",
                  "import subprocess\n",
                  "subprocess.run('pip install -q gdown biopython ml_collections dm-tree einops',\n",
                  "               shell=True, check=True)\n",
                  "subprocess.run('apt-get -qq install -y kalign', shell=True, check=True)\n"),
            _code("# 3. numpy-2 compat shim BEFORE importing src: restore the AF2-era aliases the\n",
                  "# code uses and make np.sum tolerate a generator (numpy 2 made that a hard error).\n",
                  "import types, numpy as np\n",
                  "np.string_ = np.bytes_\n",
                  "np.unicode_ = np.str_\n",
                  "np.float_ = np.float64\n",
                  "np.complex_ = np.complex128\n",
                  "np.bool8 = np.bool_\n",
                  "_orig_sum = np.sum\n",
                  "def _sum(a, *ar, **kw):\n",
                  "    return _orig_sum(list(a) if isinstance(a, types.GeneratorType) else a, *ar, **kw)\n",
                  "np.sum = _sum\n",
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
