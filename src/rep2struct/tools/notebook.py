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


def _affinetune_notebook(inputs: dict) -> dict:
    """Validated affinetune (phbradley/alphafold_finetune) adapter (calibrated live
    2026-07-09 on Colab T4: cognate GILGFVFTL pae 1.02 vs scramble TFVFGLIGL pae 2.25).

    INPUTS is {key: {mhc, b2m, peptide}} (cognate + scramble, keys prefixed by clonotype
    id). affinetune returns a PRESENTATION score, not a structure, so each record yields a
    one-float score file that qc reads via verdict_binding.

    Env lesson (shared with tcrdock): the JAX/AF wheels are cp38, uninstallable in Colab's
    py3.12 kernel. The recipe builds a side conda py3.8 env and SHELLS OUT to it (no kernel
    restart). Two CUDA gotchas: install cudatoolkit=11.1 + cudnn into the env (jaxlib-cuda111
    ships without the CUDA 11 runtime, else run_prediction aborts on libcudart.so.11.0), and
    every shell-out MUST set LD_LIBRARY_PATH=/opt/conda/envs/af/lib:/usr/lib64-nvidia (prepend,
    keep the driver dir, else cuInit 303 -> silent CPU fallback).

    Two model facts, both confirmed live: the score column is model_2_ptm_ft_pae and LOWER =
    presented, so the reader INVERTS (score = -pae) before verdict_binding, which treats
    HIGHER as more presented. And the pMHC target_chainseq is the alpha1+alpha2 groove domain
    (~175 aa) + peptide, no b2m. Chain C from the construct is the full alpha1-2-3 ectodomain,
    so the adapter takes its first 175 residues (the real allele's own groove) and threads it
    onto a shipped class I template via the length-matched alignment file with
    --ignore_identities. run_prediction asserts target_len == len(query), the fail-loud net if
    a peptide length has no wired alignment (only 9-mers are calibrated)."""
    return {
        "nbformat": 4, "nbformat_minor": 5, "metadata": {"accelerator": "GPU"},
        "cells": [
            _code("# affinetune inputs (embedded): {key: {mhc, b2m, peptide}}\n",
                  "INPUTS = ", repr(inputs)),
            _code("# 1. side conda py3.8 env + alphafold_finetune + its cp38 wheels (no kernel restart).\n",
                  "#    cudatoolkit=11.1 is THE fix: jaxlib-cuda111 ships without the CUDA 11 runtime.\n",
                  "import os, subprocess\n",
                  "def sh(c):\n",
                  "    print('>>>', c, flush=True); subprocess.run(c, shell=True, check=True)\n",
                  "if not os.path.isdir('/opt/conda'):\n",
                  "    sh('wget -qO /tmp/mf.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh && bash /tmp/mf.sh -b -p /opt/conda')\n",
                  "if not os.path.isdir('/opt/conda/envs/af'):\n",
                  "    sh('/opt/conda/bin/mamba create -y -n af python=3.8')\n",
                  "if not os.path.isdir('/content/aff'):\n",
                  "    sh('git clone --branch main https://github.com/phbradley/alphafold_finetune /content/aff')\n",
                  "sh('/opt/conda/envs/af/bin/pip install -q -r /content/aff/requirements_colab_python38_v2.txt')\n",
                  "sh('/opt/conda/bin/mamba install -y -n af -c conda-forge cudatoolkit=11.1 cudnn')\n"),
            _code("# 2. base AF params (model_2_ptm) + the fine-tune dataset bundle (params pkl + alignments)\n",
                  "import os, subprocess\n",
                  "os.makedirs('/content/alphafold_params/params', exist_ok=True)\n",
                  "if not os.path.exists('/content/alphafold_params/params/params_model_2_ptm.npz'):\n",
                  "    subprocess.run('wget -qO /content/alphafold_params/params/params_model_2_ptm.npz https://www.dropbox.com/s/e3uz9mwxkmmv35z/params_model_2_ptm.npz', shell=True, check=True)\n",
                  "if not os.path.isdir('/content/aff/datasets_alphafold_finetune'):\n",
                  "    subprocess.run('cd /content/aff && wget -q https://files.ipd.uw.edu/pub/alphafold_finetune_motmaen_pnas_2023/datasets_alphafold_finetune_v2_2023-02-20.tgz && tar -xzf datasets_alphafold_finetune_v2_2023-02-20.tgz', shell=True, check=True)\n"),
            _code("# 3. score each record: groove alpha (mhc[:175]) + peptide, length-matched shipped\n",
                  "#    alignment, --ignore_identities. Read model_2_ptm_ft_pae and INVERT (lower pae =\n",
                  "#    presented; verdict_binding wants higher = presented). Write one float per key.\n",
                  "import os, subprocess, json, pandas as pd\n",
                  "AFF, PY = '/content/aff', '/opt/conda/envs/af/bin/python'\n",
                  "env = dict(os.environ, LD_LIBRARY_PATH='/opt/conda/envs/af/lib:/usr/lib64-nvidia')\n",
                  "# shipped 9-mer class I alignment (target_len 184 = 175-aa alpha + 9-mer); calibrated live\n",
                  "ALIGN = {9: 'examples/tiny_pmhc_finetune/alignments/A0203_alignments.tsv'}\n",
                  "os.makedirs('/content/output', exist_ok=True)\n",
                  "scores = {}\n",
                  "for key, rec in INPUTS.items():\n",
                  "    pep = rec['peptide']; L = len(pep)\n",
                  "    if L not in ALIGN:\n",
                  "        raise NotImplementedError(f'no calibrated affinetune alignment for a {L}-mer (only 9-mers wired)')\n",
                  "    alpha = rec['mhc'][:175]   # alpha1+alpha2 groove domain of the real allele ectodomain\n",
                  "    tsv = f'{key}_targets.tsv'\n",
                  "    with open(f'{AFF}/{tsv}', 'w') as f:\n",
                  "        f.write('mhc\\tstart\\tpeptide\\ttargetid\\ttarget_chainseq\\ttemplates_alignfile\\n')\n",
                  "        f.write(f'na\\t0\\t{pep}\\t{key}\\t{alpha}/{pep}\\t{ALIGN[L]}\\n')\n",
                  "    cmd = (PY + f' run_prediction.py --targets {tsv} --outfile_prefix {key}'\n",
                  "           ' --model_names model_2_ptm_ft'\n",
                  "           ' --model_params_files datasets_alphafold_finetune/params/mixed_mhc_pae_run6_af_mhc_params_20640.pkl'\n",
                  "           ' --data_dir /content/alphafold_params/ --ignore_identities')\n",
                  "    r = subprocess.run(cmd, shell=True, cwd=AFF, capture_output=True, text=True, env=env)\n",
                  "    assert r.returncode == 0, r.stderr[-2000:]\n",
                  "    pae = float(pd.read_csv(f'{AFF}/{key}_final.tsv', sep='\\t')['model_2_ptm_ft_pae'].iloc[0])\n",
                  "    score = -pae\n",
                  "    sp = f'/content/output/{key}.score'\n",
                  "    open(sp, 'w').write(f'{score:.6f}')\n",
                  "    scores[key] = {'pae': pae, 'score': score, 'score_path': sp}\n",
                  "    print('SCORED', key, 'pae', round(pae, 3), '-> score', round(score, 3), flush=True)\n",
                  "json.dump(scores, open('/content/affinetune_result.json', 'w'), indent=2)\n",
                  "print('DONE', len(scores))\n"),
        ],
    }


def _tcrdock_notebook(inputs: dict) -> dict:
    """Validated TCRdock (phbradley/TCRdock) adapter (calibrated live 2026-07-09 on Colab
    A100: cognate GILGFVFTL peptide<->TCR pae ~11.2 vs scramble TFVFGLIGL ~20.8, peptide
    pLDDT 86 vs 65).

    INPUTS is {key: {row: {10-col tcrdock target}}} (cognate + scramble, keys prefixed by
    clonotype id). TCRdock builds the TCR:pMHC from gene names + CDR3 loops via its own
    templates, so a row is organism/mhc_class/mhc/peptide/va/ja/cdr3a/vb/jb/cdr3b, NOT chain
    sequences. TCRdock returns a structure, but recognition is judged by the peptide<->TCR
    interface PAE, so each record also yields a one-float score file that qc reads via
    verdict_binding (qc_metric=binding_score for tcrdock).

    Env lesson (the whole point; TCRdock declares almost none of it): same side conda py3.8 +
    shell-out pattern as affinetune. TCRdock's requirements.txt is only biopython/numpy/pandas/
    scipy/matplotlib; the ENTIRE AlphaFold stack is deferred to "the AlphaFold README". Its
    bundled AF fork is the 2.3.x line (it annotates jax.Array), so affinetune's AF-2.0 stack is
    the WRONG one. Install DeepMind's AlphaFold v2.3.2 python stack (tensorflow-cpu 2.11, dm-haiku
    0.0.9, numpy 1.21.6, biopython 1.79) + jaxlib 0.3.25 for CUDA 11 / cuDNN 8.05 (matches
    cudatoolkit=11.1). biopython 1.79 is required: TCRdock imports Bio.SubsMat, removed in >=1.80.
    Every shell-out sets LD_LIBRARY_PATH=/opt/conda/envs/tcrdock/lib:/usr/lib64-nvidia.

    Output facts, confirmed live: the PDB is a SINGLE merged chain (MHC groove + peptide + Valpha
    + Vbeta), not A..E. The block layout is recoverable from target_chainseq in {key}_final.tsv,
    four '/'-joined segments in order 0=MHC, 1=peptide, 2=TCRalpha, 3=TCRbeta. The score is the
    peptide<->TCR interface PAE = mean(model_2_ptm_pae_1_2, model_2_ptm_pae_1_3) over the best
    (lowest overall model_2_ptm_pae) row; LOWER pae = more recognized, so the reader INVERTS
    (score = -pae) before verdict_binding, which treats HIGHER as more recognized."""
    return {
        "nbformat": 4, "nbformat_minor": 5, "metadata": {"accelerator": "GPU"},
        "cells": [
            _code("# tcrdock inputs (embedded): {key: {row: {10-col tcrdock target}}}\n",
                  "INPUTS = ", repr(inputs)),
            _code("# 1. side conda py3.8 env + TCRdock + its FULL AF stack. TCRdock's requirements.txt\n",
                  "#    omits AlphaFold entirely; its fork is the AF 2.3.x line (uses jax.Array), so\n",
                  "#    install DeepMind's AF v2.3.2 stack + jaxlib cuda111. biopython==1.79 is required\n",
                  "#    (TCRdock imports Bio.SubsMat, gone in >=1.80). cudatoolkit=11.1 for the CUDA 11 runtime.\n",
                  "import os, subprocess\n",
                  "def sh(c):\n",
                  "    print('>>>', c, flush=True); subprocess.run(c, shell=True, check=True)\n",
                  "if not os.path.isdir('/opt/conda'):\n",
                  "    sh('wget -qO /tmp/mf.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh && bash /tmp/mf.sh -b -p /opt/conda')\n",
                  "if not os.path.isdir('/opt/conda/envs/tcrdock'):\n",
                  "    sh('/opt/conda/bin/mamba create -y -n tcrdock python=3.8')\n",
                  "if not os.path.isdir('/content/TCRdock'):\n",
                  "    sh('git clone https://github.com/phbradley/TCRdock /content/TCRdock')\n",
                  "PIP = '/opt/conda/envs/tcrdock/bin/pip'\n",
                  "sh(PIP + ' install -q -r /content/TCRdock/requirements.txt')\n",
                  "sh('cd /content/TCRdock && /opt/conda/envs/tcrdock/bin/python download_blast.py')\n",
                  "sh('/opt/conda/bin/mamba install -y -n tcrdock -c conda-forge cudatoolkit=11.1 cudnn')\n",
                  "sh('wget -qO /content/af232.txt https://raw.githubusercontent.com/google-deepmind/alphafold/v2.3.2/requirements.txt')\n",
                  "sh(PIP + ' install -q -r /content/af232.txt')\n",
                  "sh(PIP + ' install -q jax==0.3.25 jaxlib==0.3.25+cuda11.cudnn805 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html')\n",
                  "print('ENV_COMPLETE', flush=True)\n"),
            _code("# 2. base AF params (model_2_ptm) into {data_dir}/params/ (same dropbox file affinetune uses)\n",
                  "import os, subprocess\n",
                  "os.makedirs('/content/alphafold_params/params', exist_ok=True)\n",
                  "p = '/content/alphafold_params/params/params_model_2_ptm.npz'\n",
                  "if not os.path.exists(p):\n",
                  "    subprocess.run('wget -qO ' + p + ' https://www.dropbox.com/s/e3uz9mwxkmmv35z/params_model_2_ptm.npz', shell=True, check=True)\n",
                  "print('PARAMS_OK', os.path.getsize(p), flush=True)\n"),
            _code("# 3. fold each record: write the 10-col targets.tsv, setup_for_alphafold, run_prediction,\n",
                  "#    then read the peptide<->TCR interface PAE from {key}_final.tsv and INVERT it (lower pae\n",
                  "#    = more recognized; verdict_binding wants higher). blocks: 0=MHC 1=peptide 2=TCRa 3=TCRb.\n",
                  "import os, subprocess, json\n",
                  "TD, PY = '/content/TCRdock', '/opt/conda/envs/tcrdock/bin/python'\n",
                  "env = dict(os.environ, LD_LIBRARY_PATH='/opt/conda/envs/tcrdock/lib:/usr/lib64-nvidia')\n",
                  "COLS = ['organism','mhc_class','mhc','peptide','va','ja','cdr3a','vb','jb','cdr3b']\n",
                  "os.makedirs('/content/output', exist_ok=True)\n",
                  "def iface_pae(final_tsv):\n",
                  "    # peptide<->TCR interface PAE = mean(pae_1_2, pae_1_3) over the best (lowest overall pae) row\n",
                  "    lines = open(final_tsv).read().rstrip().split('\\n')\n",
                  "    hdr = lines[0].split('\\t'); idx = {h: i for i, h in enumerate(hdr)}\n",
                  "    best = None\n",
                  "    for ln in lines[1:]:\n",
                  "        f = ln.split('\\t')\n",
                  "        overall = float(f[idx['model_2_ptm_pae']])\n",
                  "        pae = (float(f[idx['model_2_ptm_pae_1_2']]) + float(f[idx['model_2_ptm_pae_1_3']])) / 2\n",
                  "        if best is None or overall < best[0]:\n",
                  "            best = (overall, pae)\n",
                  "    return best[1]\n",
                  "scores = {}\n",
                  "for key, rec in INPUTS.items():\n",
                  "    row = rec['row']\n",
                  "    tf = f'{key}_targets.tsv'\n",
                  "    with open(f'{TD}/{tf}', 'w') as f:\n",
                  "        f.write('\\t'.join(COLS) + '\\n')\n",
                  "        f.write('\\t'.join(str(row[c]) for c in COLS) + '\\n')\n",
                  "    s = subprocess.run(f'{PY} setup_for_alphafold.py --targets_tsvfile {tf} --output_dir setup_{key}',\n",
                  "                       shell=True, cwd=TD, capture_output=True, text=True, env=env)\n",
                  "    assert s.returncode == 0, s.stderr[-3000:]\n",
                  "    r = subprocess.run(f'{PY} run_prediction.py --targets setup_{key}/targets.tsv'\n",
                  "                       f' --outfile_prefix {key} --model_names model_2_ptm'\n",
                  "                       ' --data_dir /content/alphafold_params/',\n",
                  "                       shell=True, cwd=TD, capture_output=True, text=True, env=env)\n",
                  "    assert r.returncode == 0, r.stderr[-3000:]\n",
                  "    pae = iface_pae(f'{TD}/{key}_final.tsv')\n",
                  "    score = -pae\n",
                  "    sp = f'/content/output/{key}.score'\n",
                  "    open(sp, 'w').write(f'{score:.6f}')\n",
                  "    scores[key] = {'iface_pae': pae, 'score': score, 'score_path': sp}\n",
                  "    print('SCORED', key, 'iface_pae', round(pae, 3), '-> score', round(score, 3), flush=True)\n",
                  "json.dump(scores, open('/content/tcrdock_result.json', 'w'), indent=2)\n",
                  "print('DONE', len(scores))\n"),
        ],
    }


# Tools whose live Colab cell has been validated and wired to a real recipe.
_WIRED = {"mhcfine": _mhcfine_notebook, "affinetune": _affinetune_notebook,
          "tcrdock": _tcrdock_notebook}


def build_notebook(tool: str, inputs: dict) -> dict:
    """Return a self-contained Colab notebook for a fold tool. Validated tools get their
    real recipe; unwired tools get a fail-loud scaffold (never a faked result)."""
    builder = _WIRED.get(tool)
    if builder:
        return builder(inputs)
    return _stub_notebook(tool, inputs)
