"""Build the Protenix fold notebook for a benchmark seed batch (A100 runtime).

Reads every construct FASTA emitted by run_benchmark_arm.py (cognate + decoys +
scramble, keyed {clonotype}__{epitope}) and embeds them into the validated
MSA-based Protenix notebook (tools.notebook.build_notebook), the same recipe the
Protenix agent runs (live-validated on an A100: flu M1 pLDDT 46->95). The MSA is
deduped by unique protein chain, so the shared TCR/MHC chains are searched once
and reused across a TCR's constructs.

On top of the validated recipe this adds two robustness cells for an unattended
batch: it drops a truncated Protenix checkpoint before folding, and the fold loop
self-heals the 'failed finding central directory' zip error (a truncated weights
download) by deleting /root/checkpoint/*.pt and retrying the key once.

Target runtime: A100 (sm_80, stock torch works; do NOT use a Blackwell GPU, whose
sm_120 needs a torch cu128 rebuild). The user runs the cells and returns out/.

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


def _code(*lines):
    return {"cell_type": "code", "metadata": {}, "execution_count": None,
            "outputs": [], "source": list(lines)}


_CKPT_GUARD = _code(
    "# 2b. drop any TRUNCATED Protenix checkpoint (<500MB) so the fold re-downloads it whole.\n",
    "# guards the 'failed finding central directory' zip error from a partial weights download.\n",
    "import glob, os\n",
    "for f in glob.glob('/root/checkpoint/*.pt'):\n",
    "    mb = os.path.getsize(f) / 1e6\n",
    "    if mb < 500:\n",
    "        print('removing truncated checkpoint', f, f'{mb:.1f}MB'); os.remove(f)\n",
    "    else:\n",
    "        print('checkpoint ok', f, f'{mb:.0f}MB')\n",
)

_HARDENED_FOLD = _code(
    "# 3. Fold each record with the provided MSA (unpairedMsaPath in the JSON).\n",
    "#\n",
    "# Output STREAMS live (Popen, line by line) and is teed to out/{key}.log, so a\n",
    "# stall or a Protenix traceback is visible immediately instead of being swallowed\n",
    "# by capture_output. A watchdog kills a fold that hangs even when it prints NOTHING\n",
    "# (the silent-hang case: a plain wait(timeout=) never fires because the stdout read\n",
    "# loop blocks first). Self-heals a truncated checkpoint on the 'central directory'\n",
    "# / PytorchStreamReader zip error by deleting /root/checkpoint/*.pt and retrying.\n",
    "import glob, os, subprocess, json, signal, threading\n",
    "\n",
    "SMOKE = True        # True: fold ONLY the first construct as a cheap live test, then\n",
    "                    # flip to False and re-run to fold the whole batch.\n",
    "TIMEOUT_S = 1800    # kill a single fold after this many seconds (backstop for hangs).\n",
    "\n",
    "def _run_fold(key):\n",
    "    # --use_msa true is REQUIRED so Protenix consumes the provided unpairedMsaPath\n",
    "    # and does NOT fall back to an online MSA server search. Omitting it makes the\n",
    "    # fold hang in a server-poll loop (loads weights to GPU, then nanosleep with the\n",
    "    # GPU idle). This matches the hdm-tcr-derp production recipe, which asserts the\n",
    "    # log line 'do not need to update msa result' to confirm the server is untouched.\n",
    "    cmd = (f'protenix pred -i inputs/{key}.json -o out/{key} -s 101 '\n",
    "           '-n protenix_base_default_v1.0.0 --use_msa true --use_default_params true')\n",
    "    os.makedirs('out', exist_ok=True)\n",
    "    log = f'out/{key}.log'\n",
    "    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,\n",
    "                         stderr=subprocess.STDOUT, text=True, bufsize=1,\n",
    "                         start_new_session=True)  # own group so killpg reaps children\n",
    "    timed_out = {'v': False}\n",
    "    def _kill():\n",
    "        timed_out['v'] = True\n",
    "        try: os.killpg(os.getpgid(p.pid), signal.SIGKILL)\n",
    "        except ProcessLookupError: pass\n",
    "    watchdog = threading.Timer(TIMEOUT_S, _kill); watchdog.start()\n",
    "    lines = []\n",
    "    with open(log, 'w') as lf:\n",
    "        for line in p.stdout:              # blocks per line; EOF when the process dies\n",
    "            print(line, end='', flush=True)\n",
    "            lf.write(line); lf.flush(); lines.append(line)\n",
    "    p.wait(); watchdog.cancel()\n",
    "    if timed_out['v']:\n",
    "        print(f'\\nTIMEOUT after {TIMEOUT_S}s -> killed', key, flush=True)\n",
    "        return 124, ''.join(lines)\n",
    "    return p.returncode, ''.join(lines)\n",
    "\n",
    "os.makedirs('out', exist_ok=True)\n",
    "keys = sorted(INPUTS)\n",
    "if SMOKE:\n",
    "    keys = keys[:1]\n",
    "    print('SMOKE mode: folding only', keys[0], '(set SMOKE=False for the full batch)', flush=True)\n",
    "manifest = {}\n",
    "for key in keys:\n",
    "    if glob.glob(f'out/{key}/**/*.cif', recursive=True):\n",
    "        manifest[key] = sorted(glob.glob(f'out/{key}/**/*.cif', recursive=True))\n",
    "        print('SKIP (done)', key, flush=True); continue\n",
    "    print('RUN', key, flush=True)\n",
    "    rc, out = _run_fold(key)\n",
    "    if rc != 0 and ('central directory' in out or 'PytorchStreamReader' in out):\n",
    "        print('CHECKPOINT CORRUPT -> delete + re-download + retry', key, flush=True)\n",
    "        for f in glob.glob('/root/checkpoint/*.pt'):\n",
    "            os.remove(f)\n",
    "        rc, out = _run_fold(key)\n",
    "    if rc != 0:\n",
    "        print('FAIL', key, 'rc', rc, '(see out/%s.log)' % key, flush=True); continue\n",
    "    manifest[key] = sorted(glob.glob(f'out/{key}/**/*.cif', recursive=True))\n",
    "    print('FOLDED', key, len(manifest[key]), 'models', flush=True)\n",
    "json.dump(manifest, open('/content/protenix_result.json', 'w'), indent=2)\n",
    "print('DONE', {k: len(v) for k, v in manifest.items()})\n",
)


def _harden(nb):
    """Insert the checkpoint guard before the fold cell and replace the fold cell
    with the self-healing version. The fold cell is the one that runs `protenix pred`."""
    cells = nb["cells"]
    # match the actual CLI invocation, not the substring "protenix pred" which also
    # appears inside the INPUTS cell comment ("protenix prediction JSON").
    fold_idx = next(i for i, c in enumerate(cells)
                    if "protenix pred -i" in "".join(c.get("source", [])))
    cells[fold_idx] = _HARDENED_FOLD
    cells.insert(fold_idx, _CKPT_GUARD)
    return nb


def build_inputs(constructs_dir):
    inputs = {}
    for fasta in sorted(Path(constructs_dir).glob("*.fasta")):
        chains = parse_fasta(fasta.read_text())
        inputs[fasta.stem] = _to_protenix(fasta.stem, chains)
    return inputs


def main():
    constructs_dir, out_ipynb = sys.argv[1], sys.argv[2]
    inputs = build_inputs(constructs_dir)
    nb = _harden(build_notebook("protenix", inputs))
    Path(out_ipynb).write_text(json.dumps(nb, indent=1))
    print(f"embedded {len(inputs)} constructs -> {out_ipynb}")
    print("keys:", sorted(inputs))


if __name__ == "__main__":
    main()
