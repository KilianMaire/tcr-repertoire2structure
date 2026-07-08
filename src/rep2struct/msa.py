from __future__ import annotations
from pathlib import Path
from .schema import FoldJob


def build_msa(job: FoldJob, run_dir, local_runner=None, colab_runner=None) -> tuple[str, str]:
    """Compute an MSA artifact OUTSIDE the fold runtime and cache it.

    Tries local mmseqs2 first, then a Colab CPU step, then falls back to
    MSA-free. Returns (msa_ref, basis). Removing the MSA from the fold
    runtime is what kills the remote-MSA-server throttle failure.
    """
    for runner, basis in ((local_runner, "local"), (colab_runner, "colab_cpu")):
        if runner is None:
            continue
        try:
            a3m = runner(job.construct_fasta)
        except Exception:  # noqa: BLE001  -- any runner failure degrades to the next path
            continue
        out = Path(run_dir) / f"msa_{job.clonotype_id}.a3m"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(a3m)
        return str(out), basis
    return "", "none"
