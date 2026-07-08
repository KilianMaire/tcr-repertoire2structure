from __future__ import annotations

def _marker(run_state, cid):
    return run_state.run_dir / f"fold_{cid}.done.txt"

def run_folds(jobs, fold_fn, run_state):
    out = []
    for job in jobs:
        marker = _marker(run_state, job.clonotype_id)
        if marker.exists():
            job.status = "done"
            job.model_paths = marker.read_text().splitlines()
            out.append(job)
            continue
        try:
            paths = fold_fn(job)
            job.model_paths = list(paths)
            job.status = "done"
            marker.write_text("\n".join(job.model_paths))
        except Exception as exc:  # noqa: BLE001
            job.status = "failed"
            job.model_paths = []
        out.append(job)
    return out
