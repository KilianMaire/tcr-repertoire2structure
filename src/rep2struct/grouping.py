from __future__ import annotations
from .schema import FoldJob


def group_key(job: FoldJob) -> str:
    tcr = "tcr" if job.has_tcr else "notcr"
    return f"c{job.mhc_class}_{tcr}_{job.species}_{job.output_needed}"


def partition(jobs: list[FoldJob]) -> dict[str, list[FoldJob]]:
    groups: dict[str, list[FoldJob]] = {}
    for j in jobs:
        k = group_key(j)
        j.group_id = k
        groups.setdefault(k, []).append(j)
    return groups
