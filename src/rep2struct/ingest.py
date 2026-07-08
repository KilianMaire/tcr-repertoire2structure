from __future__ import annotations
import hashlib
from collections import defaultdict
import pandas as pd
from .schema import Clonotype

def _clon_id(trav, cdr3a, trbv, cdr3b) -> str:
    key = f"{trav}|{cdr3a}|{trbv}|{cdr3b}"
    return hashlib.sha1(key.encode()).hexdigest()[:12]

def parse_10x(path, report: bool = False):
    df = pd.read_csv(path)
    df = df[(df["productive"].astype(str) == "True") & (df["high_confidence"].astype(str) == "True")]
    per_cell = defaultdict(dict)  # barcode -> {"TRA": row, "TRB": row}
    for _, r in df.iterrows():
        chain = r["chain"]
        if chain in ("TRA", "TRB"):
            per_cell[r["barcode"]][chain] = r
    tuples = defaultdict(set)  # tuple -> set of barcodes
    dropped_unpaired = 0
    for bc, chains in per_cell.items():
        if "TRA" not in chains or "TRB" not in chains:
            dropped_unpaired += 1
            continue
        a, b = chains["TRA"], chains["TRB"]
        key = (a["v_gene"], a["cdr3"], b["v_gene"], b["cdr3"])
        tuples[key].add(bc)
    clons = [
        Clonotype(id=_clon_id(*k), trav=k[0], cdr3a=k[1], trbv=k[2], cdr3b=k[3], size=len(bcs))
        for k, bcs in tuples.items()
    ]
    clons.sort(key=lambda c: c.size, reverse=True)
    if report:
        return clons, {"dropped_unpaired": dropped_unpaired, "clonotypes": len(clons)}
    return clons
