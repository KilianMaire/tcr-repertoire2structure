"""Validation arm: measure the annotation step against dextramer ground truth.

The 10x 4-donor CD8 dextramer set gives, per cell, a binarized specificity call
(the `_binder` columns; the peptide is embedded in the column name). We collapse
cells to clonotypes, assign each clonotype the dextramer epitope its cells agree
on, then run the real TCRdist annotation and score predicted vs true epitope.

Leakage guard: the 10x dextramer study is itself a VDJdb source, so a clonotype's
identical TCR may sit in the reference (distance ~0), trivially inflating
accuracy. We report metrics both raw and after dropping near-zero-distance
(leakage-suspected) matches, and print the nearest-distance distribution.

Usage:
  python scripts/run_validation_arm.py <dextramer_dir> <out_json> [donors] [max_clonos]
  # dextramer_dir holds donor{N}_all_contig_annotations.csv + donor{N}_binarized_matrix.csv
"""
from __future__ import annotations
import json
import sys
from collections import defaultdict, Counter
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from rep2struct.ingest import _clon_id, standardize_alleles
from rep2struct.annotate import annotate, _default_sim, _tier, DEFAULT_TIERS
from rep2struct.schema import Clonotype, Annotation
from rep2struct.validate import annotation_metrics


def _binder_cols(bm):
    return [c for c in bm.columns if c.endswith("_binder") and c[0] in "AB" and c[1:5].isdigit()]


def _epitope_of(col):
    # 'A0201_GILGFVFTL_Flu-MP_Influenza_binder' -> ('GILGFVFTL', 'HLA-A*02:01')
    parts = col.split("_")
    pep = parts[1]
    hla = f"HLA-{parts[0][0]}*{parts[0][1:3]}:{parts[0][3:5]}"
    return pep, hla


def cell_labels(bm):
    """barcode -> (epitope, hla) when exactly one dextramer binds (NC excluded)."""
    cols = [c for c in _binder_cols(bm) if "_NC_" not in c]
    sub = bm[["barcode"] + cols].fillna(0)
    labels = {}
    mat = sub[cols].values
    for i, bc in enumerate(sub["barcode"].values):
        pos = [j for j in range(len(cols)) if mat[i, j] == 1]
        if len(pos) == 1:
            labels[bc] = _epitope_of(cols[pos[0]])
    return labels


def labeled_clonotypes(contig_csv, bm_csv):
    df = pd.read_csv(contig_csv)
    df = df[(df["productive"].astype(str) == "True") &
            (df["high_confidence"].astype(str) == "True") &
            (df["chain"].isin(["TRA", "TRB"]))]
    per_cell = defaultdict(dict)
    for r in df.itertuples():
        per_cell[r.barcode][r.chain] = r
    bm = pd.read_csv(bm_csv)
    cl = cell_labels(bm)

    tuples = defaultdict(set)
    jgenes = defaultdict(lambda: {"A": [], "B": []})
    votes = defaultdict(Counter)
    for bc, ch in per_cell.items():
        if "TRA" not in ch or "TRB" not in ch:
            continue
        a, b = ch["TRA"], ch["TRB"]
        key = (a.v_gene, a.cdr3, b.v_gene, b.cdr3)
        tuples[key].add(bc)
        jgenes[key]["A"].append(a.j_gene)
        jgenes[key]["B"].append(b.j_gene)
        if bc in cl:
            votes[key][cl[bc]] += 1

    clons, labels, hlas = [], {}, {}
    for key, bcs in tuples.items():
        cid = _clon_id(*key)
        clons.append(Clonotype(id=cid, trav=key[0], cdr3a=key[1], trbv=key[2],
                               cdr3b=key[3], size=len(bcs),
                               traj=Counter(jgenes[key]["A"]).most_common(1)[0][0],
                               trbj=Counter(jgenes[key]["B"]).most_common(1)[0][0]))
        v = votes[key]
        if v:
            (ep, hla), n = v.most_common(1)[0]
            if n / sum(v.values()) >= 0.5:  # dominant dextramer
                labels[cid] = ep
                hlas[cid] = hla
    return clons, labels, hlas


def nearest_cache(clons, top_k=1):
    """One similarity query per clonotype: cid -> nearest neighbour dict (or None).
    Everything downstream (annotation, distance, sweep) reads this cache."""
    cache = {}
    for i, c in enumerate(clons):
        neigh, *_ = _default_sim(c.cdr3a, c.trav, c.cdr3b, c.trbv, top_k=top_k)
        cache[c.id] = neigh[0] if neigh else None
        if (i + 1) % 100 == 0:
            print(f"  ...{i + 1}/{len(clons)} queried", flush=True)
    return cache


def annotations_from_cache(clons, cache, tiers=DEFAULT_TIERS):
    anns = []
    for c in clons:
        n = cache.get(c.id)
        if not n:
            anns.append(Annotation(c.id, annotatable=False, confidence_tier="unannotatable"))
            continue
        tier = _tier(n["distance"], tiers)
        if tier == "unannotatable":
            anns.append(Annotation(c.id, annotatable=False, confidence_tier="unannotatable",
                                   tcrdist=n["distance"]))
        else:
            anns.append(Annotation(c.id, annotatable=True, confidence_tier=tier,
                                   tcrdist=n["distance"], epitope=n.get("epitope"),
                                   hla=n.get("mhc"), antigen=n.get("antigen")))
    return anns


def main():
    dex_dir = Path(sys.argv[1])
    out_json = sys.argv[2]
    donors = [int(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 else [1, 2, 3, 4]
    max_clonos = int(sys.argv[4]) if len(sys.argv) > 4 else 0

    clons, labels, hlas = [], {}, {}
    for d in donors:
        c, l, h = labeled_clonotypes(
            dex_dir / f"donor{d}_all_contig_annotations.csv",
            dex_dir / f"donor{d}_binarized_matrix.csv")
        # dedup across donors by clonotype id
        seen = {x.id for x in clons}
        clons += [x for x in c if x.id not in seen]
        labels.update(l); hlas.update(h)
        print(f"donor{d}: {len(c)} clonotypes, {len(l)} dextramer-labeled")

    # score only clonotypes that carry a ground-truth label (the recall universe)
    labeled = [c for c in clons if c.id in labels]
    if max_clonos and len(labeled) > max_clonos:
        labeled = sorted(labeled, key=lambda c: c.size, reverse=True)[:max_clonos]
    print(f"total labeled clonotypes to score: {len(labeled)}")

    labeled = standardize_alleles(labeled)
    cache = nearest_cache(labeled, top_k=1)
    dist = {cid: (n["distance"] if n else float("inf")) for cid, n in cache.items()}
    anns = annotations_from_cache(labeled, cache)

    raw = annotation_metrics(anns, labels)
    # leakage-suspected = an exact/near-identical reference match
    LEAK = 1.0
    leaked = {cid for cid, dd in dist.items() if dd <= LEAK}
    nonleak_anns = [a for a in anns if a.clonotype_id not in leaked]
    deleaked = annotation_metrics(nonleak_anns, labels)

    # nearest-distance distribution among correctly annotated
    correct_d = sorted(round(dist[a.clonotype_id], 1) for a in anns
                       if a.annotatable and labels.get(a.clonotype_id) == a.epitope)
    import numpy as np
    pct = {p: float(np.percentile(correct_d, p)) for p in (10, 25, 50, 75, 90)} if correct_d else {}

    # threshold sweep (single cutoff), all from the cache; report de-leaked too
    sweep = []
    for thr in (6, 12, 18, 24, 36, 48, 60, 90):
        preds = [(c.id, cache[c.id]["epitope"]) for c in labeled
                 if cache[c.id] and cache[c.id]["distance"] <= thr]
        corr = sum(1 for cid, ep in preds if labels.get(cid) == ep)
        corr_nl = sum(1 for cid, ep in preds if labels.get(cid) == ep and cid not in leaked)
        n_nl = sum(1 for cid, ep in preds if cid not in leaked)
        sweep.append({"thr": thr, "n_pred": len(preds),
                      "precision": corr / len(preds) if preds else None,
                      "recall": corr / len(labeled) if labeled else None,
                      "precision_deleaked": corr_nl / n_nl if n_nl else None,
                      "recall_deleaked": corr_nl / (len(labeled) - len(leaked)) if labeled else None})

    result = {"donors": donors, "n_labeled_clonotypes": len(labeled),
              "n_leakage_suspected": len(leaked & {c.id for c in labeled}),
              "metrics_raw": raw, "metrics_deleaked": deleaked,
              "correct_distance_percentiles": pct, "threshold_sweep": sweep}
    Path(out_json).write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
