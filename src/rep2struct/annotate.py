from __future__ import annotations
from .schema import Clonotype, Annotation
from .tools.construct_io import normalize_hla

# ascending tcrdist thresholds -> tier. Calibrated on the validation arm later.
DEFAULT_TIERS = [(12.0, "high"), (24.0, "medium"), (48.0, "low")]

def _tier(distance, tiers):
    for thr, name in tiers:
        if distance <= thr:
            return name
    return "unannotatable"

def _default_sim(cdr3_a, v_a, cdr3_b, v_b, species="human", top_k=5):
    from tcr_explorer.similarity import find_similar_paired_tcrs
    neigh, engine, total, warns = find_similar_paired_tcrs(
        cdr3_a, v_a, cdr3_b, v_b, species=species, top_k=top_k)
    # find_similar_paired_tcrs returns pydantic PairedNeighbour objects with
    # fields (distance, epitope_aa, mhc_a, antigen, ...). annotate consumes a
    # dict keyed epitope/mhc/antigen/distance, so map them here. Keeping the
    # dict contract lets the offline fakes stay plain dicts.
    dicts = [{"distance": n.distance, "epitope": n.epitope_aa,
              "mhc": n.mhc_a, "antigen": n.antigen} for n in neigh]
    return dicts, engine, total, warns

def _cache_key(c, species="human"):
    return "|".join([species, c.trav or "", c.cdr3a or "", c.trbv or "", c.cdr3b or ""])


def _neighbours(clonotypes, fn, cache_path, max_workers):
    """Return {index: neighbour-list} for every clonotype, keyed by index so the caller
    reassembles in input order (determinism preserved regardless of completion order).

    The durable scalability win is the ON-DISK CACHE keyed by the (species, V/CDR3) tuple: a
    rerun does zero lookups, and the cache is flushed incrementally so a run that dies at
    clonotype 2999/3000 stays resumable. The thread pool only helps when fn is I/O-bound (the
    hosted tcr_explorer service) or releases the GIL; the bundled find_similar_paired_tcrs is a
    local CPU-bound BLOSUM/pandas scan, so threads give little speedup there (set
    R2S_ANNOTATE_WORKERS=1 to skip the overhead). The cap upstream bounds the total work."""
    import concurrent.futures
    import json
    from pathlib import Path

    def _flush(cache):
        if not cache_path:
            return
        try:
            Path(cache_path).write_text(json.dumps(cache))
        except (TypeError, ValueError, OSError):
            pass  # a serialization slip must never crash annotate after the expensive work

    cache = {}
    if cache_path and Path(cache_path).exists():
        try:
            cache = json.loads(Path(cache_path).read_text())
        except (ValueError, OSError):
            cache = {}
    results, todo = {}, []
    for i, c in enumerate(clonotypes):
        key = _cache_key(c)
        if key in cache:
            results[i] = cache[key]
        else:
            todo.append((i, c))
    if todo:
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as ex:
            fut = {ex.submit(fn, c.cdr3a, c.trav, c.cdr3b, c.trbv, "human", 5): (i, c)
                   for i, c in todo}
            for n, f in enumerate(concurrent.futures.as_completed(fut), 1):
                i, c = fut[f]
                neigh = f.result()[0]     # a raising fn propagates loudly, never caches a false empty
                results[i] = neigh
                cache[_cache_key(c)] = neigh
                if n % 256 == 0:          # incremental flush -> partial runs stay resumable
                    _flush(cache)
        _flush(cache)
    return results


def annotate(clonotypes, sim_fn=None, tiers=DEFAULT_TIERS, cache_path=None, max_workers=8):
    fn = sim_fn or _default_sim
    neigh_by_i = _neighbours(clonotypes, fn, cache_path, max_workers)
    out = []
    for i, c in enumerate(clonotypes):
        neigh = neigh_by_i.get(i) or []
        if not neigh:
            out.append(Annotation(clonotype_id=c.id, annotatable=False,
                                  confidence_tier="unannotatable"))
            continue
        best = min(neigh, key=lambda n: n["distance"])
        tier = _tier(best["distance"], tiers)
        if tier == "unannotatable":
            out.append(Annotation(clonotype_id=c.id, annotatable=False,
                                  confidence_tier="unannotatable", tcrdist=best["distance"]))
        else:
            out.append(Annotation(
                clonotype_id=c.id, annotatable=True, confidence_tier=tier,
                tcrdist=best["distance"], epitope=best.get("epitope"),
                hla=normalize_hla(best.get("mhc")), antigen=best.get("antigen")))
    return out
