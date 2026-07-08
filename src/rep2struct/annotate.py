from __future__ import annotations
from .schema import Clonotype, Annotation

# ascending tcrdist thresholds -> tier. Calibrated on the validation arm later.
DEFAULT_TIERS = [(12.0, "high"), (24.0, "medium"), (48.0, "low")]

def _tier(distance, tiers):
    for thr, name in tiers:
        if distance <= thr:
            return name
    return "unannotatable"

def _default_sim(cdr3_a, v_a, cdr3_b, v_b, species="human", top_k=5):
    from tcr_explorer.similarity import find_similar_paired_tcrs
    return find_similar_paired_tcrs(cdr3_a, v_a, cdr3_b, v_b, species=species, top_k=top_k)

def annotate(clonotypes, sim_fn=None, tiers=DEFAULT_TIERS):
    fn = sim_fn or _default_sim
    out = []
    for c in clonotypes:
        neigh, *_ = fn(c.cdr3a, c.trav, c.cdr3b, c.trbv, species="human", top_k=5)
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
                hla=best.get("mhc"), antigen=best.get("antigen")))
    return out
