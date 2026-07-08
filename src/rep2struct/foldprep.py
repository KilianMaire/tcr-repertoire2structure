from __future__ import annotations
from .schema import FoldJob

TIER_WEIGHT = {"high": 3.0, "medium": 2.0, "low": 1.0, "unannotatable": 0.0}

def select_top(clonotypes, annotations, n):
    by_id = {a.clonotype_id: a for a in annotations}
    scored = []
    for c in clonotypes:
        a = by_id.get(c.id)
        if a is None:
            continue
        score = TIER_WEIGHT.get(a.confidence_tier, 0.0) * c.size
        scored.append((score, c, a))
    scored.sort(key=lambda t: t[0], reverse=True)
    return [(c, a) for _, c, a in scored[:n]]

def build_construct(clonotype, annotation, tcr_seqs, mhc_seqs) -> FoldJob:
    t = tcr_seqs[clonotype.id]
    m = mhc_seqs[annotation.hla]
    fasta = "\n".join([
        ">A", t["A"],            # TCR alpha V domain
        ">B", t["B"],            # TCR beta V domain
        ">C", m["heavy"],        # MHC class I heavy chain
        ">D", m["b2m"],          # beta 2 microglobulin
        ">E", annotation.epitope,  # peptide
    ])
    return FoldJob(clonotype_id=clonotype.id, construct_fasta=fasta)
