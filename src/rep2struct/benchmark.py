from __future__ import annotations
from collections import defaultdict
import random as _random
from .schema import Annotation
from .foldprep import build_construct
from .qc import ensemble_contact

def is_novel(tcrdist, leak_thr: float = 1.0) -> bool:
    return tcrdist is None or tcrdist > leak_thr

def panel_epitopes(truth):
    return sorted({(ep, hla) for (ep, hla) in truth.values()})

def per_hla_novel_counts(clonotypes, truth, annotations):
    dist = {a.clonotype_id: getattr(a, "tcrdist", None) for a in annotations}
    out = {}
    for c in clonotypes:
        if c.id not in truth:
            continue
        ep, hla = truth[c.id]
        novel = is_novel(dist.get(c.id))
        h = out.setdefault(hla, {"n_total": 0, "n_novel": 0, "epitopes": {}})
        h["n_total"] += 1
        h["n_novel"] += int(novel)
        e = h["epitopes"].setdefault(ep, {"n": 0, "n_novel": 0})
        e["n"] += 1
        e["n_novel"] += int(novel)
    return out

def decoys_for(cognate, hla, panel, k):
    # same-HLA ONLY: cross-HLA decoys reintroduce the HLA-geometry confound
    same = sorted(ep for (ep, h) in panel if h == hla and ep != cognate)
    return same[:k]

def scramble_peptide(cognate, seed=0):
    # composition-preserving shuffle; deterministic; retry so it differs from cognate
    chars = list(cognate)
    rng = _random.Random(f"{cognate}:{seed}")
    for _ in range(20):
        rng.shuffle(chars)
        s = "".join(chars)
        if s != cognate:
            return s
    return "".join(chars)

def build_panel_constructs(clonotype, cognate, hla, decoys, tcr_seqs, mhc_seqs, scramble_seed=0):
    jobs = {}
    peptides = {ep: ep for ep in [cognate, *decoys]}
    peptides["__scramble__"] = scramble_peptide(cognate, scramble_seed)
    for key, pep in peptides.items():
        ann = Annotation(clonotype_id=clonotype.id, annotatable=True,
                         confidence_tier="benchmark", epitope=pep, hla=hla)
        jobs[key] = build_construct(clonotype, ann, tcr_seqs, mhc_seqs)
    return jobs

def contact_by_epitope(paths_by_epitope):
    return {ep: ensemble_contact(paths)[0] for ep, paths in paths_by_epitope.items()}

def retrieval_result(contacts, cognate):
    # exclude the scramble key from retrieval; it is a separate contrast
    scored = {e: v for e, v in contacts.items() if e != "__scramble__"}
    cval = scored.get(cognate)
    valid = [v for v in scored.values() if v is not None]
    ranked = sorted(scored, key=lambda e: (-1.0 if scored[e] is None else scored[e]),
                    reverse=True)
    if cval is None or not valid:
        top1 = False
    else:
        best = max(valid)
        n_at_best = sum(1 for v in valid if v == best)
        top1 = (cval == best) and (n_at_best == 1)
    return {"ranked": ranked, "top1": top1, "cognate_contact": cval}

def auroc(pairs):
    num = den = 0.0
    for cog, decoys in pairs:
        for d in decoys:
            den += 1
            num += 1.0 if cog > d else (0.5 if cog == d else 0.0)
    return None if den == 0 else num / den
