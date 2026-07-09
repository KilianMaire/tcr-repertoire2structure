from __future__ import annotations
from collections import defaultdict
import random as _random
from .schema import Annotation
from .foldprep import build_construct

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
