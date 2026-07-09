from __future__ import annotations
from collections import defaultdict
import random as _random
import warnings
import numpy as np
from .schema import Annotation
from .foldprep import build_construct
from .qc import ensemble_contact, mean_confidence

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

def _residue_bfactors(cif_path, chain_id):
    from Bio.PDB import MMCIFParser
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = next(MMCIFParser(QUIET=True).get_structure("x", str(cif_path)).get_models())
    for ch in m:
        if ch.id == chain_id:
            return [float(np.mean([a.get_bfactor() for a in r])) for r in ch]
    return None

def model_cdr3b_plddt(cif_path, chain_b_seq, cdr3b):
    """Mean pLDDT over the CDR3beta residues, located as a substring of chain B."""
    if not chain_b_seq or not cdr3b:
        return None
    start = chain_b_seq.find(cdr3b)
    if start < 0:
        return None
    bfs = _residue_bfactors(cif_path, "B")
    if bfs is None or start + len(cdr3b) > len(bfs):
        return None
    return mean_confidence(bfs[start:start + len(cdr3b)])

def cdr3b_plddt_by_epitope(paths_by_epitope, chain_b_seq, cdr3b):
    out = {}
    for ep, paths in paths_by_epitope.items():
        vals = [v for v in (model_cdr3b_plddt(p, chain_b_seq, cdr3b) for p in paths)
                if v is not None]
        out[ep] = float(np.median(vals)) if vals else None
    return out

def sequence_baseline_top1(annotation_epitope, cognate):
    return annotation_epitope == cognate

def bootstrap_ci(hits, n_boot: int = 2000, seed: int = 0):
    n = len(hits)
    pt = sum(hits) / n if n else 0.0
    if n == 0:
        return 0.0, 0.0, 0.0
    rng = _random.Random(seed)
    means = []
    for _ in range(n_boot):
        s = sum(hits[rng.randrange(n)] for _ in range(n)) / n
        means.append(s)
    means.sort()
    lo = means[int(0.025 * n_boot)]
    hi = means[min(int(0.975 * n_boot), n_boot - 1)]
    return pt, lo, hi

def permutation_p(hits, chance, n_perm: int = 10000, seed: int = 0):
    n = len(hits)
    obs = sum(hits)
    if n == 0:
        return 1.0
    rng = _random.Random(seed)
    ge = 0
    for _ in range(n_perm):
        draw = sum(1 for _ in range(n) if rng.random() < chance)
        if draw >= obs:
            ge += 1
    return (ge + 1) / (n_perm + 1)

def tcr_blind_prediction(per_tcr_contacts):
    sums, counts = {}, {}
    for d in per_tcr_contacts:
        for ep, v in d.items():
            if v is None:
                continue
            sums[ep] = sums.get(ep, 0.0) + v
            counts[ep] = counts.get(ep, 0) + 1
    if not sums:
        return None
    means = {ep: sums[ep] / counts[ep] for ep in sums}
    return max(means, key=means.get)

def tcr_blind_accuracy(per_tcr_contacts, cognates):
    pred = tcr_blind_prediction(per_tcr_contacts)
    if pred is None:
        return 0.0
    return sum(1 for cog in cognates if cog == pred) / len(cognates) if cognates else 0.0

def label_permutation_p(observed_top1_mean, per_tcr_contacts, cognates,
                        n_perm=10000, seed=0):
    rng = _random.Random(seed)
    cogs = list(cognates)
    ge = 0
    for _ in range(n_perm):
        perm = cogs[:]
        rng.shuffle(perm)
        hits = 0
        for d, cog in zip(per_tcr_contacts, perm):
            r = retrieval_result({**d, "__scramble__": None}, cog)
            hits += 1 if r["top1"] else 0
        if hits / len(cogs) >= observed_top1_mean:
            ge += 1
    return (ge + 1) / (n_perm + 1)

def paired_contrast(pairs, seed=0, n_boot=2000):
    vals = [(c, s) for (c, s) in pairs if c is not None and s is not None]
    if not vals:
        return {"n": 0, "frac_cognate_higher": None, "mean_delta": None, "ci_delta": [None, None]}
    deltas = [c - s for (c, s) in vals]
    frac = sum(1 for d in deltas if d > 0) / len(deltas)
    rng = _random.Random(seed)
    boots = []
    for _ in range(n_boot):
        boots.append(sum(deltas[rng.randrange(len(deltas))] for _ in range(len(deltas))) / len(deltas))
    boots.sort()
    return {"n": len(deltas), "frac_cognate_higher": frac,
            "mean_delta": sum(deltas) / len(deltas),
            "ci_delta": [boots[int(0.025 * n_boot)], boots[min(int(0.975 * n_boot), n_boot - 1)]]}
