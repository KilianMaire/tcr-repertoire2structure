from __future__ import annotations
from collections import defaultdict
from pathlib import Path as _Path
import glob as _glob
import json as _json
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

def contact_by_epitope(paths_by_epitope, cdr3b=None):
    return {ep: ensemble_contact(paths, cdr3b)[0] for ep, paths in paths_by_epitope.items()}

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

# --- Structural-confidence readouts ---------------------------------------
# Protenix writes one summary_confidence_sample_N.json per sample next to the
# CIFs. Chain order in the construct is fixed: A=TCRalpha, B=TCRbeta, C=MHC
# heavy, D=beta2m, E=peptide. chain_pair_iptm[i][j] is the interface predicted
# TM between chains i,j (higher = better); chain_pair_gpde[i][j] is a PAE-analog
# (lower = better). The readouts below rank the panel by a confidence signal
# rather than by the (refuted) CDR3beta-peptide contact count.
_A, _B, _C, _D, _E = 0, 1, 2, 3, 4

_CONF_READOUTS = {
    # confidence at the TCR-peptide interface (the mechanistically-relevant one)
    "iptm_TCRpep_max":   lambda c: max(c["chain_pair_iptm"][_A][_E], c["chain_pair_iptm"][_B][_E]),
    "iptm_TCRpep_mean":  lambda c: (c["chain_pair_iptm"][_A][_E] + c["chain_pair_iptm"][_B][_E]) / 2,
    "iptm_beta_pep":     lambda c: c["chain_pair_iptm"][_B][_E],
    "iptm_alpha_pep":    lambda c: c["chain_pair_iptm"][_A][_E],
    "neg_gpde_beta_pep": lambda c: -c["chain_pair_gpde"][_B][_E],
    # whole-complex confidence
    "iptm_global":       lambda c: c["iptm"],
    "ptm_global":        lambda c: c["ptm"],
    "ranking_score":     lambda c: c["ranking_score"],
    # NEGATIVE CONTROL: the MHC-peptide groove must NOT carry TCR specificity;
    # if this ranks the cognate it would mean the signal is peptide-in-groove,
    # not TCR-recognition. It ranks BELOW chance on the panel, as it should.
    "iptm_groove_ctrl":  lambda c: c["chain_pair_iptm"][_C][_E],
}

def _confidence_samples(fold_dir):
    out = []
    for jp in _glob.glob(str(_Path(fold_dir) / "**" / "*summary_confidence_sample_*.json"),
                         recursive=True):
        try:
            out.append(_json.loads(_Path(jp).read_text()))
        except (ValueError, OSError):
            continue
    return out

def confidence_readout_by_epitope(manifest_entry, folds_root, cid, readout):
    """Median value of one confidence readout per epitope (scramble excluded)."""
    fn = _CONF_READOUTS[readout]
    out = {}
    for ep in manifest_entry["epitopes"]:
        if ep == "__scramble__":
            continue
        vals = []
        for c in _confidence_samples(_Path(folds_root) / f"{cid}__{ep}"):
            try:
                vals.append(float(fn(c)))
            except (KeyError, IndexError, TypeError):
                continue
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

def _paths_by_epitope(manifest_entry, folds_root, cid):
    out = {}
    for ep in manifest_entry["epitopes"]:
        d = _Path(folds_root) / f"{cid}__{ep}"
        out[ep] = sorted(str(p) for p in d.rglob("*.cif")) if d.exists() else []
    return out

def evaluate(manifest, folds_root, annotations):
    seq_ep = {a.clonotype_id: getattr(a, "epitope", None) for a in annotations}
    per = {}
    for cid, ent in manifest.items():
        pbe = _paths_by_epitope(ent, folds_root, cid)
        contacts = contact_by_epitope(pbe, ent.get("cdr3b"))
        cog = ent["cognate"]
        if contacts.get(cog) is None:
            continue
        panel_contacts = {e: v for e, v in contacts.items() if e != "__scramble__"}
        plddts = cdr3b_plddt_by_epitope(
            {e: p for e, p in pbe.items() if e != "__scramble__"},
            ent.get("chain_b_seq"), ent.get("cdr3b"))
        decoy_c = [contacts[e] for e in ent["decoys"] if contacts.get(e) is not None]
        conf = {r: confidence_readout_by_epitope(ent, folds_root, cid, r)
                for r in _CONF_READOUTS}
        per[cid] = {
            "novel": ent["novel"],
            "panel_contacts": panel_contacts,
            "confidence": conf,
            "contact_top1": 1.0 if retrieval_result(contacts, cog)["top1"] else 0.0,
            "plddt_top1": 1.0 if retrieval_result(plddts, cog)["top1"] else 0.0,
            "seq_top1": 1.0 if sequence_baseline_top1(seq_ep.get(cid), cog) else 0.0,
            "cognate_contact": contacts.get(cog),
            "scramble_contact": contacts.get("__scramble__"),
            "decoy_contacts": decoy_c,
            "cognate": cog,
            "n_panel": 1 + len(ent["decoys"]),
        }

    def strata(rows):
        if not rows:
            return {"n": 0}
        chance = sum(1.0 / r["n_panel"] for r in rows) / len(rows)
        contact_hits = [r["contact_top1"] for r in rows]
        pt, lo, hi = bootstrap_ci(contact_hits)
        blind_acc = tcr_blind_accuracy([r["panel_contacts"] for r in rows],
                                       [r["cognate"] for r in rows])
        perm_p = label_permutation_p(pt, [r["panel_contacts"] for r in rows],
                                     [r["cognate"] for r in rows])
        pc = paired_contrast([(r["cognate_contact"], r["scramble_contact"]) for r in rows])
        cogs = [r["cognate"] for r in rows]

        def _confidence_stratum(readout):
            panels = [r["confidence"][readout] for r in rows]
            hits = [1.0 if retrieval_result(p, cog)["top1"] else 0.0
                    for p, cog in zip(panels, cogs)]
            cpt, clo, chi = bootstrap_ci(hits)
            return {"top1": cpt, "ci": [clo, chi],
                    "blind": tcr_blind_accuracy(panels, cogs),
                    "perm_p": label_permutation_p(cpt, panels, cogs)}

        return {
            "n": len(rows), "chance": chance, "tcr_blind_acc": blind_acc,
            "contact": {"top1": pt, "ci": [lo, hi],
                        "p_vs_chance": permutation_p(contact_hits, chance),
                        "p_vs_blind": perm_p,
                        "auroc": auroc([(r["cognate_contact"], r["decoy_contacts"])
                                        for r in rows if r["cognate_contact"] is not None])},
            "plddt": {"top1": sum(r["plddt_top1"] for r in rows) / len(rows)},
            "seq": {"top1": sum(r["seq_top1"] for r in rows) / len(rows)},
            "confidence": {r: _confidence_stratum(r) for r in _CONF_READOUTS},
            "scramble_contrast": pc,
        }

    rows = list(per.values())
    return {"per_tcr": per,
            "overall": strata(rows),
            "novel": strata([r for r in rows if r["novel"]]),
            "leaked": strata([r for r in rows if not r["novel"]])}
