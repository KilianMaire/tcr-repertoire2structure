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

def mhc_class_of(hla: str) -> int:
    """Human class II alleles are HLA-DR/DP/DQ (the gene starts with D); everything else
    (HLA-A/B/C, mouse H2-K/D/L) is treated as class I. Keeps the grouping honest instead of
    trusting the hardcoded schema default, and lets build_construct refuse class II, which the
    class I A-E+b2m construct cannot represent. Mouse class II (H2-IA/IE) is not distinguished."""
    if not hla:
        return 1
    gene = hla.upper().replace("HLA-", "").strip()
    return 2 if gene[:1] == "D" else 1


def build_construct(clonotype, annotation, tcr_seqs, mhc_seqs):
    """Build the class I A-E TCR-pMHC construct for one clonotype. Returns a FoldJob, or None
    for a class II allele: the b2m-bearing class I construct would silently mis-model it, so we
    refuse it here and the caller records that it was skipped rather than folding a wrong complex."""
    if mhc_class_of(annotation.hla) == 2:
        return None
    t = tcr_seqs[clonotype.id]
    m = mhc_seqs[annotation.hla]
    fasta = "\n".join([
        ">A", t["A"],            # TCR alpha V domain
        ">B", t["B"],            # TCR beta V domain
        ">C", m["heavy"],        # MHC class I heavy chain
        ">D", m["b2m"],          # beta 2 microglobulin
        ">E", annotation.epitope,  # peptide
    ])
    return FoldJob(clonotype_id=clonotype.id, construct_fasta=fasta,
                   mhc_class=1, tcr_reconstructed=t.get("reconstructed", True))
