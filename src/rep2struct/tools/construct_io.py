from __future__ import annotations


def parse_fasta(text):
    chains, cur = {}, None
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            cur = line[1:].strip()
            chains[cur] = ""
        elif cur:
            chains[cur] += line
    return chains


def scramble_peptide(pep):
    # deterministic non-identity shuffle: reverse then rotate by 1 (reproducible).
    s = pep[::-1]
    s = s[1:] + s[:1]
    return s if s != pep else pep[1:] + pep[:1]


def pmhc_only(chains: dict) -> dict:
    return {k: v for k, v in chains.items() if k not in ("A", "B")}


def normalize_hla(hla):
    """Collapse a redundant / multi-value HLA to one canonical allele, keeping the
    'HLA-' prefix form the IPD/IMGT lookup uses. The reference DB can return
    'HLA-A*02,HLA-A*02:01'; take the most specific (colon-bearing, longest) token
    so downstream sequence tools resolve a single allele. Idempotent on clean input.
    """
    if not hla:
        return hla
    def _pref(t):
        t = t.strip()
        return t if t.startswith("HLA-") else "HLA-" + t
    toks = [_pref(t) for t in hla.split(",") if t.strip()]
    if not toks:
        return hla
    pool = [t for t in toks if ":" in t] or toks
    return max(pool, key=len)
