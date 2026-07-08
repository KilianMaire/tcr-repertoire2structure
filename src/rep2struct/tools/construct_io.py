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
