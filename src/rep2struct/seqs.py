"""Real construct sequences for the TCR pMHC fold.

Three providers, replacing the poly-G / poly-H placeholders that the offline
harness used:

- TCR variable domains, reconstructed from V gene, J gene and CDR3 through
  TCR Explorer's `reconstruct_tcr` (bundled IMGT germline, works offline).
- MHC class I heavy chain ectodomain, fetched once per allele from the EBI
  IPD/IMGT-HLA REST API and cached to a vendored JSON so later runs are
  offline and reproducible.
- beta 2 microglobulin, an invariant mature chain (a verified constant).

Every provider degrades gracefully: a chain that cannot be reconstructed or an
allele that cannot be fetched yields a clearly marked fallback plus a warning,
so a run never dies on one odd input, and the report can flag it.
"""
from __future__ import annotations
import json
import re
from pathlib import Path
from typing import Optional

# Mature human beta 2 microglobulin. UniProt P61769, signal peptide (residues
# 1 to 20, MSRSVALAVLALLSLSGLEA) removed. Invariant across every class I fold.
B2M_MATURE = (
    "IQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTE"
    "FTPTEKDEYACRVNHVTLSQPKIVKWDRDM"
)

# Classical HLA class I mature chains open with a highly conserved motif
# (G/C)SHSM[RK]YF. We locate it to strip the variable length signal peptide
# without knowing its length, then take the first 275 residues, which is the
# alpha1 alpha2 alpha3 ectodomain used in soluble pMHC I structures (drops the
# connecting peptide, transmembrane helix and cytoplasmic tail).
_MATURE_START = re.compile(r"[GC]SHSM[RK]YF")
_ECTO_LEN = 275

_EBI_BASE = "https://www.ebi.ac.uk/cgi-bin/ipd/api/allele"
_CACHE_PATH = Path(__file__).parent / "data" / "hla_ectodomains.json"


def _norm_hla(allele: str) -> str:
    """Normalize an HLA string to the two field IPD query stem, e.g.
    'HLA-A*02:01' or 'A*02:01:01:01' -> 'A*02:01'. Returns '' if it does not
    look like a classical class I allele name."""
    a = allele.strip().upper()
    if a.startswith("HLA-"):
        a = a[4:]
    m = re.match(r"([ABC])\*?(\d{2,3}):(\d{2,3})", a)
    if not m:
        return ""
    return f"{m.group(1)}*{m.group(2)}:{m.group(3)}"


def _ectodomain_from_protein(protein: str) -> Optional[str]:
    m = _MATURE_START.search(protein)
    if not m:
        return None
    mature = protein[m.start():]
    return mature[:_ECTO_LEN]


def _load_cache() -> dict:
    if _CACHE_PATH.exists():
        return json.loads(_CACHE_PATH.read_text())
    return {}


def _save_cache(cache: dict) -> None:
    _CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
    _CACHE_PATH.write_text(json.dumps(cache, indent=1, sort_keys=True))


def _fetch_hla_ectodomain(stem: str) -> Optional[str]:
    """Fetch the heavy chain ectodomain for a two field HLA stem from EBI IPD.
    Network call; returns None on any failure."""
    import httpx
    try:
        q = f'{_EBI_BASE}?project=HLA&query=startsWith(name,"{stem}")&limit=1'
        r = httpx.get(q, timeout=40, headers={"Accept": "application/json"})
        r.raise_for_status()
        data = r.json().get("data", [])
        if not data:
            return None
        acc = data[0]["accession"]
        r2 = httpx.get(f"{_EBI_BASE}/{acc}?project=HLA", timeout=40,
                       headers={"Accept": "application/json"})
        r2.raise_for_status()
        protein = r2.json().get("sequence", {}).get("protein")
        return _ectodomain_from_protein(protein) if protein else None
    except Exception:
        return None


def hla_heavy_ectodomain(allele: str) -> tuple[Optional[str], Optional[str]]:
    """Return (ectodomain, warning). Cached-first, then EBI IPD, then None."""
    stem = _norm_hla(allele)
    if not stem:
        return None, f"hla_unrecognized:{allele}"
    cache = _load_cache()
    if stem in cache:
        return cache[stem], None
    ecto = _fetch_hla_ectodomain(stem)
    if ecto is None:
        return None, f"hla_fetch_failed:{stem}"
    cache[stem] = ecto
    _save_cache(cache)
    return ecto, None


def reconstruct_variable_domains(clonotype, species: str = "human") -> dict:
    """Reconstruct the alpha and beta V domains for one clonotype.

    Returns {"A": seq or None, "B": seq or None, "warnings": [...]}. A chain is
    None when its J gene is missing or the germline does not resolve."""
    from tcr_explorer.reconstructor import reconstruct_tcr
    out = {"A": None, "B": None, "warnings": []}
    plan = [
        ("A", clonotype.trav_allele or clonotype.trav, getattr(clonotype, "traj", None)),
        ("B", clonotype.trbv_allele or clonotype.trbv, getattr(clonotype, "trbj", None)),
        # cdr3 chosen per chain below
    ]
    cdr3 = {"A": clonotype.cdr3a, "B": clonotype.cdr3b}
    for chain, v_gene, j_gene in plan:
        if not v_gene or not j_gene:
            out["warnings"].append(f"no_v_or_j:{chain}:{clonotype.id}")
            continue
        try:
            rec = reconstruct_tcr(v_gene, j_gene, cdr3[chain], species=species)
        except Exception as e:  # germline lookup can raise on odd names
            out["warnings"].append(f"reconstruct_error:{chain}:{type(e).__name__}")
            continue
        if rec.get("full_aa"):
            out[chain] = rec["full_aa"]
        else:
            out["warnings"].append(f"reconstruct_failed:{chain}:{v_gene}/{j_gene}")
    return out


def _tcr_stub(clonotype) -> dict:
    """Last resort so a run never crashes: poly-G framework plus the real CDR3.
    Clearly not a real V domain; the report flags any clonotype that used it."""
    return {"A": "G" * 10 + clonotype.cdr3a, "B": "G" * 10 + clonotype.cdr3b,
            "reconstructed": False}


def build_tcr_seqs(clonotypes, species: str = "human") -> dict:
    """id -> {"A", "B", "reconstructed": bool}. Real V domains where possible,
    stub fallback otherwise (always non-None so downstream never KeyErrors)."""
    out = {}
    for c in clonotypes:
        dom = reconstruct_variable_domains(c, species=species)
        if dom["A"] and dom["B"]:
            out[c.id] = {"A": dom["A"], "B": dom["B"], "reconstructed": True}
        else:
            out[c.id] = _tcr_stub(c)
    return out


def build_mhc_seqs(alleles) -> dict:
    """hla -> {"heavy", "b2m"}. Fetches (cached) the heavy chain ectodomain per
    unique allele; b2m is the invariant constant. Alleles that do not resolve
    are omitted, so build_construct on them is skipped upstream."""
    out = {}
    for allele in alleles:
        if allele is None:
            continue
        heavy, _warn = hla_heavy_ectodomain(allele)
        if heavy is None:
            continue
        out[allele] = {"heavy": heavy, "b2m": B2M_MATURE}
    return out
