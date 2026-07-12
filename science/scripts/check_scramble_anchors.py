"""Does anchor loss in the scramble explain the A*02:01 vs A*11:01 presentation gap?

Read-only. Tests the intuition that A*11:01 scrambles separate cleanly because a
shuffle destroys the strict C-terminal K/R anchor, while A*02:01 scrambles stay
bindable. Two checks:

  1. Anchor motif retention: fraction of scrambles that still satisfy the allele
     anchor motif (A*02:01: hydrophobic P2 and C-terminus; A*11:01: basic
     C-terminus). If similar across alleles, anchor loss is not the driver.
  2. Absolute groove confidence: median iptm_groove for binders vs scrambles per
     allele. This shows where the separation actually comes from.

Conclusion from the shipped runs: retention is similar (60% vs 50%), so the gap
is not anchor loss; it is that the permissive A*02:01 groove keeps model
confidence high even for a scramble, compressing the binder-vs-scramble gap.

Usage: python scripts/check_scramble_anchors.py runs/panel1 runs/hla_a1101
"""
from __future__ import annotations
import glob, json, statistics as st, sys
from pathlib import Path

A0201_P2 = set("LMIVAT"); A0201_PO = set("VLIAMT")
A1101_PO = set("KR")
MOTIF = {
    "HLA-A*02:01": lambda p: len(p) >= 3 and p[1] in A0201_P2 and p[-1] in A0201_PO,
    "HLA-A*11:01": lambda p: len(p) >= 3 and p[-1] in A1101_PO,
}


def _pep_chain_E(fasta_path):
    cur, seq = None, {}
    for line in Path(fasta_path).read_text().splitlines():
        if line.startswith(">"):
            cur = line[1:].strip(); seq[cur] = ""
        elif cur:
            seq[cur] += line.strip()
    return seq.get("E")


def _median_groove(cid, ep, folds_root):
    vals = []
    for jp in glob.glob(str(Path(folds_root) / f"{cid}__{ep}" / "**"
                             / "*summary_confidence_sample_*.json"), recursive=True):
        try:
            vals.append(json.loads(Path(jp).read_text())["chain_pair_iptm"][2][4])
        except (ValueError, OSError, KeyError, IndexError, TypeError):
            continue
    return st.median(vals) if vals else None


def analyze(run_dir):
    run = Path(run_dir)
    manifest = json.loads((run / "manifest.json").read_text())
    hla = next(iter(manifest.values()))["hla"]
    motif = MOTIF.get(hla, lambda p: None)
    cog_ok = scr_ok = n = 0
    binders, scrambles = [], []
    for cid, ent in manifest.items():
        scr_pep = _pep_chain_E(run / "constructs" / f"{cid}____scramble__.fasta")
        if scr_pep is None:
            continue
        n += 1
        cog_ok += bool(motif(ent["cognate"]))
        scr_ok += bool(motif(scr_pep))
        for ep in [ent["cognate"], *ent["decoys"]]:
            v = _median_groove(cid, ep, run / "folds")
            if v is not None:
                binders.append(v)
        sv = _median_groove(cid, "__scramble__", run / "folds")
        if sv is not None:
            scrambles.append(sv)
    print(f"=== {run.name}  ({hla}, n={n}) ===")
    print(f"  anchor motif retained:  cognate {cog_ok}/{n} ({cog_ok/n:.0%})   "
          f"scramble {scr_ok}/{n} ({scr_ok/n:.0%})")
    bm, sm = st.median(binders), st.median(scrambles)
    print(f"  iptm_groove median:     binder {bm:.3f}   scramble {sm:.3f}   drop {bm-sm:+.3f}")
    print(f"    binder IQR   [{st.quantiles(binders)[0]:.3f}, {st.quantiles(binders)[2]:.3f}]")
    print(f"    scramble IQR [{st.quantiles(scrambles)[0]:.3f}, {st.quantiles(scrambles)[2]:.3f}]")
    print()


if __name__ == "__main__":
    for d in (sys.argv[1:] or ["runs/panel1", "runs/hla_a1101"]):
        analyze(d)
