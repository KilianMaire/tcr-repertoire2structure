"""Structure-vs-sequence retrieval benchmark driver.

Pre-fold half: select seed TCRs in one HLA (novel-first, seeded-random within
stratum), build cognate+decoy+scramble constructs, emit constructs + manifest
for the user-driven Colab/H100 notebook. Scoring half is added later.

Usage:
  python scripts/run_benchmark_arm.py emit <dextramer_dir> <out_dir> \
      --hla 'HLA-A*02:01' --n 4 --k 3 --samples 5
"""
from __future__ import annotations
import argparse, json, sys
import random as _random
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from rep2struct import benchmark as bm
from rep2struct.seqs import build_tcr_seqs, build_mhc_seqs
from run_validation_arm import labeled_clonotypes, nearest_cache, annotations_from_cache

def select_seed_tcrs(clonotypes, truth, annotations, hla, n, prefer_novel=True, seed=0,
                     unannotatable_only=False, balance_epitopes=False):
    dist = {a.clonotype_id: getattr(a, "tcrdist", None) for a in annotations}
    annot = {a.clonotype_id: a.annotatable for a in annotations}
    cands = [c for c in clonotypes if truth.get(c.id, (None, None))[1] == hla]
    if unannotatable_only:
        cands = [c for c in cands if not annot.get(c.id, False)]
    novel = [c for c in cands if bm.is_novel(dist.get(c.id))]
    leaked = [c for c in cands if not bm.is_novel(dist.get(c.id))]
    rng = _random.Random(seed)
    rng.shuffle(novel); rng.shuffle(leaked)
    ordered = (novel + leaked) if prefer_novel else (leaked + novel)
    if not balance_epitopes:
        return [c.id for c in ordered[:n]]
    from collections import defaultdict
    by_ep = defaultdict(list)
    for c in ordered:
        by_ep[truth[c.id][0]].append(c)
    eps = sorted(by_ep, key=lambda e: (-len(by_ep[e]), e))
    picked = []
    while len(picked) < n and any(by_ep[e] for e in eps):
        for e in eps:
            if by_ep[e]:
                picked.append(by_ep[e].pop(0))
                if len(picked) >= n:
                    break
    return [c.id for c in picked]

def emit_manifest(out_dir, selected, truth, annotations, panel,
                  tcr_seqs, mhc_seqs, k, samples, clono_by_id):
    out_dir = Path(out_dir)
    (out_dir / "constructs").mkdir(parents=True, exist_ok=True)
    dist = {a.clonotype_id: getattr(a, "tcrdist", None) for a in annotations}
    manifest = {}
    for cid in selected:
        cognate, hla = truth[cid]
        decoys = bm.decoys_for(cognate, hla, panel, k)
        clono = clono_by_id[cid]
        jobs = bm.build_panel_constructs(clono, cognate, hla, decoys,
                                         tcr_seqs, mhc_seqs)
        eps = {}
        for key, job in jobs.items():
            fp = out_dir / "constructs" / f"{cid}__{key}.fasta"
            fp.write_text(job.construct_fasta)
            eps[key] = str(fp)
        manifest[cid] = {"cognate": cognate, "hla": hla, "decoys": decoys,
                         "epitopes": eps, "novel": bm.is_novel(dist.get(cid)),
                         "tcrdist": dist.get(cid), "samples": samples,
                         "cdr3b": clono.cdr3b, "chain_b_seq": tcr_seqs[cid]["B"]}
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest

def _load_truth_and_anns(dextramer_dir, hla):
    clons, labels, hlas = labeled_clonotypes(
        f"{dextramer_dir}/donor1_all_contig_annotations.csv",
        f"{dextramer_dir}/donor1_binarized_matrix.csv")
    truth = {cid: (labels[cid], hlas[cid]) for cid in labels}
    in_hla = [c for c in clons if truth.get(c.id, (None, None))[1] == hla]
    cache = nearest_cache(in_hla)
    anns = annotations_from_cache(in_hla, cache)
    return clons, truth, anns, in_hla

def _emit_cmd(args):
    from rep2struct.ingest import standardize_alleles
    clonotypes, truth, anns, in_hla = _load_truth_and_anns(args.dextramer_dir, args.hla)
    counts = bm.per_hla_novel_counts(in_hla, truth, anns)
    print(json.dumps(counts.get(args.hla, {}), indent=2))
    selected = select_seed_tcrs(in_hla, truth, anns, args.hla, args.n, unannotatable_only=args.unannotatable_only, balance_epitopes=args.balance_epitopes)
    panel = bm.panel_epitopes(truth)
    sel_clonos = [c for c in clonotypes if c.id in set(selected)]
    standardize_alleles(sel_clonos)
    tcr_seqs = build_tcr_seqs(sel_clonos)
    mhc_seqs = build_mhc_seqs([args.hla])
    clono_by_id = {c.id: c for c in sel_clonos}
    emit_manifest(args.out_dir, selected, truth, anns, panel,
                  tcr_seqs, mhc_seqs, args.k, args.samples, clono_by_id)
    print(f"emitted {len(selected)} TCRs x (1+{args.k}+scramble) constructs to {args.out_dir}")

def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("emit")
    e.add_argument("dextramer_dir"); e.add_argument("out_dir")
    e.add_argument("--hla", required=True); e.add_argument("--n", type=int, default=4)
    e.add_argument("--k", type=int, default=3); e.add_argument("--samples", type=int, default=5)
    e.add_argument("--unannotatable-only", action="store_true")
    e.add_argument("--balance-epitopes", action="store_true")
    e.set_defaults(func=_emit_cmd)
    args = ap.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
