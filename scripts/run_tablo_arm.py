"""Application arm: run the pipeline on a real human repertoire (TABLO).

TABLO ships a derived paired-clonotype table, not 10x contig rows, so we build
Clonotype objects directly and drive the real stage functions: allele
standardization, TCRdist annotation against VDJdb/IEDB (never a forced label),
real reconstructed class I constructs. The fold is left pending here (no GPU in
this environment); QC runs once the fold agent produces structures on Colab.

Usage:
  python scripts/run_tablo_arm.py <clonotypes.csv> <run_dir> [top_n] [donor]
"""
from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from rep2struct.ingest import _clon_id, standardize_alleles
from rep2struct.annotate import annotate
from rep2struct.foldprep import select_top, build_construct
from rep2struct.seqs import build_tcr_seqs, build_mhc_seqs
from rep2struct.report import render_report
from rep2struct.runstate import RunState
from rep2struct.schema import Clonotype, QCResult


def load_tablo_cd8(csv_path, min_public=3):
    """Screen the public CD8 compartment. Public (donor-shared) clones are
    enriched for known antigen specificities, so this is where TCRdist
    annotation is informative; the vast private repertoire stays unannotatable
    by design. Dedup by the clonotype tuple (a public clone appears once per
    donor), summing cells and keeping the max donor breadth."""
    df = pd.read_csv(csv_path)
    df = df[(df["lineage"].astype(str) == "CD8") &
            (df["public_n_donors"] >= min_public)]
    df = df.dropna(subset=["v_a_gene", "cdr3_a_aa", "v_b_gene", "cdr3_b_aa"])
    agg = (df.groupby(["v_a_gene", "cdr3_a_aa", "v_b_gene", "cdr3_b_aa",
                       "j_a_gene", "j_b_gene"], dropna=False)
             .agg(n_cells=("n_cells", "sum"),
                  public_n_donors=("public_n_donors", "max"))
             .reset_index())
    agg = agg.sort_values(["public_n_donors", "n_cells"], ascending=False)
    clons = [
        Clonotype(
            id=_clon_id(r.v_a_gene, r.cdr3_a_aa, r.v_b_gene, r.cdr3_b_aa),
            trav=r.v_a_gene, cdr3a=r.cdr3_a_aa, trbv=r.v_b_gene, cdr3b=r.cdr3_b_aa,
            size=int(r.n_cells), traj=r.j_a_gene, trbj=r.j_b_gene)
        for r in agg.itertuples()
    ]
    return clons


def main():
    csv_path, run_dir = sys.argv[1], sys.argv[2]
    top_n = int(sys.argv[3]) if len(sys.argv) > 3 else 12
    min_public = int(sys.argv[4]) if len(sys.argv) > 4 else 3

    clons = load_tablo_cd8(csv_path, min_public=min_public)
    print(f"screened {len(clons)} unique public CD8 clonotypes "
          f"(shared across >= {min_public} donors)")

    clons = standardize_alleles(clons)
    anns = annotate(clons)
    n_ann = sum(a.annotatable for a in anns)
    print(f"annotation: {n_ann}/{len(anns)} annotatable; tiers "
          f"{ {t: sum(a.confidence_tier==t for a in anns) for t in ('high','medium','low','unannotatable')} }")

    top = select_top(clons, anns, n=top_n)
    foldable = [(c, a) for c, a in top if a.annotatable and a.hla]
    seqs = build_tcr_seqs([c for c, _ in foldable])
    mhc = build_mhc_seqs(sorted({a.hla for _, a in foldable}))
    jobs = [build_construct(c, a, seqs, mhc) for c, a in foldable if a.hla in mhc]
    n_real = sum(seqs[c.id].get("reconstructed", False) for c, _ in foldable)
    print(f"fold prep: {len(jobs)} constructs built; "
          f"{n_real}/{len(foldable)} with fully reconstructed V domains; "
          f"HLAs resolved: {sorted(mhc)}")

    rs = RunState(run_dir)
    rs.write_stage("ingest", clons)
    rs.write_stage("annotate", anns)
    rs.write_stage("foldjobs", jobs)
    for j in jobs:  # persist the real construct FASTA for the fold agent
        (Path(run_dir) / f"construct_{j.clonotype_id}.fasta").write_text(j.construct_fasta)

    # No GPU here: folds are pending, so QC is honestly "not folded".
    qcs = [QCResult(j.clonotype_id, "qc_failed", "fold pending (run the fold agent on Colab)")
           for j in jobs]
    html = render_report(clons, anns, qcs, metrics=None)
    out = Path(run_dir) / "report.html"
    out.write_text(html)
    print(f"report -> {out}")


if __name__ == "__main__":
    main()
