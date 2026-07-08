from pathlib import Path
from rep2struct.ingest import parse_10x, standardize_alleles
from rep2struct.schema import Clonotype

FIX = Path(__file__).parent / "fixtures" / "tenx_tiny.csv"

def test_collapse_identical_clonotypes():
    clons = parse_10x(FIX)
    # AAAC and AAAD share the exact alpha+beta tuple -> one clonotype, size 2
    top = [c for c in clons if c.cdr3b == "CASSIRSSYEQYF"]
    assert len(top) == 1
    assert top[0].size == 2
    assert top[0].trav == "TRAV1-2"

def test_unpaired_rows_dropped_with_reason():
    clons, report = parse_10x(FIX, report=True)
    # AAAF has only a beta (unpaired); AAAG has only an alpha (unpaired) -> both dropped as unpaired
    ids = {(c.trav, c.cdr3b) for c in clons}
    assert ("TRAV8-1", "CASSIRSSYEQYF") not in ids
    assert report["dropped_unpaired"] >= 2

def test_productive_filter_drops_paired_cell():
    clons, report = parse_10x(FIX, report=True)
    # AAAH-1 has both alpha and beta but productive=False -> filtered before unpaired check
    alphas = {c.cdr3a for c in clons}
    assert "CAVINNDYKLSF" not in alphas
    # dropped_unpaired should still be 2 (only AAAF and AAAG), not counting AAAH
    assert report["dropped_unpaired"] == 2

def test_standardize_alleles_fills_allele_or_keeps_none():
    clons = [Clonotype(id="x", trav="TRAV1-2", cdr3a="CAVMDSSYKLIF",
                       trbv="TRBV19", cdr3b="CASSIRSSYEQYF", size=2)]
    def fake_assign(gene, species, chain):
        return {"TRAV1-2": "TRAV1-2*01", "TRBV19": None}.get(gene)
    out = standardize_alleles(clons, assign_fn=fake_assign)
    assert out[0].trav_allele == "TRAV1-2*01"
    assert out[0].trbv_allele is None  # failure keeps None, clonotype kept
