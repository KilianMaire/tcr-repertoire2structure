from pathlib import Path
from rep2struct.ingest import parse_10x

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
    # AAAF has only a non-productive beta; AAAG has only an alpha -> both dropped
    ids = {(c.trav, c.cdr3b) for c in clons}
    assert ("TRAV8-1", "CASSIRSSYEQYF") not in ids
    assert report["dropped_unpaired"] >= 2
