from rep2struct.foldprep import select_top, build_construct, TIER_WEIGHT
from rep2struct.schema import Clonotype, Annotation

def _c(cid, size): return Clonotype(cid, "TRAV1-2", "CAVA", "TRBV19", "CASSB", size)

def test_ranking_prefers_confident_and_expanded():
    clons = [_c("c1", 2), _c("c2", 100)]
    anns = [Annotation("c1", True, "high", 3.0, epitope="GILGFVFTL"),
            Annotation("c2", False, "unannotatable")]
    top = select_top(clons, anns, n=1)
    assert top[0][0].id == "c1"   # high tier beats big-but-unannotatable

def test_build_construct_has_five_chains():
    c = _c("c1", 2)
    a = Annotation("c1", True, "high", 3.0, epitope="GILGFVFTL", hla="HLA-A*02:01")
    job = build_construct(c, a,
        tcr_seqs={"c1": {"A": "AAAA", "B": "BBBB"}},
        mhc_seqs={"HLA-A*02:01": {"heavy": "HHHH", "b2m": "MMMM"}})
    chains = [l for l in job.construct_fasta.splitlines() if l.startswith(">")]
    assert chains == [">A", ">B", ">C", ">D", ">E"]
    assert "GILGFVFTL" in job.construct_fasta
