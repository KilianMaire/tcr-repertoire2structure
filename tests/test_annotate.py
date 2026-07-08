from rep2struct.annotate import annotate, DEFAULT_TIERS
from rep2struct.schema import Clonotype

def _clon(cid, cdr3b): return Clonotype(id=cid, trav="TRAV1-2", cdr3a="CAVA",
                                        trbv="TRBV19", cdr3b=cdr3b, size=1)

def test_close_neighbour_is_annotated_high():
    def sim(cdr3_a, v_a, cdr3_b, v_b, species="human", top_k=5):
        return ([{"epitope": "GILGFVFTL", "mhc": "HLA-A*02:01",
                  "antigen": "Flu M1", "distance": 3.0}], "tcrdist", 100, [])
    a = annotate([_clon("c1", "CASSIRSSYEQYF")], sim_fn=sim)[0]
    assert a.annotatable and a.confidence_tier == "high"
    assert a.epitope == "GILGFVFTL" and a.tcrdist == 3.0

def test_far_neighbour_is_unannotatable():
    def sim(*args, **kw):
        return ([{"epitope": "X", "mhc": "Y", "antigen": "Z", "distance": 999.0}], "tcrdist", 100, [])
    a = annotate([_clon("c2", "CASSNOMATCH")], sim_fn=sim)[0]
    assert a.annotatable is False
    assert a.confidence_tier == "unannotatable"
    assert a.epitope is None

def test_no_neighbour_is_unannotatable():
    def sim(*args, **kw): return ([], "tcrdist", 0, ["no candidates"])
    a = annotate([_clon("c3", "CASSEMPTY")], sim_fn=sim)[0]
    assert a.annotatable is False and a.epitope is None
