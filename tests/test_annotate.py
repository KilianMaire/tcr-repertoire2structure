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

def test_tier_boundary_at_thresholds():
    """Test exact tier boundaries: 12.0 (high), 48.0 (low), 48.01 (unannotatable)."""
    # Distance exactly 12.0 -> tier "high"
    def sim_12(cdr3_a, v_a, cdr3_b, v_b, species="human", top_k=5):
        return ([{"epitope": "HIGH_EP", "mhc": "HLA-A*02:01",
                  "antigen": "Ag1", "distance": 12.0}], "tcrdist", 100, [])
    a = annotate([_clon("c_12", "CDR3_AT_12")], sim_fn=sim_12)[0]
    assert a.annotatable is True
    assert a.confidence_tier == "high"
    assert a.tcrdist == 12.0
    assert a.epitope == "HIGH_EP"

    # Distance exactly 48.0 -> tier "low"
    def sim_48(cdr3_a, v_a, cdr3_b, v_b, species="human", top_k=5):
        return ([{"epitope": "LOW_EP", "mhc": "HLA-B*07:02",
                  "antigen": "Ag2", "distance": 48.0}], "tcrdist", 100, [])
    a = annotate([_clon("c_48", "CDR3_AT_48")], sim_fn=sim_48)[0]
    assert a.annotatable is True
    assert a.confidence_tier == "low"
    assert a.tcrdist == 48.0
    assert a.epitope == "LOW_EP"

    # Distance 48.01 (just above last threshold) -> unannotatable
    def sim_48_01(cdr3_a, v_a, cdr3_b, v_b, species="human", top_k=5):
        return ([{"epitope": "UNREACHABLE", "mhc": "HLA-C*06:02",
                  "antigen": "Ag3", "distance": 48.01}], "tcrdist", 100, [])
    a = annotate([_clon("c_48_01", "CDR3_ABOVE_48")], sim_fn=sim_48_01)[0]
    assert a.annotatable is False
    assert a.confidence_tier == "unannotatable"
    assert a.epitope is None
    assert a.tcrdist == 48.01

def test_min_distance_neighbour_selected():
    """Test that annotate selects the closest neighbour even if not first in list."""
    def sim_multi(cdr3_a, v_a, cdr3_b, v_b, species="human", top_k=5):
        # Return multiple neighbours with different distances; closest is not first
        return ([
            {"epitope": "FIRST_FAR", "mhc": "HLA-A*02:01", "antigen": "Ag_first", "distance": 25.0},
            {"epitope": "CLOSEST", "mhc": "HLA-B*07:02", "antigen": "Ag_close", "distance": 5.0},
            {"epitope": "MIDDLE", "mhc": "HLA-C*06:02", "antigen": "Ag_mid", "distance": 15.0},
        ], "tcrdist", 100, [])
    a = annotate([_clon("c_multi", "CDR3_MULTI_CHOICE")], sim_fn=sim_multi)[0]
    # Should pick the closest (distance=5.0), not the first (distance=25.0)
    assert a.annotatable is True
    assert a.confidence_tier == "high"
    assert a.tcrdist == 5.0
    assert a.epitope == "CLOSEST"
    assert a.antigen == "Ag_close"
