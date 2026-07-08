from rep2struct.validate import annotation_metrics
from rep2struct.schema import Annotation

def test_metrics_basic():
    anns = [
        Annotation("c1", True, "high", 3.0, epitope="GILGFVFTL"),   # correct
        Annotation("c2", True, "medium", 20.0, epitope="NLVPMVATV"),# wrong
        Annotation("c3", False, "unannotatable"),                   # missed
    ]
    labels = {"c1": "GILGFVFTL", "c2": "GILGFVFTL", "c3": "KLGGALQAK"}
    m = annotation_metrics(anns, labels)
    assert m["precision"] == 0.5           # 1 of 2 annotated correct
    assert round(m["recall"], 3) == 0.333  # 1 of 3 labeled correct
    assert round(m["unannotatable_rate"], 3) == 0.333
