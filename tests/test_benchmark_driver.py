import json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import run_benchmark_arm as drv
from rep2struct.schema import Clonotype, Annotation

def test_select_seed_novel_first():
    clonos = [Clonotype("a","TRAV1","CAA","TRBV1","CBB",2),
              Clonotype("b","TRAV1","CAC","TRBV1","CBD",9)]
    truth = {"a": ("GILGFVFTL","HLA-A*02:01"), "b": ("GILGFVFTL","HLA-A*02:01")}
    anns = [Annotation("a", False, "unannotatable", tcrdist=None),
            Annotation("b", True, "high", tcrdist=0.0)]
    sel = drv.select_seed_tcrs(clonos, truth, anns, "HLA-A*02:01", n=2)
    assert sel[0] == "a"

def test_select_unannotatable_only_excludes_annotatable():
    clonos = [Clonotype("a","TRAV1","CAA","TRBV1","CBB",2),
              Clonotype("b","TRAV1","CAC","TRBV1","CBD",9)]
    truth = {"a": ("GILGFVFTL","HLA-A*02:01"), "b": ("GILGFVFTL","HLA-A*02:01")}
    anns = [Annotation("a", False, "unannotatable", tcrdist=None),
            Annotation("b", True, "low", tcrdist=18.0)]
    sel = drv.select_seed_tcrs(clonos, truth, anns, "HLA-A*02:01", n=2, unannotatable_only=True)
    assert sel == ["a"]

def test_emit_manifest_writes_constructs(tmp_path):
    clonos = [Clonotype("a","TRAV1","CAA","TRBV1","CBB",2)]
    truth = {"a": ("GILGFVFTL","HLA-A*02:01")}
    anns = [Annotation("a", False, "unannotatable", tcrdist=None)]
    panel = [("GILGFVFTL","HLA-A*02:01"), ("NLVPMVATV","HLA-A*02:01")]
    tcr_seqs = {"a": {"A": "AAAA", "B": "BBBB"}}
    mhc_seqs = {"HLA-A*02:01": {"heavy": "HHHH", "b2m": "MMMM"}}
    clono_by_id = {"a": clonos[0]}
    man = drv.emit_manifest(tmp_path, ["a"], truth, anns, panel,
                            tcr_seqs, mhc_seqs, k=1, samples=5, clono_by_id=clono_by_id)
    assert man["a"]["novel"] is True
    assert set(man["a"]["epitopes"]) == {"GILGFVFTL", "NLVPMVATV", "__scramble__"}
    assert man["a"]["cdr3b"] == "CBB"
    assert man["a"]["chain_b_seq"] == "BBBB"
    assert (tmp_path / "manifest.json").exists()
    for p in man["a"]["epitopes"].values():
        assert Path(p).exists()

def test_select_balance_epitopes_spreads_across_cognates():
    # 3 TCR for epitope P, 1 for Q; balanced pick of 2 must include Q, not two P
    clonos = [Clonotype(x,"TRAV1","C"+x,"TRBV1","CB"+x,1) for x in ["a","b","c","d"]]
    truth = {"a":("P","HLA-A*02:01"),"b":("P","HLA-A*02:01"),
             "c":("P","HLA-A*02:01"),"d":("Q","HLA-A*02:01")}
    anns = [Annotation(x, False, "unannotatable", tcrdist=None) for x in ["a","b","c","d"]]
    sel = drv.select_seed_tcrs(clonos, truth, anns, "HLA-A*02:01", n=2,
                               balance_epitopes=True)
    cogs = sorted(truth[cid][0] for cid in sel)
    assert cogs == ["P", "Q"]   # one from each epitope, not P+P
