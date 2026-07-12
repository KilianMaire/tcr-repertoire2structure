import json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "science" / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import run_benchmark_arm as drv
from rep2struct.schema import Clonotype, Annotation
from rep2struct import benchmark as bm2

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

def _make_folds(root, cid, cognate, decoy, fix):
    for ep, src in [(cognate, fix/"cognate_min.cif"), (decoy, fix/"scramble_min.cif"),
                    ("__scramble__", fix/"scramble_min.cif")]:
        d = root / "folds" / f"{cid}__{ep}"
        d.mkdir(parents=True)
        (d / "sample_0.cif").write_bytes(src.read_bytes())

def test_evaluate_stratifies_and_scores(tmp_path):
    fix = Path(__file__).parent / "fixtures"
    cid = "a"
    _make_folds(tmp_path, cid, "GILGFVFTL", "NLVPMVATV", fix)
    manifest = {cid: {"cognate": "GILGFVFTL", "hla": "HLA-A*02:01",
                      "decoys": ["NLVPMVATV"],
                      "epitopes": {"GILGFVFTL": "", "NLVPMVATV": "", "__scramble__": ""},
                      "novel": True, "tcrdist": None, "samples": 1,
                      "cdr3b": "AAAA", "chain_b_seq": "AAAA"}}
    anns = [Annotation(cid, False, "unannotatable", tcrdist=None)]
    out = bm2.evaluate(manifest, tmp_path / "folds", anns)
    assert out["novel"]["n"] == 1
    assert "scramble_contrast" in out["novel"]
    assert out["novel"]["tcr_blind_acc"] in (0.0, 1.0)
    assert out["novel"]["contact"]["top1"] in (0.0, 1.0)
    assert out["novel"]["seq"]["top1"] == 0.0
    # confidence readouts are wired even when no summary_confidence JSONs exist:
    # every readout resolves to a stratum (None panels -> Top-1 0.0), never a crash
    conf = out["novel"]["confidence"]
    assert "iptm_TCRpep_max" in conf and "iptm_groove_ctrl" in conf
    assert conf["iptm_TCRpep_max"]["top1"] == 0.0


def test_confidence_readout_reads_summary_json(tmp_path):
    # A construct dir with two samples; iptm_TCRpep_max = max(iptm[A][E], iptm[B][E]).
    # Median over the two samples must be the wired value; missing keys are skipped.
    d = tmp_path / "folds" / "a__GILGFVFTL" / "seed0"
    d.mkdir(parents=True)
    def _cpi(ae, be, ce):
        m = [[0.0] * 5 for _ in range(5)]
        m[0][4] = ae; m[1][4] = be; m[2][4] = ce
        return m
    for i, (ae, be, ce) in enumerate([(0.3, 0.5, 0.1), (0.4, 0.6, 0.2)]):
        (d / f"x_summary_confidence_sample_{i}.json").write_text(json.dumps(
            {"chain_pair_iptm": _cpi(ae, be, ce),
             "chain_pair_gpde": _cpi(1.0, 2.0, 3.0),
             "iptm": 0.4, "ptm": 0.5, "ranking_score": 0.6}))
    ent = {"epitopes": {"GILGFVFTL": "", "__scramble__": ""}}
    vals = bm2.confidence_readout_by_epitope(ent, tmp_path / "folds", "a", "iptm_TCRpep_max")
    assert "__scramble__" not in vals
    # samples give max(0.3,0.5)=0.5 and max(0.4,0.6)=0.6 -> median 0.55
    assert abs(vals["GILGFVFTL"] - 0.55) < 1e-9
    groove = bm2.confidence_readout_by_epitope(ent, tmp_path / "folds", "a", "iptm_groove_ctrl")
    assert abs(groove["GILGFVFTL"] - 0.15) < 1e-9  # median of 0.1, 0.2
