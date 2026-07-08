from pathlib import Path
from rep2struct.pipeline import run_pipeline

FIX = Path(__file__).parent / "fixtures"


def _sim(cdr3_a, v_a, cdr3_b, v_b, species="human", top_k=5):
    if cdr3_b == "CASSIRSSYEQYF":
        return ([{"epitope": "GILGFVFTL", "mhc": "HLA-A*02:01",
                  "antigen": "Flu M1", "distance": 3.0}], "tcrdist", 1, [])
    return ([], "tcrdist", 0, [])


def _assign_fn(gene, species="human", chain=None):
    return gene + "*01"


def _fold_fn(job):
    return [str(FIX / "cognate_min.cif")]


def test_offline_end_to_end(tmp_path):
    report = run_pipeline(
        csv_path=FIX / "tenx_tiny.csv", run_dir=tmp_path / "run",
        top_n=1, sim_fn=_sim, assign_fn=_assign_fn, fold_fn=_fold_fn,
        tcr_seqs=None, mhc_seqs={"HLA-A*02:01": {"heavy": "H" * 20, "b2m": "M" * 20}},
        scramble_threshold=1.0)
    html = Path(report).read_text()
    assert "GILGFVFTL" in html and "reliable" in html


def test_resume_is_idempotent(tmp_path):
    # running twice must not raise and must return the same report path
    run_dir = tmp_path / "run"
    kwargs = dict(
        csv_path=FIX / "tenx_tiny.csv", run_dir=run_dir,
        top_n=1, sim_fn=_sim, assign_fn=_assign_fn, fold_fn=_fold_fn,
        tcr_seqs=None, mhc_seqs={"HLA-A*02:01": {"heavy": "H" * 20, "b2m": "M" * 20}},
        scramble_threshold=1.0)
    report1 = run_pipeline(**kwargs)
    report2 = run_pipeline(**kwargs)
    assert report1 == report2
    html = Path(report2).read_text()
    assert "GILGFVFTL" in html and "reliable" in html


def test_unannotatable_selection_does_not_crash(tmp_path):
    # top_n exceeds the number of annotatable clonotypes: select_top will
    # include an unannotatable clonotype (hla=None), which must be filtered
    # out before build_construct/fold rather than raising a KeyError on
    # mhc_seqs[None].
    report = run_pipeline(
        csv_path=FIX / "tenx_tiny.csv", run_dir=tmp_path / "run",
        top_n=10, sim_fn=_sim, assign_fn=_assign_fn, fold_fn=_fold_fn,
        tcr_seqs=None, mhc_seqs={"HLA-A*02:01": {"heavy": "H" * 20, "b2m": "M" * 20}},
        scramble_threshold=1.0)
    html = Path(report).read_text()
    assert "GILGFVFTL" in html
