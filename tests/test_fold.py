from rep2struct.fold import run_folds
from rep2struct.schema import FoldJob
from rep2struct.runstate import RunState

def _job(cid): return FoldJob(clonotype_id=cid, construct_fasta=">A\nAAAA")

def test_resume_skips_done(tmp_path):
    rs = RunState(tmp_path / "r")
    calls = []
    def fold_fn(job):
        calls.append(job.clonotype_id)
        return [f"{job.clonotype_id}.cif"]
    run_folds([_job("c1")], fold_fn, rs)
    run_folds([_job("c1")], fold_fn, rs)  # second run resumes
    assert calls == ["c1"]  # folded once, skipped the second time

def test_failure_does_not_abort_batch(tmp_path):
    rs = RunState(tmp_path / "r2")
    def fold_fn(job):
        if job.clonotype_id == "bad":
            raise RuntimeError("colab wedged")
        return ["ok.cif"]
    out = run_folds([_job("bad"), _job("good")], fold_fn, rs)
    status = {j.clonotype_id: j.status for j in out}
    assert status["bad"] == "failed" and status["good"] == "done"
