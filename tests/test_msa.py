from rep2struct.schema import FoldJob
from rep2struct import msa


def _job():
    return FoldJob(clonotype_id="c1", construct_fasta=">A\nAAAA")


def test_local_runner_wins_and_caches_a3m(tmp_path):
    ref, basis = msa.build_msa(_job(), tmp_path, local_runner=lambda f: "A3M-LOCAL")
    assert basis == "local"
    assert open(ref).read() == "A3M-LOCAL"


def test_falls_back_to_colab_when_local_fails(tmp_path):
    def boom(f): raise RuntimeError("no local DB")
    ref, basis = msa.build_msa(_job(), tmp_path, local_runner=boom,
                               colab_runner=lambda f: "A3M-COLAB")
    assert basis == "colab_cpu"
    assert open(ref).read() == "A3M-COLAB"


def test_falls_back_to_msa_free_when_both_fail(tmp_path):
    def boom(f): raise RuntimeError("down")
    ref, basis = msa.build_msa(_job(), tmp_path, local_runner=boom, colab_runner=boom)
    assert basis == "none" and ref == ""
