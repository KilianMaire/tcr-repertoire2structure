from rep2struct import intake
from rep2struct.intake import IntakeSpec


def test_round_trip(tmp_path):
    spec = IntakeSpec("10x_vdj", "/data/contigs.csv", "which epitope?",
                      "local_gpu", {"working_path": "/scratch/run"})
    p = intake.save_intake(str(tmp_path), spec)
    assert p.endswith("intake.json")
    got = intake.load_intake(str(tmp_path))
    assert got == spec


def test_secret_never_written_to_disk(tmp_path):
    # A password slipped into route_params must NOT reach intake.json.
    spec = IntakeSpec("10x_vdj", "/data/c.csv", "q", "ssh",
                      {"host": "hpc.uni.dk", "user": "kilian",
                       "remote_path": "/work", "password": "hunter2"})
    p = intake.save_intake(str(tmp_path), spec)
    with open(p) as fh:
        raw = fh.read()
    assert "hunter2" not in raw
    assert "password" not in raw
    got = intake.load_intake(str(tmp_path))
    assert "password" not in got.route_params
    assert got.route_params["host"] == "hpc.uni.dk"


def test_next_phase_is_intake_when_no_file_then_run(tmp_path):
    assert intake.next_phase(str(tmp_path)) == "intake"
    intake.save_intake(str(tmp_path),
                       IntakeSpec("d", "/i.csv", "q", "colab", {}))
    assert intake.next_phase(str(tmp_path)) == "run"


def test_load_returns_none_when_absent(tmp_path):
    assert intake.load_intake(str(tmp_path)) is None
