import asyncio
from pathlib import Path

from rep2struct import agent_tools, intake


def test_record_intake_persists_and_strips_secrets(tmp_path):
    rd = str(tmp_path)
    r = asyncio.run(agent_tools.record_intake.handler({
        "run_dir": rd,
        "data_type": "10x_vdj",
        "input_path": "/data/contigs.csv",
        "question": "which epitope?",
        "compute_route": "ssh",
        "route_params": {"host": "hpc.uni.dk", "user": "kilian",
                         "remote_path": "/work", "password": "hunter2"},
    }))
    sc = r["structuredContent"]
    assert sc["intake_path"].endswith("intake.json")
    assert sc["compute_route"] == "ssh"
    assert Path(sc["intake_path"]).exists()

    raw = Path(sc["intake_path"]).read_text()
    assert "hunter2" not in raw
    assert "password" not in raw

    got = intake.load_intake(rd)
    assert got.data_type == "10x_vdj"
    assert got.input_path == "/data/contigs.csv"
    assert got.question == "which epitope?"
    assert got.compute_route == "ssh"
    assert got.route_params["host"] == "hpc.uni.dk"
    assert "password" not in got.route_params
