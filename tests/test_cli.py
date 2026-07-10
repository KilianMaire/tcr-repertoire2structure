from rep2struct import cli, intake
from rep2struct.intake import IntakeSpec


def test_plan_selects_intake_then_run(tmp_path):
    assert cli.plan_from_run_dir(str(tmp_path)) == "intake"
    intake.save_intake(str(tmp_path), IntakeSpec("d", "/i.csv", "q", "colab", {}))
    assert cli.plan_from_run_dir(str(tmp_path)) == "run"
