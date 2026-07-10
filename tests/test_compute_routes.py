import pytest

from rep2struct import compute_routes as cr


def test_default_is_colab_the_simplest():
    d = cr.get_default()
    assert d.name == "colab" and d.is_default


def test_registry_has_the_four_v1_routes():
    assert {r.name for r in cr.REGISTRY} == {"colab", "local_gpu", "ssh", "server"}


def test_colab_and_local_gpu_are_wired_ssh_and_server_are_not():
    assert cr.is_wired("colab") and cr.is_wired("local_gpu")
    assert not cr.is_wired("ssh") and not cr.is_wired("server")


def test_artifact_kind_per_route():
    assert cr.artifact_kind_for("colab") == "colab_notebook"
    assert cr.artifact_kind_for("local_gpu") == "bash_script"
    assert cr.artifact_kind_for("ssh") == "bash_script"
    assert cr.artifact_kind_for("server") == "bash_script"
    assert cr.by_name("local_gpu").required_fields == ("working_path",)
    assert cr.by_name("server").required_fields == ("address", "path")


def test_ssh_requires_connection_fields_and_marks_password_secret():
    ssh = cr.by_name("ssh")
    assert ssh.required_fields == ("host", "user", "remote_path")
    assert ssh.secret_fields == ("password",)
    # colab needs nothing extra
    assert cr.by_name("colab").required_fields == ()


def test_recommend_falls_back_to_default_when_user_does_not_know():
    assert cr.recommend("I don't know").name == "colab"
    assert cr.recommend("").name == "colab"


def test_as_dicts_exposes_fields_but_never_a_secret_value():
    d = {x["name"]: x for x in cr.as_dicts()}
    assert d["ssh"]["required_fields"] == ["host", "user", "remote_path"]
    assert d["ssh"]["secret_fields"] == ["password"]
    assert d["ssh"]["wired"] is False and d["colab"]["is_default"] is True


def test_by_name_unknown_raises_valueerror():
    with pytest.raises(ValueError):
        cr.by_name("nope")
