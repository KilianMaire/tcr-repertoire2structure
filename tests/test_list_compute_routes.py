import asyncio
from rep2struct import agent_tools


def test_lists_the_routes_with_required_and_wired_flags(tmp_path):
    r = asyncio.run(agent_tools.list_compute_routes.handler({"run_dir": str(tmp_path)}))
    routes = {x["name"]: x for x in r["structuredContent"]["routes"]}
    assert set(routes) == {"colab", "local_gpu", "ssh", "server"}
    assert routes["colab"]["is_default"] is True
    assert routes["ssh"]["required_fields"] == ["host", "user", "remote_path"]
    assert routes["ssh"]["wired"] is False


def test_text_carries_fields_so_a_stub_reader_still_sees_them(tmp_path):
    # An agent that only sees the text content (not structuredContent) must still learn a
    # route's required fields; colab's empty list must read as "none", never as unknown.
    r = asyncio.run(agent_tools.list_compute_routes.handler({"run_dir": str(tmp_path)}))
    text = r["content"][0]["text"]
    assert "colab" in text and "required_fields=[none]" in text
    assert "host, user, remote_path" in text
    assert "runner_wired=False" in text
