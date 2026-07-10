import asyncio
from rep2struct import agent_tools


def test_lists_the_routes_with_required_and_wired_flags(tmp_path):
    r = asyncio.run(agent_tools.list_compute_routes.handler({"run_dir": str(tmp_path)}))
    routes = {x["name"]: x for x in r["structuredContent"]["routes"]}
    assert set(routes) == {"colab", "local_gpu", "ssh", "server"}
    assert routes["colab"]["is_default"] is True
    assert routes["ssh"]["required_fields"] == ["host", "user", "remote_path"]
    assert routes["ssh"]["wired"] is False
