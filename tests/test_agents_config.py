from rep2struct.agents import build_agents, build_options, orchestrator_prompt


def test_specialist_agents_present():
    agents = build_agents()
    assert set(agents) >= {"fold-agent", "qc-agent", "report-agent"}
    # the fold agent must be allowed to drive the browser
    assert any("playwright" in t for t in agents["fold-agent"].tools)


def test_options_wire_tools_and_agents(tmp_path):
    opts = build_options(str(tmp_path / "run"))
    assert "rep2struct" in opts.mcp_servers
    assert "playwright" in opts.mcp_servers
    assert any(t.startswith("mcp__rep2struct__") or t == "Agent" for t in opts.allowed_tools)
    assert set(opts.agents) >= {"fold-agent", "qc-agent", "report-agent"}
    assert opts.permission_mode == "bypassPermissions"


def test_prompt_names_the_stages(tmp_path):
    p = orchestrator_prompt("x.csv", str(tmp_path), 8)
    for kw in ["ingest", "annotate", "fold", "qc", "report"]:
        assert kw in p.lower()
