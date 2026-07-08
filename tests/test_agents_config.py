from rep2struct.agents import build_agents, build_options, orchestrator_prompt


def test_strategist_and_executors_present():
    agents = build_agents()
    assert "structure-strategist" in agents
    assert {"protenix-agent", "af3-agent", "mhcfine-agent",
            "tcrdock-agent", "affinetune-agent"} <= set(agents)
    assert {"qc-agent", "report-agent"} <= set(agents)


def test_strategist_reads_registry_and_can_delegate():
    a = build_agents()["structure-strategist"]
    assert any("list_structure_tools" in t for t in a.tools)
    assert "Agent" in a.tools


def test_strategist_exact_tools_and_model():
    a = build_agents()["structure-strategist"]
    assert set(a.tools) == {"mcp__rep2struct__list_structure_tools",
                            "mcp__rep2struct__list_fold_jobs", "Agent"}
    assert a.model == "opus"


def test_executors_can_drive_browser_and_record():
    ex = build_agents()["protenix-agent"]
    assert any("playwright" in t for t in ex.tools)
    assert any("record_fold_result" in t for t in ex.tools)


def test_every_executor_exact_tools_and_model():
    agents = build_agents()
    for name in ("protenix-agent", "af3-agent", "mhcfine-agent",
                 "tcrdock-agent", "affinetune-agent"):
        ex = agents[name]
        assert ex.model == "sonnet"
        assert set(ex.tools) == {"mcp__rep2struct__list_fold_jobs",
                                 "mcp__rep2struct__record_fold_result",
                                 "mcp__playwright__*"}


def test_options_wire_new_tool_and_agents(tmp_path):
    opts = build_options(str(tmp_path / "run"))
    assert any("list_structure_tools" in t for t in opts.allowed_tools)
    assert "structure-strategist" in opts.agents
    assert opts.permission_mode == "bypassPermissions"


def test_prompt_supports_optional_question(tmp_path):
    base = orchestrator_prompt("x.csv", str(tmp_path), 8)
    for kw in ["ingest", "annotate", "fold", "qc", "report"]:
        assert kw in base.lower()
    steered = orchestrator_prompt("x.csv", str(tmp_path), 8, question="which clones are presented in class II")
    assert "class ii" in steered.lower()
