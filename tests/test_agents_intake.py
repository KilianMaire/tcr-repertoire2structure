from rep2struct import agents


def test_intake_agent_exists_without_playwright():
    a = agents.build_agents()
    assert "intake-agent" in a
    assert not any("playwright" in t for t in a["intake-agent"].tools)


def test_handoff_options_drop_playwright_and_add_artifact_tool():
    opts = agents.build_options("/tmp/run", mode="handoff")
    assert not any("playwright" in t for t in opts.allowed_tools)
    assert "mcp__rep2struct__build_fold_artifact" in opts.allowed_tools
    assert "mcp__rep2struct__record_local_folds" in opts.allowed_tools


def test_auto_options_keep_playwright_unchanged():
    opts = agents.build_options("/tmp/run")  # default mode="auto"
    assert any("playwright" in t for t in opts.allowed_tools)


def test_handoff_executor_prompt_builds_artifact_and_does_not_drive_colab():
    ex = agents._executor("protenix-agent", "protenix", mode="handoff")
    assert "build_group_artifact" in ex.prompt  # one artifact per group, batched
    assert "playwright" not in ex.prompt.lower()
    assert "Ctrl+Enter" not in ex.prompt


def test_intake_orchestrator_prompt_threads_the_spec():
    from rep2struct.intake import IntakeSpec
    spec = IntakeSpec("10x", "/data/c.csv", "which epitope?", "local_gpu",
                      {"working_path": "/scratch"})
    p = agents.intake_orchestrator_prompt("/tmp/run", spec)
    assert "/data/c.csv" in p and "which epitope?" in p and "local_gpu" in p
    assert "top_n 8" in p  # default selection depth


def test_intake_orchestrator_prompt_respects_top_n():
    from rep2struct.intake import IntakeSpec
    spec = IntakeSpec("10x", "/data/c.csv", "q?", "colab", {})
    p = agents.intake_orchestrator_prompt("/tmp/run", spec, top_n=25)
    assert "top_n 25" in p and "top_n 8" not in p
