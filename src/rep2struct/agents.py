from __future__ import annotations
from claude_agent_sdk import ClaudeAgentOptions, AgentDefinition
from .agent_tools import build_server
from . import structure_tools

_EXEC_TOOLS = ["mcp__rep2struct__list_fold_jobs", "mcp__rep2struct__build_fold_notebook",
               "mcp__rep2struct__record_fold_result", "mcp__playwright__*"]


_HANDOFF_EXEC_TOOLS = ["mcp__rep2struct__list_fold_jobs",
                       "mcp__rep2struct__build_group_artifact"]


def _handoff_executor(name, tool):
    return AgentDefinition(
        description=f"{name}: builds ONE {tool} fold artifact per group and stops "
                    f"for the user to run it.",
        prompt=(
            f"You prepare the {tool} folds for the jobs assigned to your group; you do NOT "
            f"run them. Call list_fold_jobs. For each DISTINCT group_id whose jobs are "
            f"tool='{tool}' and not done, call build_group_artifact ONCE with that group_id, "
            f"tool='{tool}', and the run's compute_route: it writes a SINGLE self-contained "
            f"artifact (a Colab notebook or a local bash script, chosen by the route) that "
            f"folds every pending clonotype in the group in one run, so the user uploads/runs "
            f"one artifact per group, not one per TCR. Report each artifact path with one line "
            f"on how to run it: a Colab notebook gets uploaded and its cells run; a bash script "
            f"runs on the target machine. If build_group_artifact reports route_wired false, "
            f"say plainly the route runner is not wired and hand over the script for the user "
            f"to run. Never fabricate a model or a score, and never open a browser."),
        tools=list(_HANDOFF_EXEC_TOOLS),
        model="sonnet",
    )


def _executor(name, tool, mode="auto"):
    if mode == "handoff":
        return _handoff_executor(name, tool)
    return AgentDefinition(
        description=f"{name}: folds the {tool} group by driving its Colab notebook through the browser.",
        prompt=(
            f"You run the {tool} structure tool for the jobs assigned to your group. "
            f"Call list_fold_jobs, and for each job whose tool is '{tool}' AND whose done is false: call "
            f"build_fold_notebook with tool='{tool}' to write that job's self-contained "
            f".ipynb and get its path, then drive it with the mcp__playwright tools: "
            f"upload the notebook to Colab, connect a GPU runtime, then run the cells by "
            f"KEYBOARD SHORTCUT (focus a cell and press Ctrl+Enter, or Shift+Enter to run "
            f"and advance through all cells) because the Run all toolbar button does not "
            f"fire reliably through Playwright. Wait for the fold to finish, download the "
            f"produced model, then call record_fold_result with tool='{tool}' and the "
            f"downloaded path. For a score-based tool (tcrdock, affinetune) the output is the "
            f"per-construct .score files; the notebook's LAST cell zips them and calls "
            f"files.download('{tool}_scores.zip'). Capture that zip, UNZIP it into the run dir so "
            f"{{cid}}_cognate.score and {{cid}}_scramble.score land in the SAME folder, and record "
            f"the cognate .score path, so QC reads the group's own scramble null next to it. "
            f"For Protenix (a structure tool folded as a "
            f"cognate+scramble pair, {tool}=='protenix'), the notebook writes CIFs under "
            f"{{cid}}_cognate and {{cid}}_scramble, and its LAST cell repatriates them: it zips the "
            f"CIFs (keeping the out/{{cid}}_cognate and out/{{cid}}_scramble layout) and calls "
            f"files.download('protenix_folds.zip'). Capture that browser download, UNZIP it into the "
            f"run dir so every extracted CIF path still carries {{cid}}_cognate or {{cid}}_scramble, "
            f"then record ALL of them in one record_fold_result call, so QC can ensemble the cognate "
            f"samples and calibrate against the scramble null. Do not consider the job done until the "
            f"CIFs are actually on disk locally. The loop is resumable by checkpoint: list_fold_jobs "
            f"marks done=true once a job's result is recorded (the folds stage is the on-disk save that "
            f"survives a crash between clonotypes), so a rerun re-folds only the pending jobs. If "
            f"build_fold_notebook returns a fail-loud scaffold (the {tool} adapter is not "
            f"yet wired) or the Colab run errors, report the job as not-run and never "
            f"fabricate a model or a score. ALWAYS, once every job for your group is folded "
            f"and recorded (or on a fatal error), release the GPU: open the Runtime menu and "
            f"choose 'Disconnect and delete runtime' so no runtime is left running."),
        tools=list(_EXEC_TOOLS),
        model="sonnet",
    )


def build_agents(mode="auto"):
    agents = {
        "intake-agent": AgentDefinition(
            description="Conversational intake: asks data, question, and compute environment, "
                        "then hands a structured brief to the orchestrator.",
            prompt=(
                "You run the R2S intake interview. Ask the user, one question at a time and "
                "branching on their answers: (1) what kind of data they have, (2) where the "
                "input file is (a path, or a file they dropped into the run folder), (3) what "
                "question or task they want answered, (4) which compute environment they can "
                "run predictions on. For (4) call list_compute_routes and collect EVERY "
                "required_field of the route the user names, never skipping one even if the "
                "route seems to 'need no setup': local_gpu requires working_path (the "
                "directory on that machine where the fold outputs must be written); ssh "
                "requires host, user, an SSH key first with a password only as a last resort, "
                "and remote_path; server requires address and path; colab requires nothing "
                "extra. If the user does "
                "not know, propose the default route (Colab), the simplest option. Never "
                "write a password anywhere. When you have all four, summarize the brief back "
                "and confirm before the run proceeds. When you have all four answers, call "
                "record_intake to persist the brief, then confirm before the run proceeds."),
            tools=["mcp__rep2struct__list_compute_routes",
                   "mcp__rep2struct__ingest_repertoire",
                   "mcp__rep2struct__record_intake", "Agent"],
            model="opus",
        ),
        "structure-strategist": AgentDefinition(
            description="Reasons over the tool registry and construct tags; routes each group to a tool.",
            prompt=(
                "You choose structure tools. Call list_structure_tools to read each tool's "
                "validity domain, and list_fold_jobs to see the jobs (each carries group_id, "
                "mhc_class, has_tcr, species, output_needed). For each homogeneous group pick "
                "ONE tool. Protenix is the default and it folds the FULL TCR-pMHC (all chains "
                "A-E), so has_tcr=True does NOT on its own justify a switch: Protenix already "
                "handles TCR groups. Start from Protenix and keep it unless the group gives a "
                "positive, specific reason to switch, which you must name:\n"
                "- af3 only if the user has the gated weights and they help;\n"
                "- affinetune when the question is literally is-this-peptide-presented "
                "(output_needed=binding_score), class I or II;\n"
                "- mhcfine only for a precise class I pose with NO TCR (needs_tcr=False);\n"
                "- tcrdock only when the interface PAE is the actual point AND you can say why "
                "Protenix's full fold plus beta V-domain to peptide QC is insufficient for this group; never "
                "pick tcrdock merely because a TCR is present, and never justify it by claiming "
                "Protenix cannot model a TCR (it can).\n"
                "Never Boltz. If no tool's validity domain covers a group, fall back to Protenix "
                "and state plainly that an un-wired tool would fit better. Justify each choice in "
                "one sentence that names the group field driving it. Then, for each group, call "
                "assign_group_tool(group_id, tool) to persist the choice onto its jobs BEFORE "
                "delegating, so the executor filters on the tag by fact and the report shows the "
                "right tool per group. Finally delegate the group to that tool's executor agent "
                "(protenix-agent, af3-agent, mhcfine-agent, tcrdock-agent, affinetune-agent)."),
            tools=["mcp__rep2struct__list_structure_tools", "mcp__rep2struct__list_fold_jobs",
                   "mcp__rep2struct__assign_group_tool", "Agent"],
            model="opus",
        ),
        "qc-agent": AgentDefinition(
            description="Skeptical QC per group; calibration is per tool, output-type aware.",
            prompt=(
                "You are a skeptical structural referee. For each folded clonotype call "
                "qc_structure with that group's OWN scramble_threshold (never a global one) "
                "and its output_type and tool. A clean fold does NOT confirm specificity. For "
                "binding_score tools you judge predicted presentation, not geometry. For a "
                "peptide_groove (pose) tool like mhcfine the verdict is pose-only: it seats "
                "any peptide in the groove, so a pose is placement, never proof of "
                "recognition. Report reliable only when the beta V-domain to peptide contact beats the "
                "group's scramble calibration; otherwise suspect. Never upgrade a verdict to "
                "please the caller."),
            tools=["mcp__rep2struct__qc_structure"],
            model="opus",
        ),
        "report-agent": AgentDefinition(
            description="Renders the final HTML report tying clonotype to specificity to structure to QC.",
            prompt="Call render_final_report and return the report path. Add no unsupported claims.",
            tools=["mcp__rep2struct__render_final_report"],
            model="sonnet",
        ),
    }
    for t in structure_tools.REGISTRY:
        agents[f"{t.name}-agent"] = _executor(f"{t.name}-agent", t.name, mode=mode)
    return agents


def build_options(run_dir, mode="auto"):
    base = [
        "Agent",
        "mcp__rep2struct__ingest_repertoire", "mcp__rep2struct__annotate_specificity",
        "mcp__rep2struct__prep_and_select", "mcp__rep2struct__list_structure_tools",
        "mcp__rep2struct__list_fold_jobs", "mcp__rep2struct__assign_group_tool",
        "mcp__rep2struct__list_compute_routes",
        "mcp__rep2struct__record_fold_result", "mcp__rep2struct__record_local_folds",
        "mcp__rep2struct__qc_structure", "mcp__rep2struct__render_final_report",
        "mcp__rep2struct__record_intake",
    ]
    if mode == "handoff":
        allowed = base + ["mcp__rep2struct__build_fold_artifact",
                          "mcp__rep2struct__build_group_artifact"]
        servers = {"rep2struct": build_server()}
    else:
        allowed = base + ["mcp__rep2struct__build_fold_notebook", "mcp__playwright__*"]
        servers = {"rep2struct": build_server(),
                   "playwright": {"command": "npx", "args": ["@playwright/mcp@latest"]}}
    return ClaudeAgentOptions(
        mcp_servers=servers,
        agents=build_agents(mode=mode),
        allowed_tools=allowed,
        # Autonomous multi-agent demo run; no interactive approval needed during live orchestration.
        permission_mode="bypassPermissions",
    )


def orchestrator_prompt(csv_path, run_dir, top_n, question=None):
    steer = (f"\nUser question steering the routing: {question}\n" if question else
             "\nNo user question: route each group from its construct tags.\n")
    return (
        f"Run the repertoire to structure pipeline on {csv_path} with run_dir {run_dir}.{steer}"
        f"1. Call ingest_repertoire, then annotate_specificity (honest annotation, keep unannotatable as is).\n"
        f"2. Call prep_and_select with top_n {top_n}.\n"
        f"3. Delegate to the structure-strategist to route each group to a tool and drive its executor to fold.\n"
        f"4. Delegate to the qc-agent to QC each folded clonotype with that group's calibration.\n"
        f"5. Delegate to the report-agent to render the final HTML report, and return its path.")


def intake_orchestrator_prompt(run_dir, spec, top_n=8):
    return (
        f"Run the repertoire to structure pipeline on {spec.input_path} with run_dir "
        f"{run_dir}.\nUser question steering the routing: {spec.question}\n"
        f"Compute route: {spec.compute_route} (route params {spec.route_params}).\n"
        f"1. Call ingest_repertoire, then annotate_specificity (honest annotation, keep "
        f"unannotatable as is).\n"
        f"2. Call prep_and_select with top_n {top_n}.\n"
        f"3. Delegate to the structure-strategist to route each group to a tool; each tool's "
        f"executor builds the fold artifact for compute_route '{spec.compute_route}' and "
        f"STOPS. Present the plan (which group to which tool, which artifact files) and let "
        f"the user confirm before generating artifacts.\n"
        f"4. Tell the user to run the artifacts and rerun R2S on this run_dir to resume.")
