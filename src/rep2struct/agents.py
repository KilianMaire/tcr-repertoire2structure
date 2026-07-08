from __future__ import annotations
from claude_agent_sdk import ClaudeAgentOptions, AgentDefinition
from .agent_tools import build_server
from . import structure_tools

_EXEC_TOOLS = ["mcp__rep2struct__list_fold_jobs", "mcp__rep2struct__build_fold_notebook",
               "mcp__rep2struct__record_fold_result", "mcp__playwright__*"]


def _executor(name, tool):
    return AgentDefinition(
        description=f"{name}: folds the {tool} group by driving its Colab notebook through the browser.",
        prompt=(
            f"You run the {tool} structure tool for the jobs assigned to your group. "
            f"Call list_fold_jobs, and for each job whose tool is '{tool}': call "
            f"build_fold_notebook with tool='{tool}' to write that job's self-contained "
            f".ipynb and get its path, then drive it with the mcp__playwright tools: "
            f"upload the notebook to Colab, connect a GPU runtime, run the cells (inputs "
            f"are embedded), wait, download the produced model, then call "
            f"record_fold_result with tool='{tool}' and the downloaded path. The loop is "
            f"resumable; skip jobs already recorded. If build_fold_notebook returns a "
            f"fail-loud scaffold (the {tool} adapter is not yet wired) or the Colab run "
            f"errors, report the job as not-run. Never fabricate a model or a score."),
        tools=list(_EXEC_TOOLS),
        model="sonnet",
    )


def build_agents():
    agents = {
        "structure-strategist": AgentDefinition(
            description="Reasons over the tool registry and construct tags; routes each group to a tool.",
            prompt=(
                "You choose structure tools. Call list_structure_tools to read each tool's "
                "validity domain, and list_fold_jobs to see the jobs (each carries group_id, "
                "mhc_class, has_tcr, species, output_needed). For each homogeneous group pick "
                "ONE tool: Protenix is the default workhorse; switch to a specialized tool "
                "only when the group justifies it (af3 if the user has weights and it helps, "
                "affinetune for is-it-presented in class I or II, mhcfine for a precise class I "
                "pose, tcrdock for the TCR interface). Never Boltz. If no tool's validity "
                "domain covers a group, fall back to Protenix and state plainly that an "
                "un-wired tool would fit better. Justify each choice in one sentence, then "
                "delegate the group to that tool's executor agent (protenix-agent, af3-agent, "
                "mhcfine-agent, tcrdock-agent, affinetune-agent)."),
            tools=["mcp__rep2struct__list_structure_tools", "mcp__rep2struct__list_fold_jobs", "Agent"],
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
                "recognition. Report reliable only when the CDR3 to peptide contact beats the "
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
        agents[f"{t.name}-agent"] = _executor(f"{t.name}-agent", t.name)
    return agents


def build_options(run_dir):
    return ClaudeAgentOptions(
        mcp_servers={
            "rep2struct": build_server(),
            "playwright": {"command": "npx", "args": ["@playwright/mcp@latest"]},
        },
        agents=build_agents(),
        allowed_tools=[
            "Agent",
            "mcp__rep2struct__ingest_repertoire", "mcp__rep2struct__annotate_specificity",
            "mcp__rep2struct__prep_and_select", "mcp__rep2struct__list_structure_tools",
            "mcp__rep2struct__list_fold_jobs", "mcp__rep2struct__build_fold_notebook",
            "mcp__rep2struct__record_fold_result",
            "mcp__rep2struct__qc_structure", "mcp__rep2struct__render_final_report",
            "mcp__playwright__*",
        ],
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
