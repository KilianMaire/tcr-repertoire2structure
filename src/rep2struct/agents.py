from __future__ import annotations
from claude_agent_sdk import ClaudeAgentOptions, AgentDefinition
from .agent_tools import build_server


def build_agents():
    return {
        "fold-agent": AgentDefinition(
            description="Folds TCR pMHC constructs by driving Protenix on Colab through the browser.",
            prompt=(
                "You fold structures. Call list_fold_jobs to get constructs. For each job, "
                "drive the Protenix Colab notebook with the mcp__playwright tools: open the "
                "notebook, submit the construct FASTA at 5 seeds, wait for completion, download "
                "the CIF models, then call record_fold_result with the local model paths. The "
                "loop is resumable; skip a job that already has recorded models."),
            tools=["mcp__rep2struct__list_fold_jobs", "mcp__rep2struct__record_fold_result",
                   "mcp__playwright__*"],
            model="sonnet",
        ),
        "qc-agent": AgentDefinition(
            description="Skeptical QC of predicted TCR pMHC structures. Flags geometry hallucinations.",
            prompt=(
                "You are a skeptical structural referee. For each folded clonotype call "
                "qc_structure. A clean fold does NOT confirm specificity: Protenix imposes "
                "canonical docking geometry even on non binders. Report reliable only when the "
                "CDR3 to peptide contact beats the scramble calibration; otherwise suspect. Never "
                "upgrade a verdict to please the caller."),
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
            "mcp__rep2struct__prep_and_select", "mcp__rep2struct__list_fold_jobs",
            "mcp__rep2struct__record_fold_result", "mcp__rep2struct__qc_structure",
            "mcp__rep2struct__render_final_report", "mcp__playwright__*",
        ],
        # Autonomous multi-agent demo run; no interactive approval needed during live orchestration.
        permission_mode="bypassPermissions",
    )


def orchestrator_prompt(csv_path, run_dir, top_n):
    return (
        f"Run the repertoire to structure pipeline on {csv_path} with run_dir {run_dir}.\n"
        f"1. Call ingest_repertoire, then annotate_specificity (honest annotation, keep unannotatable as is).\n"
        f"2. Call prep_and_select with top_n {top_n}.\n"
        f"3. Delegate to the fold-agent to fold the prepared jobs.\n"
        f"4. Delegate to the qc-agent to QC each folded clonotype (scramble_threshold from calibration).\n"
        f"5. Delegate to the report-agent to render the final HTML report, and return its path.")
