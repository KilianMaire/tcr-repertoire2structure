# Intake agent and plastic compute routes

Date: 2026-07-10
Status: approved (sections 1 and 2 approved by user, spec pending review)
Branch: `feat/intake-agent-compute-routes`, created off `main` to stay isolated from the in flight benchmark branch

## Problem

R2S today runs one shot and non interactively. `app.run(csv_path, run_dir, top_n)` feeds a fixed CSV path to `orchestrator_prompt`, then `query()` streams a fully autonomous multi agent pass under `bypassPermissions`, and the tool executors drive Colab themselves through Playwright.

That does not match how the project is actually run. On live folds the user drives the Colab cells (see the recorded practice in the project memory), and the user also has a local H100 in addition to the Colab A100, so a browser driven Colab notebook is not always the right target. There is no way to tell R2S "here is my data, here is my question, here is the machine I can run predictions on" and have it adapt.

## Goal

Add a conversational intake phase in front of the existing orchestration, built as the plastic core of the package. When the user launches R2S, an intake agent asks the right questions, branching by the answers:

1. What kind of data.
2. Where the input document is (a path, or a file dropped into a watched folder).
3. What question or task the user wants answered.
4. What compute environment is available to run the predictions. If the user does not know, the agent proposes the simplest option for the case.

The agent then ingests these answers and proposes the downstream tasks to the following agents (strategist, executors, QC, report). Execution of the fold stays user driven: the orchestrator builds the right artifact for the chosen route and stops. The user folds, then reruns R2S on the same run directory and a checkpoint resume carries the run through QC and report.

## Non goals for v1

1. No agent driven SSH execution. The intake collects SSH details so the interview is complete and plastic, but the SSH runner is a stubbed, honest fallback that hands the user a job script. Wiring real SSH execution is a later, separate piece.
2. No change to the existing autonomous Playwright path. It stays intact and demoable. The new behavior is additive.
3. No new structure tools. Protenix stays the default workhorse; the existing `structure_tools` registry and QC are reused unchanged.

## Architecture

```
r2s (CLI, ClaudeSDKClient multi turn)
  |
  |-- PHASE 0  intake-agent (opus, plastic)
  |       asks the right questions, branched by answer:
  |         data type ? . where is the input ? . your question/task ?
  |         which environment ?  (unknown -> recommend the simplest, Colab)
  |       writes a non secret IntakeSpec into run_dir
  |
  |-- PHASE A  orchestrator (parameterized by the IntakeSpec)
  |       ingest -> annotate -> prep -> strategist routes each group
  |       -> executors in HANDOFF mode: build the artifact (notebook OR bash script), STOP
  |       -> propose then confirm the plan before generating artifacts
  |       [ user folds: Colab driven by the user, or the local H100 ]
  |
  |-- PHASE B  resume (rerun r2s on the same run_dir)
        checkpoint sees the done jobs / detects the CIFs -> record -> qc-agent -> report-agent
```

Two new modules (`intake.py`, `compute_routes.py`), one new agent (`intake-agent`), an interactive CLI, and a handoff mode on the executors (build the artifact, no Playwright). The prior full auto Playwright path is left intact.

## Section 2: the plastic core, `compute_routes` (open registry)

Same pattern as `structure_tools.py`: an open registry that is honest about what is wired. Each route declares the questions it requires and the artifact it produces.

| route | branched questions | artifact | runner wired in v1 |
|---|---|---|---|
| `colab` (default, "the simplest") | none extra | `.ipynb` (the `_protenix_notebook` recipe) | yes (user drives) |
| `local_gpu` (the H100) | working path on the machine | bash script running Protenix (`protenix_inputs.build` without the notebook shell) | yes (user runs) |
| `ssh` | host, user, SSH key first then password as last resort (ephemeral, never on disk), remote_path | bash script plus `scp`/`sbatch` instructions | no (honest stubbed runner) |
| `server` | address, path | same as ssh | no |

The agent reads this registry through `list_compute_routes`, so adding a route means adding one entry, not rewriting the interview. That is the plasticity requirement. Secrets never reach disk: the persisted `IntakeSpec` keeps only host, user, and path.

### ComputeRoute data shape

A pure data record, mirroring the `StructureTool` registry entry:

```
ComputeRoute(
  name: str,                 # "colab" | "local_gpu" | "ssh" | "server"
  description: str,          # one line the agent reads to reason
  required_fields: list[str],# e.g. ["host","user","remote_path"] ; [] for colab
  secret_fields: list[str],  # e.g. ["password"] ; never persisted
  artifact_kind: str,        # "colab_notebook" | "bash_script"
  wired: bool,               # True colab/local_gpu, False ssh/server
  is_default: bool,          # True for colab
)
```

Helpers, tested only, mirror `structure_tools`: `get_default()`, `by_name(name)`, `recommend(context)` (returns the default when the user does not know). Honesty rule: when a route is not wired, any downstream artifact builder must produce the fallback bash script and the executor must state plainly that the route runner is not wired, never fabricate a remote run.

## Section 3: IntakeSpec and intake MCP tools

`IntakeSpec` is a pure dataclass:

```
IntakeSpec(
  data_type: str,
  input_path: str,
  question: str,
  compute_route: str,
  route_params: dict,   # non secret only: host, user, path, working_path...
)
```

Round trips to `run_dir/intake.json`. Secret fields collected during the interview (an SSH password) live only in session memory and are never written into this file or any log.

New MCP tools, exposed to the intake agent only (no Playwright in its tool set):

1. `list_compute_routes()` returns the registry entries so the agent reasons over required fields and wired status.
2. `recommend_route(context)` returns the simplest route for the described case (the default, Colab) when the user does not know.
3. `record_intake(run_dir, spec)` validates and persists the non secret `IntakeSpec` to `run_dir/intake.json`.

The orchestrator then reads the IntakeSpec to parameterize `orchestrator_prompt`: `input_path` becomes the ingest CSV, `question` steers the strategist, and `compute_route` selects which artifact the executors build.

## Section 4: handoff mode on the orchestration

`build_options(run_dir, mode="handoff")` gains a mode. In handoff mode the executor prompt changes: for each pending job it calls the artifact builder for the chosen route and records the artifact path, then stops. It does not call any `mcp__playwright__*` tool. The strategist still routes groups to tools exactly as today.

Artifact dispatch: extend the existing `build_fold_notebook` seam into a route aware builder. For `artifact_kind == "colab_notebook"` it returns the `.ipynb` as today. For `artifact_kind == "bash_script"` it writes a self contained shell script built from `protenix_inputs.build` plus the proven Protenix invocation from `_protenix_notebook`, without the notebook shell, so the CIFs land directly under `run_dir/out/{cid}_cognate` and `{cid}_scramble` on the local machine. Unwired routes (ssh, server) reuse the bash script plus a printed `scp`/`sbatch` instruction block, and the executor reports the route runner as not wired.

Propose then confirm: after intake, before generating artifacts, the orchestrator prints its plan (which groups route to which tool on which compute route, and which artifact files it will write) and waits for the user to confirm in the CLI.

## Section 5: resume by checkpoint

Resume reuses the existing on disk checkpoint. `list_fold_jobs` marks a job done once its result is recorded. On a rerun over the same `run_dir`:

1. The `local_gpu` bash route wrote CIFs directly into `run_dir/out`, so they are already on disk. The resume path detects them and records them.
2. The `colab` route produced a downloaded zip. The resume interview asks the user where the download is, unzips it into `run_dir/out` keeping the `{cid}_cognate` and `{cid}_scramble` layout, and records the paths.

Once folds are recorded, the run proceeds to the `qc-agent` (per group scramble calibration, unchanged) and the `report-agent`.

## Testing

All deterministic layers are offline testable, no network and no live fold required to merge:

1. `compute_routes`: required fields per route, `is_default`, `recommend` falling back to the default, and the honesty of the ssh/server stub (not wired).
2. `IntakeSpec`: round trip to `run_dir/intake.json` with no secret field ever written (assert a collected password does not appear in the file).
3. Artifact dispatch: route to artifact_kind mapping produces a notebook for colab and a bash script for local_gpu, built without network, with the correct `{cid}_cognate` and `{cid}_scramble` output layout.
4. Resume: jobs marked done cause the fold stage to be skipped and the run to advance to QC.

The interview itself (ClaudeSDKClient multi turn) is validated in a manual run. No live fold is required to merge this layer.

## CLI entry point

`python -m rep2struct` (or a `r2s` console script) opens a `ClaudeSDKClient` session. The agent asks a question, reads the user answer from stdin, loops through the intake, then runs Phase A. A rerun on the same `run_dir` enters Phase B. The prior non interactive `app.run` stays available for the autonomous demo.

## Open decisions deferred to the plan

1. Exact console script name (`r2s` vs `python -m rep2struct`).
2. Whether `recommend_route` is a separate tool or folded into the agent prompt reasoning over `list_compute_routes`.
3. Whether the resume interview is a distinct small agent or the same intake agent detecting an existing `run_dir/intake.json` plus pending recorded folds.
