from __future__ import annotations
from dataclasses import asdict
from pathlib import Path
from claude_agent_sdk import tool, create_sdk_mcp_server
from .runstate import RunState
from .ingest import parse_10x, standardize_alleles
from .annotate import annotate
from .foldprep import select_top, build_construct
from .seqs import build_tcr_seqs, build_mhc_seqs
from .qc import score_model, verdict
from .report import render_report
from .schema import Clonotype, Annotation, QCResult
from . import structure_tools
from . import compute_routes
from . import intake
from .grouping import partition
from .msa import build_msa

_CFG = {"sim_fn": None, "assign_fn": None}


def configure(sim_fn=None, assign_fn=None):
    """Inject offline fakes for the similarity search and allele assignment
    calls that stage functions would otherwise make against tcr_explorer /
    the network. Call with no arguments to reset to the defaults."""
    _CFG["sim_fn"] = sim_fn
    _CFG["assign_fn"] = assign_fn


def _txt(s):
    return {"content": [{"type": "text", "text": s}]}


def _scramble_null(cognate_score_path):
    """The group's OWN scramble score, read from the sibling file the fold writes next to the
    cognate one ({cid}_cognate.score and {cid}_scramble.score). That score IS the per-group
    binding-score null, so it sets the threshold (beat-the-null: a real cognate must score
    above its own shuffled peptide). Returns None if this is not a cognate score path or the
    scramble sibling was not downloaded, so the caller can fall back to an explicit threshold."""
    p = Path(cognate_score_path)
    if "_cognate" not in p.name:
        return None
    sib = p.with_name(p.name.replace("_cognate", "_scramble"))
    if not sib.exists():
        return None
    try:
        return float(sib.read_text().strip())
    except (ValueError, OSError):
        return None


def scan_recorded_folds(run_dir, tool="protenix"):
    """Find fold outputs already on disk under <run_dir>/out and group them per clonotype,
    so the handoff/local resume path (record_local_folds -> QC) works for every tool, not
    just Protenix. The on-disk shape is set by the tool's qc_metric: a structure tool leaves
    CIFs, a binding-score tool leaves a per-construct .score file, and the two are laid out
    differently, so dispatch on the metric. Resume path: the local_gpu bash route wrote the
    outputs here directly, and a Colab download unzipped here has the same layout."""
    out = Path(run_dir) / "out"
    if not out.exists():
        return {}
    if structure_tools.qc_metric_for(tool) == "binding_score":
        return _scan_score_folds(out, tool)
    return _scan_structure_folds(out, tool)


def _scan_structure_folds(out, tool):
    """Structure tools (Protenix, af3): the construct is marked by the DIRECTORY
    ({cid}_cognate / {cid}_scramble), not the filename, so parse the cid from the top-level
    out/ dir and collect the CIFs inside it."""
    found = {}
    for d in sorted(out.iterdir()):
        if not d.is_dir():
            continue
        for suffix in ("_cognate", "_scramble"):
            if d.name.endswith(suffix):
                cid = d.name[: -len(suffix)]
                cifs = sorted(str(p) for p in d.rglob("*.cif"))
                if cifs:
                    found.setdefault(cid, {"paths": [], "tool": tool})
                    found[cid]["paths"].extend(cifs)
    return found


def _scan_score_folds(out, tool):
    """Binding-score tools (tcrdock, affinetune): the construct is marked by the FILE
    ({cid}_cognate.score, with a {cid}_scramble.score sibling in the same folder), not by a
    _cognate/_scramble directory. Parse the cid from the filename and record the cognate
    score path per clonotype; QC's _scramble_null reads the scramble sibling next to it."""
    found = {}
    for p in sorted(out.rglob("*_cognate.score")):
        cid = p.name[: -len("_cognate.score")]
        found.setdefault(cid, {"paths": [], "tool": tool})
        found[cid]["paths"].append(str(p))
    return found


def _load(rd, name, cls):
    rs = RunState(rd)
    return [cls(**d) for d in rs.read_stage(name)] if rs.stage_done(name) else []


@tool("ingest_repertoire", "Parse a 10x contig CSV into paired clonotypes and persist them.",
      {"run_dir": str, "csv_path": str})
async def ingest_repertoire(args):
    if not Path(args["csv_path"]).exists():
        return _txt(f"input CSV not found: {args['csv_path']} (give the intake agent a real "
                    f"path, or drop the file into the run folder)")
    clons = standardize_alleles(parse_10x(args["csv_path"]), assign_fn=_CFG["assign_fn"])
    RunState(args["run_dir"]).write_stage("ingest", clons)
    r = _txt(f"{len(clons)} clonotypes ingested")
    r["structuredContent"] = {"clonotypes": len(clons)}
    return r


@tool("annotate_specificity", "Annotate persisted clonotypes with candidate epitopes by TCRdist. Never forces a label.",
      {"run_dir": str})
async def annotate_specificity(args):
    import os
    clons = _load(args["run_dir"], "ingest", Clonotype)
    # Scale guard: annotation is one TCRdist call per clonotype, so on a 20k-clonotype 10x sample
    # annotating everything before select_top is the ingest-time hang. Clonotypes are size-sorted
    # (parse_10x) and select_top ranks on tier*size, so annotating only the largest `cap` spares
    # the O(N) calls with essentially no change to which top_n get selected. The skipped tail is
    # reported, never silently dropped (raise R2S_ANNOTATE_CAP to annotate deeper).
    cap = int(os.environ.get("R2S_ANNOTATE_CAP", "3000"))
    head = clons[:cap]
    skipped = len(clons) - len(head)
    anns = annotate(head, sim_fn=_CFG["sim_fn"])
    RunState(args["run_dir"]).write_stage("annotate", anns)
    tiers = {}
    for a in anns:
        tiers[a.confidence_tier] = tiers.get(a.confidence_tier, 0) + 1
    msg = f"annotated {len(anns)} clonotypes: {tiers}"
    if skipped:
        msg += f" ({skipped} smaller clones beyond the top {cap} by size not annotated)"
    r = _txt(msg)
    r["structuredContent"] = {"tiers": tiers, "annotated": len(anns), "skipped": skipped}
    return r


@tool("prep_and_select", "Rank clonotypes and build class I fold constructs for the top N.",
      {"run_dir": str, "top_n": int})
async def prep_and_select(args):
    rs = RunState(args["run_dir"])
    # Checkpoint: once foldjobs exist, this stage is immutable (like list_fold_jobs treats the
    # folds stage). A resume run must NOT recompute it: doing so would wipe the strategist's
    # persisted per-group tool tags and the per-clonotype MSA, and a different top_n would
    # re-select a different set after folds already exist, silently dropping folded clonotypes
    # from the report. Return the existing selection unchanged.
    if rs.stage_done("foldjobs"):
        jobs = rs.read_stage("foldjobs")
        r = _txt(f"foldjobs already prepared ({len(jobs)} jobs); keeping the existing selection "
                 f"and the strategist's tool assignments (rerun in a fresh run_dir to re-select)")
        r["structuredContent"] = {"jobs": [j["clonotype_id"] for j in jobs], "reused": True}
        return r
    clons = _load(args["run_dir"], "ingest", Clonotype)
    anns = _load(args["run_dir"], "annotate", Annotation)
    top = select_top(clons, anns, n=args["top_n"])
    # Guard (same as Task 10's pipeline.run_pipeline): select_top can return
    # clonotypes whose annotation is not foldable (unannotatable, hla=None).
    # build_construct would KeyError on mhc_seqs[None] for those, so filter
    # to foldable entries before building constructs. Filtered-out
    # clonotypes simply get no fold job.
    foldable = [(c, a) for c, a in top if a.annotatable and a.hla]
    # Real reconstructed V domains and fetched/cached MHC ectodomains (same
    # providers as pipeline.run_pipeline). Clonotypes whose HLA cannot be
    # resolved to a heavy chain get no fold job.
    seqs = build_tcr_seqs([c for c, _ in foldable])
    mhc = build_mhc_seqs(sorted({a.hla for _, a in foldable}))
    # build_construct returns None for a class II allele (the class I A-E+b2m construct cannot
    # represent it), so filter those out rather than folding a wrong complex.
    jobs = [j for j in (build_construct(c, a, seqs, mhc)
                        for c, a in foldable if a.hla in mhc) if j is not None]
    partition(jobs)  # stamps group_id on each job in place
    for j in jobs:   # MSA is a pre-fold artifact; no runners here = MSA-free default
        j.msa_ref, j.msa_basis = build_msa(j, args["run_dir"])
    RunState(args["run_dir"]).write_stage("foldjobs", jobs)
    r = _txt(f"prepared {len(jobs)} fold jobs")
    r["structuredContent"] = {"jobs": [j.clonotype_id for j in jobs]}
    return r


@tool("list_fold_jobs", "List fold jobs and their construct FASTA for the fold agent.",
      {"run_dir": str})
async def list_fold_jobs(args):
    rs = RunState(args["run_dir"])
    jobs = rs.read_stage("foldjobs") if rs.stage_done("foldjobs") else []
    # Deterministic resume checkpoint: the folds stage is the on-disk record that survives a
    # crash between clonotypes. Mark each job done iff its result is already recorded, so the
    # executor skips it by fact, not by the LLM's judgement, and only folds the pending ones.
    done = rs.read_stage("folds") if rs.stage_done("folds") else {}
    for j in jobs:
        j["done"] = j["clonotype_id"] in done
    n_pending = sum(1 for j in jobs if not j["done"])
    # The per-job routing tags go in the TEXT too, not only structuredContent: the strategist
    # routes each group on mhc_class/has_tcr/species/output_needed, and if it only sees the
    # summary count it cannot route and (correctly) refuses to guess. The bulky construct_fasta
    # stays out of the text; the executor reads it back from disk by clonotype_id.
    lines = [
        f"- {j['clonotype_id']} (group {j.get('group_id')}): mhc_class={j.get('mhc_class')}, "
        f"has_tcr={j.get('has_tcr')}, species={j.get('species')}, "
        f"output_needed={j.get('output_needed')}, tool={j.get('tool')}, done={j['done']}"
        for j in jobs
    ]
    head = f"{len(jobs)} jobs ({n_pending} pending, {len(jobs) - n_pending} done)"
    r = _txt(head + (":\n" + "\n".join(lines) if lines else ""))
    r["structuredContent"] = {"jobs": jobs, "pending": n_pending}
    return r


@tool("assign_group_tool",
      "Persist the tool the strategist chose for a group onto its fold jobs, so the executor "
      "filters on the tag by fact and the report shows the right tool per group.",
      {"run_dir": str, "group_id": str, "tool": str})
async def assign_group_tool(args):
    rs = RunState(args["run_dir"])
    jobs = rs.read_stage("foldjobs") if rs.stage_done("foldjobs") else []
    n = 0
    for j in jobs:
        if j.get("group_id") == args["group_id"]:
            j["tool"] = args["tool"]
            n += 1
    rs.write_stage("foldjobs", jobs)
    r = _txt(f"assigned {args['tool']} to {n} jobs in group {args['group_id']}")
    r["structuredContent"] = {"assigned": n, "tool": args["tool"], "group_id": args["group_id"]}
    return r


@tool("list_structure_tools", "List the structure tools and their validity domains for the strategist.",
      {"run_dir": str})
async def list_structure_tools(args):
    # The actionable payload goes in the TEXT too, not only structuredContent: a nested
    # sub-agent may only see the text, and routing on the tool NAME alone (a prior) instead
    # of its validity domain is exactly the failure this avoids.
    lines = []
    for t in structure_tools.as_dicts():
        v = t["validity"]
        dflt = " [default]" if t["is_default"] else ""
        lines.append(
            f"- {t['name']}{dflt}: mhc_class={v['mhc_class']}, needs_tcr={v['needs_tcr']}, "
            f"species={v['species']}, output={t['output_type']}, qc_metric={t['qc_metric']}")
    r = _txt("structure tool registry (route each group to ONE tool by its validity "
             "domain, not its name):\n" + "\n".join(lines))
    r["structuredContent"] = {"tools": structure_tools.as_dicts()}
    return r


@tool("list_compute_routes",
      "List the compute routes (Colab, local GPU, SSH, server), the fields each needs, and "
      "whether its runner is wired, so the intake agent asks the right questions.",
      {"run_dir": str})
async def list_compute_routes(args):
    # Put the fields in the TEXT, not only structuredContent: the intake agent asks the user
    # for exactly a route's required_fields, so if it only sees a one-line stub it invents
    # fields (e.g. a Colab notebook_url that does not exist) and deadlocks. An empty
    # required_fields list must read as "no extra fields", never as "unknown".
    lines = []
    for x in compute_routes.as_dicts():
        req = ", ".join(x["required_fields"]) if x["required_fields"] else "none"
        sec = ", ".join(x["secret_fields"]) if x["secret_fields"] else "none"
        dflt = " [default]" if x["is_default"] else ""
        lines.append(
            f"- {x['name']}{dflt}: required_fields=[{req}], secret_fields=[{sec}], "
            f"artifact={x['artifact_kind']}, runner_wired={x['wired']}")
    r = _txt("compute route registry (collect EVERY required_field of the chosen route; "
             "required_fields=[none] means the route needs no extra fields, so record it "
             "as soon as it is chosen):\n" + "\n".join(lines))
    r["structuredContent"] = {"routes": compute_routes.as_dicts()}
    return r


def _fold_inputs(tool: str, job: dict, clonotype_id: str, clon=None, ann=None) -> dict:
    """Shape a fold job's construct into the tool's Colab inputs. mhcfine, affinetune and
    tcrdock each take a cognate + scramble pair (keys prefixed by clonotype id so per-clonotype
    output files never collide): mhcfine gets MHC heavy + peptide, affinetune gets
    {mhc, b2m, peptide}, tcrdock gets the gene-level TSV row. tcrdock is built from the
    Clonotype+Annotation (gene names + CDR3 + peptide + HLA), which the FoldJob does NOT carry
    (it keeps only construct_fasta), so build_fold_notebook reads them back from the ingest and
    annotate stages and passes them here. Unwired tools embed the raw construct for their
    fail-loud scaffold."""
    from .tools import mhcfine_inputs, affinetune_inputs, tcrdock_inputs, protenix_inputs
    if tool == "tcrdock":
        if clon is None or ann is None:
            raise ValueError(
                f"tcrdock needs the clonotype and annotation (gene-level fields) for "
                f"{clonotype_id}; not found in the ingest/annotate stages")
        built = tcrdock_inputs.build(clon, ann)
        return {f"{clonotype_id}_{k}": v for k, v in built.items()}
    fasta = job["construct_fasta"]
    if tool == "protenix":
        built = protenix_inputs.build(fasta)
        return {f"{clonotype_id}_{k}": v for k, v in built.items()}
    if tool == "mhcfine":
        built = mhcfine_inputs.build(fasta)
        return {f"{clonotype_id}_{k}": v for k, v in built.items()}
    if tool == "affinetune":
        built = affinetune_inputs.build(fasta)
        return {f"{clonotype_id}_{k}": v for k, v in built.items()}
    return {"construct_fasta": fasta}


def _tcrdock_ctx(run_dir, cid, tool):
    """tcrdock needs the gene-level Clonotype+Annotation (the FoldJob keeps only construct_fasta),
    so read them back from the persisted stages. (None, None) for every other tool."""
    if tool != "tcrdock":
        return None, None
    clon = next((c for c in _load(run_dir, "ingest", Clonotype) if c.id == cid), None)
    ann = next((a for a in _load(run_dir, "annotate", Annotation) if a.clonotype_id == cid), None)
    return clon, ann


def _write_artifact(run_dir, stem, tool, inputs, route):
    """Write ONE fold artifact (Colab notebook or bash script, chosen by the route) for the given
    inputs and return (path, kind, wired). Shared by the per-clonotype and per-group tools; the
    notebook/script builders already loop over the inputs dict, so merged multi-clonotype inputs
    fold in a single run."""
    import json as _json
    from .tools.notebook import build_notebook
    from .tools.protenix_script import build as build_script
    kind = compute_routes.artifact_kind_for(route)
    wired = compute_routes.is_wired(route)
    if kind == "colab_notebook":
        from .tools.notebook import is_tool_wired
        nb_dir = Path(run_dir) / "notebooks"
        nb_dir.mkdir(parents=True, exist_ok=True)
        out = nb_dir / f"{stem}.ipynb"
        out.write_text(_json.dumps(build_notebook(tool, inputs), indent=1))
        # A gated tool (af3) yields a fail-loud stub notebook, not a runnable recipe; report it
        # as not wired so the executor hands it over honestly instead of claiming it will run.
        wired = wired and is_tool_wired(tool)
    else:  # bash_script (local_gpu, and the honest ssh/server handoff)
        sc_dir = Path(run_dir) / "scripts"
        sc_dir.mkdir(parents=True, exist_ok=True)
        out = sc_dir / f"{stem}.sh"
        if tool != "protenix":
            # Only the Protenix bash runner is wired. Never emit a Protenix script fed another
            # tool's inputs (it would fold the wrong model on wrong-shaped JSON and claim success);
            # write a fail-loud stub and report the route as not wired for this tool.
            out.write_text(f"#!/usr/bin/env bash\necho 'no bash runner wired for {tool}; use the "
                           f"colab route for {tool} folds' >&2\nexit 1\n")
            return out, kind, False
        spec = intake.load_intake(run_dir)
        working = (spec.route_params.get("working_path") if spec else None) or "."
        out.write_text(build_script(inputs, working_path=working))
    return out, kind, wired


@tool("build_fold_notebook",
      "Build the Colab notebook for one clonotype's fold job, write it under the run dir, and return its path.",
      {"run_dir": str, "clonotype_id": str, "tool": str})
async def build_fold_notebook(args):
    import json as _json
    from .tools.notebook import build_notebook
    rs = RunState(args["run_dir"])
    jobs = rs.read_stage("foldjobs") if rs.stage_done("foldjobs") else []
    cid = args["clonotype_id"]
    job = next((j for j in jobs if j["clonotype_id"] == cid), None)
    if job is None:
        return _txt(f"no fold job for {cid}")
    tool = args["tool"]
    # tcrdock needs the gene-level Clonotype+Annotation, which the FoldJob does not carry;
    # read them back from the persisted ingest/annotate stages (tcr_explorer already enriched
    # them upstream, so nothing is re-derived here).
    clon = ann = None
    if tool == "tcrdock":
        clon = next((c for c in _load(args["run_dir"], "ingest", Clonotype) if c.id == cid), None)
        ann = next((a for a in _load(args["run_dir"], "annotate", Annotation)
                    if a.clonotype_id == cid), None)
    nb = build_notebook(tool, _fold_inputs(tool, job, cid, clon, ann))
    nb_dir = Path(args["run_dir"]) / "notebooks"
    nb_dir.mkdir(parents=True, exist_ok=True)
    out = nb_dir / f"{cid}_{tool}.ipynb"
    out.write_text(_json.dumps(nb, indent=1))
    r = _txt(f"notebook for {cid} ({tool}) written to {out}")
    r["structuredContent"] = {"notebook_path": str(out), "clonotype_id": cid, "tool": tool}
    return r


@tool("build_fold_artifact",
      "Build the fold artifact (Colab notebook or local bash script) for one clonotype, "
      "chosen by the compute route, write it under the run dir, and return its path.",
      {"run_dir": str, "clonotype_id": str, "tool": str, "compute_route": str})
async def build_fold_artifact(args):
    rs = RunState(args["run_dir"])
    jobs = rs.read_stage("foldjobs") if rs.stage_done("foldjobs") else []
    cid = args["clonotype_id"]
    job = next((j for j in jobs if j["clonotype_id"] == cid), None)
    if job is None:
        return _txt(f"no fold job for {cid}")
    tool, route = args["tool"], args["compute_route"]
    clon, ann = _tcrdock_ctx(args["run_dir"], cid, tool)
    inputs = _fold_inputs(tool, job, cid, clon, ann)
    out, kind, wired = _write_artifact(args["run_dir"], f"{cid}_{tool}", tool, inputs, route)
    note = "" if wired else f" (route '{route}' runner not wired; run the script yourself)"
    r = _txt(f"{kind} for {cid} ({tool}) via {route} written to {out}{note}")
    r["structuredContent"] = {"artifact_path": str(out), "artifact_kind": kind,
                              "route_wired": wired, "clonotype_id": cid, "tool": tool}
    return r


# Cap clonotypes per artifact so one Colab session's wall-clock/memory stays sane; a bigger
# group is sharded into {gid}_{tool}_partK artifacts a researcher runs independently.
_MAX_BATCH = 16


@tool("build_group_artifact",
      "Build fold artifacts (Colab notebook or bash script) that fold a WHOLE group's pending "
      "clonotypes, batched into as few artifacts as possible (sharded at 16 per artifact), chosen "
      "by the compute route. Use this instead of one artifact per clonotype so a researcher runs "
      "one artifact per group (or a few shards for a big group), not one per TCR. The tool is "
      "taken from the strategist's persisted choice; a mismatching passed tool fails loud.",
      {"run_dir": str, "group_id": str, "tool": str, "compute_route": str})
async def build_group_artifact(args):
    rs = RunState(args["run_dir"])
    jobs = rs.read_stage("foldjobs") if rs.stage_done("foldjobs") else []
    done = rs.read_stage("folds") if rs.stage_done("folds") else {}
    gid, route = args["group_id"], args["compute_route"]
    pend = [j for j in jobs if j.get("group_id") == gid and j["clonotype_id"] not in done]
    if not pend:
        return _txt(f"no pending fold jobs in group {gid}")
    # The strategist's persisted choice (assign_group_tool) is the source of truth; the passed
    # tool is only validated against it, so a prompt slip can never fold a group with a different
    # model than the report will show.
    persisted = {j.get("tool") for j in pend if j.get("tool")}
    if len(persisted) > 1:
        return _txt(f"group {gid} has mixed persisted tools {sorted(persisted)}; "
                    f"assign_group_tool one tool per group first")
    tool = next(iter(persisted)) if persisted else args["tool"]
    if persisted and args.get("tool") and args["tool"] != tool:
        return _txt(f"tool mismatch for group {gid}: strategist persisted '{tool}', executor "
                    f"passed '{args['tool']}'; not building")
    # Only Protenix has a bash runner. For a bash route + another tool, short-circuit to one
    # fail-loud stub before building any per-clonotype inputs (which we could not run anyway).
    if compute_routes.artifact_kind_for(route) == "bash_script" and tool != "protenix":
        out, kind, _ = _write_artifact(args["run_dir"], f"{gid}_{tool}", tool, {}, route)
        r = _txt(f"no bash runner wired for {tool} in group {gid}; wrote fail-loud stub {out}")
        r["structuredContent"] = {
            "artifacts": [{"artifact_path": str(out),
                           "clonotypes": [j["clonotype_id"] for j in pend]}],
            "artifact_kind": kind, "route_wired": False, "tool": tool, "group_id": gid,
            "n_clonotypes": len(pend), "n_artifacts": 1}
        return r
    shards = [pend[i:i + _MAX_BATCH] for i in range(0, len(pend), _MAX_BATCH)]
    artifacts, kind, wired = [], None, True
    for k, chunk in enumerate(shards):
        merged = {}
        for j in chunk:
            cid = j["clonotype_id"]
            clon, ann = _tcrdock_ctx(args["run_dir"], cid, tool)
            merged.update(_fold_inputs(tool, j, cid, clon, ann))
        stem = f"{gid}_{tool}" if len(shards) == 1 else f"{gid}_{tool}_part{k + 1}"
        out, kind, wired = _write_artifact(args["run_dir"], stem, tool, merged, route)
        artifacts.append({"artifact_path": str(out),
                          "clonotypes": [j["clonotype_id"] for j in chunk]})
    note = "" if wired else f" (route '{route}' runner not wired for {tool}; run it yourself)"
    r = _txt(f"{kind} for group {gid}: {len(pend)} clonotypes in {len(shards)} artifact(s) "
             f"({tool}) via {route}{note}")
    r["structuredContent"] = {"artifacts": artifacts, "artifact_kind": kind, "route_wired": wired,
                              "tool": tool, "group_id": gid, "n_clonotypes": len(pend),
                              "n_artifacts": len(shards)}
    return r


@tool("record_fold_result", "Record the model paths a fold produced for one clonotype, with the tool used.",
      {"run_dir": str, "clonotype_id": str, "model_paths": list, "tool": str})
async def record_fold_result(args):
    rs = RunState(args["run_dir"])
    done = rs.read_stage("folds") if rs.stage_done("folds") else {}
    done[args["clonotype_id"]] = {"paths": args["model_paths"],
                                  "tool": args.get("tool", "protenix")}
    rs.write_stage("folds", done)
    return _txt(f"recorded {len(args['model_paths'])} models for {args['clonotype_id']} via {args.get('tool', 'protenix')}")


@tool("record_local_folds",
      "Scan <run_dir>/out for fold CIFs already on disk (local_gpu run, or a Colab download "
      "unzipped there) and record them per clonotype so QC can proceed.",
      {"run_dir": str, "tool": str})
async def record_local_folds(args):
    rs = RunState(args["run_dir"])
    done = rs.read_stage("folds") if rs.stage_done("folds") else {}
    found = scan_recorded_folds(args["run_dir"], args.get("tool", "protenix"))
    added = {cid: rec for cid, rec in found.items() if cid not in done}
    done.update(added)
    rs.write_stage("folds", done)
    r = _txt(f"recorded {len(added)} clonotypes from disk: {sorted(added)}")
    r["structuredContent"] = {"recorded": len(added), "clonotypes": sorted(added)}
    return r


@tool("qc_structure", "Score a fold (per-group threshold) and return a skeptical verdict; output-type aware.",
      {"run_dir": str, "clonotype_id": str, "scramble_threshold": float,
       "output_type": str, "tool": str})
async def qc_structure(args):
    from .qc import verdict_binding, load_chains, common_checks, score_pose, verdict_groove
    rs = RunState(args["run_dir"])
    done = rs.read_stage("folds") if rs.stage_done("folds") else {}
    rec = done.get(args["clonotype_id"], {})
    if isinstance(rec, list):                      # back-compat with the old list shape
        rec = {"paths": rec, "tool": "protenix"}
    paths = rec.get("paths", [])
    tool = args.get("tool", rec.get("tool", "protenix"))
    metric = structure_tools.qc_metric_for(tool)   # tool decides, not the agent
    if not paths:
        res = QCResult(args["clonotype_id"], "qc_failed", "no model recorded", tool=tool)
        validity_summary = "no model"
    elif metric == "binding_score":
        try:
            score = float(Path(paths[0]).read_text().strip())
        except (ValueError, OSError):
            # A corrupt/empty .score (a truncated download) must fail this ONE clonotype, not
            # crash the QC tool mid-run and abort every remaining clonotype.
            res = QCResult(args["clonotype_id"], "qc_failed",
                           f"unreadable score file {Path(paths[0]).name}", tool=tool)
            validity_summary = "no score"
        else:
            # Set the threshold from this group's OWN scramble null (the sibling scramble score the
            # fold wrote), never a global number; fall back to the caller's explicit threshold only
            # if that scramble file is absent. tcrdock's validated flu M1 null: cognate score
            # -11.219 beats scramble -20.574 (interface PAE 11.2 vs 20.6).
            threshold = _scramble_null(paths[0])
            if threshold is None:
                threshold = args["scramble_threshold"]
            res = verdict_binding(score, threshold, args["clonotype_id"], tool=tool)
            validity_summary = "n/a (binding score)"
    else:
        # cdr3_peptide needs the full TCR-pMHC (A-E); peptide_groove is an mhcfine
        # pose = MHC heavy (C) + peptide (E) only, with no b2m/TCR modelled.
        expected = {"A", "B", "C", "D", "E"} if metric == "cdr3_peptide" else {"C", "E"}
        chains = load_chains(paths[0])
        cc = common_checks(chains, expected)
        validity_summary = "valid" if cc["ok"] else "; ".join(cc["issues"])
        if not cc["ok"]:
            failed = "pose_failed" if metric == "peptide_groove" else "qc_failed"
            res = QCResult(args["clonotype_id"], failed,
                           "; ".join(cc["issues"]), tool=tool)
        elif metric == "peptide_groove":
            res = verdict_groove(score_pose(chains), args["clonotype_id"], tool=tool)
        else:  # cdr3_peptide
            from .qc import ensemble_contact
            # Protenix emits several samples whose docking pose varies a lot, and it writes a
            # cognate + scramble pair per clonotype ({cid}_cognate / {cid}_scramble in the path).
            # Ensemble the CDR3-peptide contact over the cognate samples and calibrate it against
            # this clonotype's OWN scramble ensemble (the per-clonotype null), never a global
            # number. Fall back to the caller's explicit threshold + a single model only when the
            # recorded paths carry no cognate/scramble split (legacy single-model records).
            # Protenix marks the construct in the DIRECTORY ({cid}_cognate/...), not the
            # filename (which is cognate_sample_N.cif), so match on the full path, not the
            # basename. Learned from the first real repatriated fold.
            cognate = [p for p in paths if "_cognate" in str(p)]
            scramble = [p for p in paths if "_scramble" in str(p)]
            if cognate:
                cog, _, n_valid = ensemble_contact(cognate)
                scr, _, _ = ensemble_contact(scramble) if scramble else (None, 0, 0)
                thr = scr if scr is not None else args["scramble_threshold"]
                if cog is None:  # no cognate sample parsed to a full 5-chain complex
                    res = QCResult(args["clonotype_id"], "qc_failed",
                                   "no cognate sample parsed to a full 5-chain TCR-pMHC")
                else:
                    s = {"cdr3_pep_atoms": cog, "n_chains": 5,
                         "clonotype_id": args["clonotype_id"]}
                    res = verdict(s, thr)
            else:
                s = score_model(paths[0])
                s["clonotype_id"] = args["clonotype_id"]
                res = verdict(s, args["scramble_threshold"])
            res.tool = tool
            res.calibration_basis = "scramble_null"
    # Honesty guard: a clonotype folded on a poly-G stub V-domain (germline reconstruction
    # failed) is not a real TCR, so it can never earn a positive verdict no matter how the
    # fold scored. Downgrade reliable/presented to suspect and say why.
    job = next((j for j in (rs.read_stage("foldjobs") if rs.stage_done("foldjobs") else [])
                if j["clonotype_id"] == args["clonotype_id"]), None)
    if job is not None and not job.get("tcr_reconstructed", True) and \
            res.qc_verdict in ("reliable", "presented"):
        res.qc_verdict = "suspect"
        res.reason = "V-domain was a poly-G stub (germline reconstruction failed), not a real TCR"
    qcs = rs.read_stage("qc") if rs.stage_done("qc") else []
    qcs = [q for q in qcs if q["clonotype_id"] != res.clonotype_id]
    qcs.append(asdict(res))
    rs.write_stage("qc", qcs)
    validity = rs.read_stage("validity") if rs.stage_done("validity") else {}
    validity[args["clonotype_id"]] = validity_summary
    rs.write_stage("validity", validity)
    r = _txt(f"{res.clonotype_id}: {res.qc_verdict} ({res.reason})")
    r["structuredContent"] = {"qc_verdict": res.qc_verdict, "reason": res.reason}
    return r


@tool("record_intake",
      "Persist the intake brief (data type, input path, question, compute route, route params) "
      "to run_dir/intake.json so the run can proceed and later resume. Secrets are stripped.",
      {"run_dir": str, "data_type": str, "input_path": str, "question": str,
       "compute_route": str, "route_params": dict})
async def record_intake(args):
    from .intake import IntakeSpec
    spec = IntakeSpec(args["data_type"], args["input_path"], args["question"],
                      args["compute_route"], args.get("route_params", {}))
    path = intake.save_intake(args["run_dir"], spec)
    r = _txt(f"intake recorded to {path}")
    r["structuredContent"] = {"intake_path": path, "compute_route": spec.compute_route}
    return r


@tool("render_final_report", "Render the self contained HTML report for the run.",
      {"run_dir": str})
async def render_final_report(args):
    import json as _json
    from .report import msa_basis_from_manifest
    clons = _load(args["run_dir"], "ingest", Clonotype)
    anns = _load(args["run_dir"], "annotate", Annotation)
    rs = RunState(args["run_dir"])
    qcs = [QCResult(**d) for d in (rs.read_stage("qc") if rs.stage_done("qc") else [])]
    fjs = rs.read_stage("foldjobs") if rs.stage_done("foldjobs") else []
    manifests = {}
    for mp in Path(args["run_dir"]).rglob("*_msa_manifest.json"):
        cid = mp.name[:-len("_msa_manifest.json")]
        try:
            manifests[cid] = _json.loads(mp.read_text())
        except (ValueError, OSError):
            pass
    msa_basis = {
        j["clonotype_id"]: (msa_basis_from_manifest(manifests[j["clonotype_id"]])
                            if j["clonotype_id"] in manifests else j.get("msa_basis"))
        for j in fjs
    }
    validity = rs.read_stage("validity") if rs.stage_done("validity") else {}
    # A clonotype folded on a poly-G stub V-domain is flagged so the report never presents a
    # fabricated framework as a real TCR (QC already caps its verdict).
    tcr_stub = {j["clonotype_id"]: not j.get("tcr_reconstructed", True) for j in fjs}
    # Scale guard: show only the selected/folded clonotypes, not one row per ingested clonotype,
    # so the HTML does not blow up to 20k rows on a large repertoire. Fall back to all when
    # nothing has been selected yet (small/early runs).
    selected = {j["clonotype_id"] for j in fjs} | {q.clonotype_id for q in qcs}
    shown = [c for c in clons if c.id in selected] or clons
    html = render_report(shown, anns, qcs, msa_basis=msa_basis, validity=validity,
                         tcr_stub=tcr_stub)
    out = Path(args["run_dir"]) / "report.html"
    out.write_text(html)
    r = _txt(f"report written to {out}")
    r["structuredContent"] = {"report_path": str(out)}
    return r


def build_server():
    return create_sdk_mcp_server(name="rep2struct", version="0.1.0", tools=[
        ingest_repertoire, annotate_specificity, prep_and_select, list_fold_jobs,
        assign_group_tool, list_structure_tools, list_compute_routes, build_fold_notebook,
        build_fold_artifact, build_group_artifact, record_fold_result,
        record_local_folds, qc_structure, render_final_report, record_intake,
    ])
