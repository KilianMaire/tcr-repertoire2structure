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


def _load(rd, name, cls):
    rs = RunState(rd)
    return [cls(**d) for d in rs.read_stage(name)] if rs.stage_done(name) else []


@tool("ingest_repertoire", "Parse a 10x contig CSV into paired clonotypes and persist them.",
      {"run_dir": str, "csv_path": str})
async def ingest_repertoire(args):
    clons = standardize_alleles(parse_10x(args["csv_path"]), assign_fn=_CFG["assign_fn"])
    RunState(args["run_dir"]).write_stage("ingest", clons)
    r = _txt(f"{len(clons)} clonotypes ingested")
    r["structuredContent"] = {"clonotypes": len(clons)}
    return r


@tool("annotate_specificity", "Annotate persisted clonotypes with candidate epitopes by TCRdist. Never forces a label.",
      {"run_dir": str})
async def annotate_specificity(args):
    clons = _load(args["run_dir"], "ingest", Clonotype)
    anns = annotate(clons, sim_fn=_CFG["sim_fn"])
    RunState(args["run_dir"]).write_stage("annotate", anns)
    tiers = {}
    for a in anns:
        tiers[a.confidence_tier] = tiers.get(a.confidence_tier, 0) + 1
    r = _txt(f"annotated {len(anns)} clonotypes: {tiers}")
    r["structuredContent"] = {"tiers": tiers}
    return r


@tool("prep_and_select", "Rank clonotypes and build class I fold constructs for the top N.",
      {"run_dir": str, "top_n": int})
async def prep_and_select(args):
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
    jobs = [build_construct(c, a, seqs, mhc) for c, a in foldable if a.hla in mhc]
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
    r = _txt(f"{len(jobs)} jobs")
    r["structuredContent"] = {"jobs": jobs}
    return r


@tool("list_structure_tools", "List the structure tools and their validity domains for the strategist.",
      {"run_dir": str})
async def list_structure_tools(args):
    r = _txt("structure tool registry")
    r["structuredContent"] = {"tools": structure_tools.as_dicts()}
    return r


def _fold_inputs(tool: str, job: dict, clonotype_id: str) -> dict:
    """Shape a fold job's construct into the tool's Colab inputs. mhcfine takes the MHC
    heavy chain + peptide (keys prefixed by clonotype id so per-clonotype output files
    never collide). Unwired tools embed the raw construct for their fail-loud scaffold."""
    from .tools import mhcfine_inputs
    fasta = job["construct_fasta"]
    if tool == "mhcfine":
        built = mhcfine_inputs.build(fasta)
        return {f"{clonotype_id}_{k}": v for k, v in built.items()}
    return {"construct_fasta": fasta}


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
    nb = build_notebook(tool, _fold_inputs(tool, job, cid))
    nb_dir = Path(args["run_dir"]) / "notebooks"
    nb_dir.mkdir(parents=True, exist_ok=True)
    out = nb_dir / f"{cid}_{tool}.ipynb"
    out.write_text(_json.dumps(nb, indent=1))
    r = _txt(f"notebook for {cid} ({tool}) written to {out}")
    r["structuredContent"] = {"notebook_path": str(out), "clonotype_id": cid, "tool": tool}
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
        score = float(Path(paths[0]).read_text().strip())
        res = verdict_binding(score, args["scramble_threshold"], args["clonotype_id"], tool=tool)
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
            s = score_model(paths[0])
            s["clonotype_id"] = args["clonotype_id"]
            res = verdict(s, args["scramble_threshold"])
            res.tool = tool
            res.calibration_basis = "scramble_null"
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


@tool("render_final_report", "Render the self contained HTML report for the run.",
      {"run_dir": str})
async def render_final_report(args):
    clons = _load(args["run_dir"], "ingest", Clonotype)
    anns = _load(args["run_dir"], "annotate", Annotation)
    rs = RunState(args["run_dir"])
    qcs = [QCResult(**d) for d in (rs.read_stage("qc") if rs.stage_done("qc") else [])]
    fjs = rs.read_stage("foldjobs") if rs.stage_done("foldjobs") else []
    msa_basis = {j["clonotype_id"]: j.get("msa_basis") for j in fjs}
    validity = rs.read_stage("validity") if rs.stage_done("validity") else {}
    html = render_report(clons, anns, qcs, msa_basis=msa_basis, validity=validity)
    out = Path(args["run_dir"]) / "report.html"
    out.write_text(html)
    r = _txt(f"report written to {out}")
    r["structuredContent"] = {"report_path": str(out)}
    return r


def build_server():
    return create_sdk_mcp_server(name="rep2struct", version="0.1.0", tools=[
        ingest_repertoire, annotate_specificity, prep_and_select, list_fold_jobs,
        list_structure_tools, build_fold_notebook, record_fold_result, qc_structure,
        render_final_report,
    ])
