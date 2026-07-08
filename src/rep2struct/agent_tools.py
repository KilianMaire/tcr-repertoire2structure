from __future__ import annotations
from dataclasses import asdict
from pathlib import Path
from claude_agent_sdk import tool, create_sdk_mcp_server
from .runstate import RunState
from .ingest import parse_10x, standardize_alleles
from .annotate import annotate
from .foldprep import select_top, build_construct
from .qc import score_model, verdict
from .report import render_report
from .schema import Clonotype, Annotation, QCResult

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


def _tcr_seq_stub(clonotype):
    return {"A": "G" * 10 + clonotype.cdr3a, "B": "G" * 10 + clonotype.cdr3b}


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
    seqs = {c.id: _tcr_seq_stub(c) for c, _ in foldable}
    mhc = {a.hla: {"heavy": "H" * 20, "b2m": "M" * 20} for _, a in foldable}
    jobs = [build_construct(c, a, seqs, mhc) for c, a in foldable]
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


@tool("record_fold_result", "Record the model paths a fold produced for one clonotype.",
      {"run_dir": str, "clonotype_id": str, "model_paths": list})
async def record_fold_result(args):
    rs = RunState(args["run_dir"])
    done = rs.read_stage("folds") if rs.stage_done("folds") else {}
    done[args["clonotype_id"]] = args["model_paths"]
    rs.write_stage("folds", done)
    return _txt(f"recorded {len(args['model_paths'])} models for {args['clonotype_id']}")


@tool("qc_structure", "Score a predicted structure and return a skeptical reliable or suspect verdict.",
      {"run_dir": str, "clonotype_id": str, "scramble_threshold": float})
async def qc_structure(args):
    rs = RunState(args["run_dir"])
    done = rs.read_stage("folds") if rs.stage_done("folds") else {}
    paths = done.get(args["clonotype_id"], [])
    if not paths:
        res = QCResult(args["clonotype_id"], "qc_failed", "no model recorded")
    else:
        s = score_model(paths[0])
        s["clonotype_id"] = args["clonotype_id"]
        res = verdict(s, args["scramble_threshold"])
    qcs = rs.read_stage("qc") if rs.stage_done("qc") else []
    qcs = [q for q in qcs if q["clonotype_id"] != res.clonotype_id]
    qcs.append(asdict(res))
    rs.write_stage("qc", qcs)
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
    html = render_report(clons, anns, qcs)
    out = Path(args["run_dir"]) / "report.html"
    out.write_text(html)
    r = _txt(f"report written to {out}")
    r["structuredContent"] = {"report_path": str(out)}
    return r


def build_server():
    return create_sdk_mcp_server(name="rep2struct", version="0.1.0", tools=[
        ingest_repertoire, annotate_specificity, prep_and_select, list_fold_jobs,
        record_fold_result, qc_structure, render_final_report,
    ])
