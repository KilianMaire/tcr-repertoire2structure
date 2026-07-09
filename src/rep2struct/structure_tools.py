from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class StructureTool:
    name: str
    validity: dict            # {"mhc_class": set[int], "needs_tcr": Optional[bool], "species": str}
    output_type: str          # "structure" | "binding_score"
    strengths: str
    limits: str
    colab_adapter: str
    qc_metric: str = "cdr3_peptide"     # cdr3_peptide | peptide_groove | binding_score
    is_default: bool = False


REGISTRY: list[StructureTool] = [
    StructureTool(
        name="protenix",
        validity={"mhc_class": {1, 2}, "needs_tcr": None, "species": "any"},
        output_type="structure",
        strengths="full 3-chain TCR-pMHC fold, default workhorse",
        limits="imposes canonical geometry even on non-binders (basis of skeptical QC)",
        colab_adapter="protenix_colab",
        qc_metric="cdr3_peptide",
        is_default=True,
    ),
    StructureTool(
        name="af3",
        validity={"mhc_class": {1, 2}, "needs_tcr": None, "species": "any"},
        output_type="structure",
        strengths="AlphaFold3-class accuracy when weights are available",
        limits="gated model weights; only if the user has them",
        colab_adapter="af3_colab",
        qc_metric="cdr3_peptide",
    ),
    StructureTool(
        name="mhcfine",
        validity={"mhc_class": {1}, "needs_tcr": False, "species": "any"},
        output_type="structure",
        strengths="most precise class I peptide pose (live-validated 0.54A peptide "
                  "backbone RMSD vs the 6VRN crystal, matching the paper's ~0.66A)",
        limits="class I only, no TCR. The pose does NOT discriminate binders from "
               "non-binders: live calibration showed a scrambled peptide seats in the "
               "groove just as deeply and neither plddt nor groove-contact separates "
               "cognate from scramble, so a mhcfine pose is placement only, never "
               "evidence of recognition",
        colab_adapter="mhcfine_colab",
        qc_metric="peptide_groove",
    ),
    StructureTool(
        name="tcrdock",
        validity={"mhc_class": {1, 2}, "needs_tcr": True, "species": "any"},
        output_type="structure",
        strengths="TCR:pMHC interface and V-domain anchoring",
        limits="template-coverage limited; class II not systematically benchmarked. "
               "Live-validated class I 2026-07-09: folds a docking structure (a single "
               "merged chain, not A-E) but recognition is judged by the peptide<->TCR "
               "interface PAE, which discriminated a cognate flu M1 TCR (pae ~11, peptide "
               "pLDDT 86) from its scramble (pae ~21, pLDDT 65). QC threshold is set per group "
               "from its OWN folded scramble null (score = -interface_pae): validated flu M1 "
               "null cognate -11.219 beats scramble -20.574",
        colab_adapter="tcrdock_colab",
        # output_type stays "structure" (tcrdock IS the TCR:pMHC docking tool, must remain
        # selectable for structure jobs), but the single-chain output cannot feed the
        # cdr3_peptide chain-contact path, so QC scores the interface PAE via verdict_binding.
        qc_metric="binding_score",
    ),
    StructureTool(
        name="affinetune",
        validity={"mhc_class": {1, 2}, "needs_tcr": False, "species": "any"},
        output_type="binding_score",
        strengths="is-this-peptide-presented classifier, class I and II",
        limits="returns a presentation score, not a structure",
        colab_adapter="affinetune_colab",
        qc_metric="binding_score",
    ),
]


def output_type_for(name: str) -> str:
    """Return a tool's output_type from the registry, defaulting to
    "structure" for an unknown tool name."""
    for t in REGISTRY:
        if t.name == name:
            return t.output_type
    return "structure"


def qc_metric_for(name: str) -> str:
    """Return a tool's qc_metric from the registry, defaulting to
    "cdr3_peptide" for an unknown tool name."""
    for t in REGISTRY:
        if t.name == name:
            return t.qc_metric
    return "cdr3_peptide"


def get_default() -> StructureTool:
    return next(t for t in REGISTRY if t.is_default)


def _species_ok(tool: StructureTool, species: str) -> bool:
    return tool.validity["species"] == "any" or tool.validity["species"] == species


def _tcr_ok(tool: StructureTool, has_tcr: bool) -> bool:
    need = tool.validity["needs_tcr"]
    return need is None or need == has_tcr


def tools_for(mhc_class: int, has_tcr: bool, species: str, output_needed: str) -> list[StructureTool]:
    return [
        t for t in REGISTRY
        if t.output_type == output_needed
        and mhc_class in t.validity["mhc_class"]
        and _tcr_ok(t, has_tcr)
        and _species_ok(t, species)
    ]


def is_covered(mhc_class: int, has_tcr: bool, species: str, output_needed: str) -> bool:
    return len(tools_for(mhc_class, has_tcr, species, output_needed)) > 0


def as_dicts() -> list[dict]:
    out = []
    for t in REGISTRY:
        v = dict(t.validity)
        v["mhc_class"] = sorted(v["mhc_class"])
        out.append({
            "name": t.name, "validity": v, "output_type": t.output_type,
            "strengths": t.strengths, "limits": t.limits, "qc_metric": t.qc_metric, "is_default": t.is_default,
        })
    return out
