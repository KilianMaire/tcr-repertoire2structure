from __future__ import annotations
from dataclasses import dataclass, field, asdict
from typing import Optional

@dataclass
class Clonotype:
    id: str                      # content hash of the defining tuple
    trav: str
    cdr3a: str
    trbv: str
    cdr3b: str
    size: int                    # cell count
    trav_allele: Optional[str] = None
    trbv_allele: Optional[str] = None
    traj: Optional[str] = None   # J genes, needed to reconstruct the V domain
    trbj: Optional[str] = None

@dataclass
class Annotation:
    clonotype_id: str
    annotatable: bool
    confidence_tier: str         # high, medium, low, unannotatable
    tcrdist: Optional[float] = None
    epitope: Optional[str] = None
    hla: Optional[str] = None
    antigen: Optional[str] = None
    neighbour_id: Optional[str] = None

@dataclass
class FoldJob:
    clonotype_id: str
    construct_fasta: str         # A..E chains
    msa_ref: Optional[str] = None
    status: str = "pending"      # pending, done, failed
    model_paths: list[str] = field(default_factory=list)
    mhc_class: int = 1           # 1 or 2
    has_tcr: bool = True
    species: str = "human"
    output_needed: str = "structure"   # structure | binding_score
    tool: Optional[str] = None   # tool the strategist assigned
    group_id: Optional[str] = None
    msa_basis: Optional[str] = None   # local | colab_cpu | none

@dataclass
class QCResult:
    clonotype_id: str
    qc_verdict: str              # reliable, suspect, qc_failed, presented, not_presented
    reason: str
    dockq: Optional[float] = None
    cdr3_pep_atoms: Optional[float] = None
    crossing_angle: Optional[float] = None
    tool: Optional[str] = None
    calibration_basis: Optional[str] = None

def to_jsonable(obj):
    return asdict(obj) if hasattr(obj, "__dataclass_fields__") else obj
