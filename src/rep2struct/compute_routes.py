from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ComputeRoute:
    name: str
    description: str
    required_fields: tuple[str, ...]   # non secret connection fields the agent must collect
    secret_fields: tuple[str, ...]     # collected in session memory only, never persisted
    artifact_kind: str                 # "colab_notebook" | "bash_script"
    wired: bool                        # True: R2S produces a runnable artifact the user runs
    is_default: bool = False


REGISTRY: list[ComputeRoute] = [
    ComputeRoute(
        name="colab",
        description="Google Colab. The simplest option: R2S writes a self contained "
                    "notebook the user uploads and runs on a Colab GPU.",
        required_fields=(),
        secret_fields=(),
        artifact_kind="colab_notebook",
        wired=True,
        is_default=True,
    ),
    ComputeRoute(
        name="local_gpu",
        description="A GPU machine the user has a shell on (for example a local H100). "
                    "R2S writes a bash script that folds with Protenix and leaves the "
                    "CIFs on that machine.",
        required_fields=("working_path",),
        secret_fields=(),
        artifact_kind="bash_script",
        wired=True,
    ),
    ComputeRoute(
        name="ssh",
        description="A remote host reached over SSH. R2S collects the connection details "
                    "and hands the user a job script plus scp/sbatch instructions. The SSH "
                    "runner itself is not wired yet.",
        required_fields=("host", "user", "remote_path"),
        secret_fields=("password",),
        artifact_kind="bash_script",
        wired=False,
    ),
    ComputeRoute(
        name="server",
        description="A shared server the user names by address. Same handoff as ssh; the "
                    "server runner is not wired yet.",
        required_fields=("address", "path"),
        secret_fields=(),
        artifact_kind="bash_script",
        wired=False,
    ),
]


def get_default() -> ComputeRoute:
    return next(r for r in REGISTRY if r.is_default)


def by_name(name: str) -> ComputeRoute:
    try:
        return next(r for r in REGISTRY if r.name == name)
    except StopIteration:
        raise ValueError(f"unknown compute route {name!r}") from None


def recommend(context: str = "") -> ComputeRoute:
    """The simplest route for a user who does not know: the default (Colab)."""
    return get_default()


def artifact_kind_for(name: str) -> str:
    return by_name(name).artifact_kind


def is_wired(name: str) -> bool:
    return by_name(name).wired


def required_fields_for(name: str) -> tuple[str, ...]:
    return by_name(name).required_fields


def as_dicts() -> list[dict]:
    return [
        {
            "name": r.name,
            "description": r.description,
            "required_fields": list(r.required_fields),
            "secret_fields": list(r.secret_fields),
            "artifact_kind": r.artifact_kind,
            "wired": r.wired,
            "is_default": r.is_default,
        }
        for r in REGISTRY
    ]
