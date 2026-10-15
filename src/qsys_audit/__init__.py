"""qsys_audit: offline cross-layer auditing of concurrent crosstalk in shared
quantum processors.

The scientific classes are imported lazily so that ``import qsys_audit`` and
``qsys_audit.__version__`` work before the heavy dependencies (QuTiP, Qiskit,
MQT Bench) are installed.
"""

from __future__ import annotations

from importlib import import_module
from typing import TYPE_CHECKING, Any

from .paths import default_config, repo_root

__version__ = "1.0.0"

_LAZY: dict[str, str] = {
    "PhysicsEngine": ".physics_engine",
    "StressProfile": ".physics_engine",
    "ReconstructionEngine": ".reconstruction",
    "WorkloadReplay": ".workload_replay",
}


def __getattr__(name: str) -> Any:
    target = _LAZY.get(name)
    if target is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(import_module(target, __name__), name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(_LAZY))


if TYPE_CHECKING:  # pragma: no cover
    from .physics_engine import PhysicsEngine, StressProfile
    from .reconstruction import ReconstructionEngine
    from .workload_replay import WorkloadReplay

__all__ = [
    "PhysicsEngine",
    "StressProfile",
    "ReconstructionEngine",
    "WorkloadReplay",
    "default_config",
    "repo_root",
    "__version__",
]
