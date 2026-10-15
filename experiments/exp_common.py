"""Shared helpers for the VeriXTalk experiments.

This module ONLY imports from the qsys_audit pipeline (src/qsys_audit). It adds
no new physics: channels come from `perform_qpt_pipeline`, workload fidelity
from `benchmark_workload`, CPTP projection / Kraus / verdict logic from the
existing reconstruction + CLI code.

Path conventions used by every experiment script:
    ROOT    repository root
    HERE    this directory (experiments/)
    OUTDIR  where result CSVs are written (defaults to HERE)
    FIGDIR  where figures are written (defaults to HERE/figs)
Both OUTDIR and FIGDIR can be redirected with the environment variables
QVSH_OUTDIR and QVSH_FIGDIR.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np
import qutip as qt
import yaml

GLOBAL_SEED = 1234

ROOT = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent
OUTDIR = Path(os.environ.get("QVSH_OUTDIR", HERE))
FIGDIR = Path(os.environ.get("QVSH_FIGDIR", HERE / "figs"))
OUTDIR.mkdir(parents=True, exist_ok=True)
FIGDIR.mkdir(parents=True, exist_ok=True)

SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from qsys_audit.physics_engine import PhysicsEngine, StressProfile  # noqa: E402
from qsys_audit.reconstruction import ReconstructionEngine  # noqa: E402
from qsys_audit.workload_replay import WorkloadReplay  # noqa: E402


def load_config() -> Dict:
    """Load config.yaml from the repository root, or from the installed package."""
    cfg_path = ROOT / "config.yaml"
    if not cfg_path.is_file():
        from qsys_audit.paths import default_config

        cfg_path = default_config()
    with open(cfg_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def build_engines():
    cfg = load_config()
    physics = PhysicsEngine(cfg)
    recon = ReconstructionEngine(physics, cfg)
    replay = WorkloadReplay(physics, cfg)
    return cfg, physics, recon, replay


def classify_verdict(d_val: float, dlog: float, dinc: float, cfg: Dict) -> str:
    """Identical to audit_cli.classify_verdict (thresholds from config.yaml)."""
    t = cfg.get("verification_thresholds", {})
    if (
        d_val <= t.get("eps_d_pass", 0.02)
        and dlog <= t.get("eps_log_pass", 0.05)
        and dinc <= t.get("eps_inc_pass", 0.05)
    ):
        return "PASS"
    if (
        d_val >= t.get("d_hard", 0.30)
        or dlog >= t.get("log_hard", 0.90)
        or dinc >= t.get("inc_hard", 0.50)
    ):
        return "HARD FAIL"
    return "SOFT FAIL"


def channel_D_and_verdict(recon: ReconstructionEngine, physics: PhysicsEngine,
                          cfg: Dict, stress: StressProfile, theta_gate: float):
    """Run the exact Ry(theta) channel audit used by the CLI -> (D, dlog, dinc, verdict)."""
    c_phys = recon.perform_qpt_pipeline(theta_gate, stress)
    c_ideal = recon.ideal_choi_ry(theta_gate)
    d_val, c_proj = recon.process_distance_d_physical(c_ideal, c_phys)
    fit = recon.fit_unitary_plus_noise_from_choi(c_proj)
    dlog = recon.delta_log(physics.ry(theta_gate), fit["U_eff"])
    dinc = recon.delta_inc(fit["pauli_probs"])
    return d_val, dlog, dinc, classify_verdict(d_val, dlog, dinc, cfg)


def stressed_ry_kraus(recon: ReconstructionEngine, stress: StressProfile,
                      theta_gate: float) -> List[qt.Qobj]:
    """CPTP-projected Kraus operators of the physically-implemented Ry(theta) channel
    under `stress`. Reuses perform_qpt_pipeline + project_choi_to_cptp_1q +
    kraus_operators_from_choi (all existing functions)."""
    c_phys = recon.perform_qpt_pipeline(theta_gate, stress)
    c_proj = recon.project_choi_to_cptp_1q(c_phys)
    return recon.kraus_operators_from_choi(c_proj)


def apply_kraus(rho: qt.Qobj, kraus: List[qt.Qobj]) -> qt.Qobj:
    out = qt.Qobj(np.zeros((2, 2), dtype=complex))
    for k in kraus:
        out = out + k * rho * k.dag()
    return out


# ---- Single-qubit Clifford group (24 elements, generated from H and S) ----
def _h() -> qt.Qobj:
    return qt.Qobj(np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2.0))


def _s() -> qt.Qobj:
    return qt.Qobj(np.array([[1, 0], [0, 1j]], dtype=complex))


def _canon(u: np.ndarray) -> tuple:
    """Canonical key for a 2x2 unitary up to global phase (round to grid)."""
    # remove global phase by fixing first nonzero entry to be real positive
    flat = u.flatten()
    idx = np.argmax(np.abs(flat) > 1e-9)
    phase = np.exp(-1j * np.angle(flat[idx]))
    v = (u * phase).round(6)
    return tuple(v.flatten())


def single_qubit_cliffords() -> List[qt.Qobj]:
    """Return the 24 single-qubit Clifford unitaries (up to global phase)."""
    H, S = _h(), _s()
    gens = [H, S]
    found: Dict[tuple, qt.Qobj] = {}
    I = qt.qeye(2)
    frontier = [I]
    found[_canon(I.full())] = I
    while frontier:
        new_frontier = []
        for u in frontier:
            for g in gens:
                w = g * u
                key = _canon(w.full())
                if key not in found:
                    found[key] = qt.Qobj(w.full())
                    new_frontier.append(found[key])
        frontier = new_frontier
        if len(found) >= 24:
            break
    cliffs = list(found.values())
    assert len(cliffs) == 24, f"expected 24 Cliffords, got {len(cliffs)}"
    return cliffs
