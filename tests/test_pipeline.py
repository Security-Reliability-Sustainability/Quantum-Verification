"""Smoke tests for the physical pipeline.

These need the scientific stack and are skipped when it is not installed.
The reference numbers are the calibrated-regime row of docs/VERIFY.md.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

pytest.importorskip("qutip", reason="QuTiP not installed")
pytest.importorskip("scipy", reason="SciPy not installed")

import yaml  # noqa: E402

from qsys_audit.cli import classify_verdict  # noqa: E402


@pytest.fixture(scope="module")
def engines():
    from qsys_audit.physics_engine import PhysicsEngine
    from qsys_audit.reconstruction import ReconstructionEngine

    cfg = yaml.safe_load((ROOT / "config.yaml").read_text(encoding="utf-8"))
    physics = PhysicsEngine(cfg)
    return cfg, physics, ReconstructionEngine(physics, cfg)


def test_calibrated_channel_passes(engines):
    cfg, physics, recon = engines
    theta = float(cfg["single_qubit_gate"]["theta_gate"])
    stress = physics.stress_profiles()["calibrated"]

    c_phys = recon.perform_qpt_pipeline(theta, stress)
    c_ideal = recon.ideal_choi_ry(theta)
    d_val, c_proj = recon.process_distance_d_physical(c_ideal, c_phys)
    fit = recon.fit_unitary_plus_noise_from_choi(c_proj)
    dlog = recon.delta_log(physics.ry(theta), fit["U_eff"])
    dinc = recon.delta_inc(fit["pauli_probs"])

    assert d_val < 0.01, f"calibrated regime should be near ideal, got D={d_val}"
    assert classify_verdict(d_val, dlog, dinc, cfg) == "PASS"


def test_aggressive_channel_hard_fails(engines):
    cfg, physics, recon = engines
    theta = float(cfg["single_qubit_gate"]["theta_gate"])
    stress = physics.stress_profiles()["aggressive"]

    c_phys = recon.perform_qpt_pipeline(theta, stress)
    c_ideal = recon.ideal_choi_ry(theta)
    d_val, c_proj = recon.process_distance_d_physical(c_ideal, c_phys)
    fit = recon.fit_unitary_plus_noise_from_choi(c_proj)
    dlog = recon.delta_log(physics.ry(theta), fit["U_eff"])
    dinc = recon.delta_inc(fit["pauli_probs"])

    assert d_val == pytest.approx(0.3377, abs=5e-3)
    assert classify_verdict(d_val, dlog, dinc, cfg) == "HARD FAIL"


def test_workload_replay_worst_case_fidelity():
    pytest.importorskip("qiskit", reason="Qiskit not installed")
    pytest.importorskip("mqt.bench", reason="MQT Bench not installed")

    from qsys_audit.physics_engine import PhysicsEngine
    from qsys_audit.workload_replay import WorkloadReplay

    cfg = yaml.safe_load((ROOT / "config.yaml").read_text(encoding="utf-8"))
    physics = PhysicsEngine(cfg)
    replay = WorkloadReplay(physics, cfg)
    metrics = replay.benchmark_workload("ghz", physics.stress_profiles()["aggressive"])

    assert metrics.f_worst == pytest.approx(0.393, abs=5e-3)
    assert 0.0 <= metrics.f_worst <= metrics.f_expected <= 1.0
