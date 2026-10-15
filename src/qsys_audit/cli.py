"""Command-line entry point for qsys_audit (console script: ``qsys-audit``)."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Optional

import yaml

from . import __version__
from .paths import default_config

WORKLOAD_ALIASES = (
    "ghz",
    "graph_state",
    "w_state",
    "qft",
    "qaoa",
    "vqe_su2",
    "vqe_two_local",
)


def load_config(config_path: Optional[str] = None) -> Dict:
    """Load the YAML configuration, resolving the default location if needed."""
    path = Path(config_path) if config_path else default_config()
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="qsys-audit",
        description=(
            "Offline cross-layer workload auditor for quantum operating systems "
            "and schedulers."
        ),
    )
    parser.add_argument(
        "--config",
        default=None,
        help=(
            "Path to the YAML configuration file. Defaults to $QSYS_AUDIT_CONFIG, "
            "then ./config.yaml, then the config shipped with the package."
        ),
    )
    parser.add_argument(
        "--workload",
        required=True,
        help="Workload alias, for example: " + ", ".join(WORKLOAD_ALIASES) + ".",
    )
    parser.add_argument(
        "--stress",
        required=True,
        choices=["calibrated", "moderate", "aggressive"],
        help="Stress profile for concurrent hardware interference.",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Emit JSON output instead of text.",
    )
    parser.add_argument("--version", action="version", version=f"qsys_audit {__version__}")
    return parser


def classify_verdict(d_val: float, dlog: float, dinc: float, cfg: Dict) -> str:
    """Three-way channel verdict from the thresholds in config.yaml."""
    thresholds = cfg.get("verification_thresholds", {})
    if (
        d_val <= thresholds.get("eps_d_pass", 0.02)
        and dlog <= thresholds.get("eps_log_pass", 0.05)
        and dinc <= thresholds.get("eps_inc_pass", 0.05)
    ):
        return "PASS"

    if (
        d_val >= thresholds.get("d_hard", 0.30)
        or dlog >= thresholds.get("log_hard", 0.90)
        or dinc >= thresholds.get("inc_hard", 0.50)
    ):
        return "HARD FAIL"

    return "SOFT FAIL"


def run_cli(args: argparse.Namespace) -> Dict:
    # Defer the heavy scientific imports so --help stays available before install.
    from .physics_engine import PhysicsEngine
    from .reconstruction import ReconstructionEngine
    from .workload_replay import WorkloadReplay

    cfg = load_config(args.config)
    physics = PhysicsEngine(cfg)
    recon = ReconstructionEngine(physics, cfg)
    replay = WorkloadReplay(physics, cfg)

    profiles = physics.stress_profiles()
    if args.stress not in profiles:
        raise ValueError(f"Unknown stress profile: {args.stress}")

    metrics = replay.benchmark_workload(args.workload, profiles[args.stress])

    # The single-qubit channel reconstruction is included for cross-layer reporting.
    theta_gate = float(
        cfg.get("single_qubit_gate", {}).get("theta_gate", 1.5707963267948966)
    )
    c_phys = recon.perform_qpt_pipeline(theta_gate, profiles[args.stress])
    c_ideal = recon.ideal_choi_ry(theta_gate)
    d_val, c_proj = recon.process_distance_d_physical(c_ideal, c_phys)
    fit = recon.fit_unitary_plus_noise_from_choi(c_proj)
    dlog = recon.delta_log(physics.ry(theta_gate), fit["U_eff"])
    dinc = recon.delta_inc(fit["pauli_probs"])

    verdict = classify_verdict(d_val, dlog, dinc, cfg)

    return {
        "workload": metrics.workload,
        "stress": metrics.stress,
        "expected_fidelity": metrics.f_expected,
        "worst_case_fidelity": metrics.f_worst,
        "fidelity_spread": metrics.f_spread,
        "n_cx": metrics.n_cx,
        "n_1q": metrics.n_1q,
        "depth": metrics.depth,
        "size": metrics.size,
        "process_distance_d": d_val,
        "delta_log": dlog,
        "delta_inc": dinc,
        "verdict": verdict,
    }


def main(argv: Optional[list] = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    result = run_cli(args)

    if args.json:
        print(json.dumps(result, indent=2, sort_keys=True))
        return

    print(f"Workload: {result['workload']}")
    print(f"Stress: {result['stress']}")
    print(f"Expected Fidelity: {result['expected_fidelity']:.6f}")
    print(f"Worst-Case Fidelity: {result['worst_case_fidelity']:.6f}")
    print(f"Fidelity Spread: {result['fidelity_spread']:.6f}")
    print(f"Process Infidelity D: {result['process_distance_d']:.6f}")
    print(f"delta_log: {result['delta_log']:.6f}")
    print(f"delta_inc: {result['delta_inc']:.6f}")
    print(f"Verdict: {result['verdict']}")


if __name__ == "__main__":  # pragma: no cover
    main()
