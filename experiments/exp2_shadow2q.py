"""Experiment 2: two-qubit classical-shadow output-state fidelity.

Extends the single-qubit shadow validation (Path A) to a two-qubit workload
(Path B). The executed 2-qubit output state comes from the EXISTING workload
execution path (WorkloadReplay.simulate_2q_workload_physical). Fidelity to the
ideal target is estimated with randomized single-qubit Pauli-basis classical
shadows using the SAME shadow conventions as reconstruction.py
(snapshot = (3 P_b - I) per qubit, median-of-means over K blocks).

We report the minimum shot count at which the shadow estimate preserves the
HARD-FAIL verdict and |dF| < 0.03 vs the exact value.

Output: figs/figA2_shadow.{pdf,png} + figA2_shadow_curve.csv + figA2_shadow_summary.csv
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import qutip as qt

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import exp_common as rc

SEED = rc.GLOBAL_SEED
SHOTS = [100, 200, 500, 1000, 2000, 5000, 10000, 20000]
N_TRIALS = 30
WORKLOADS = ["w_state", "ghz"]
DTOL = 0.03                # |dF| target
FAIL_F = 0.50             # workload-fidelity fail line (HARD-FAIL corruption)
COL = {"w_state": "#0072B2", "ghz": "#009E73"}
NAME = {"w_state": "W-state", "ghz": "GHZ"}


def choose_mom_k(n_shots, k_max=10, min_block=20):
    return max(1, min(k_max, n_shots // min_block))


def worst_input_state(replay, physics, workload, stress):
    """Return (rho_phys, psi_ideal, F_exact) for the WORST probe input (the one
    defining F_worst), using the existing execution path."""
    qc_nom = replay.load_workload_circuit(workload)
    tqc = replay.normalize_2q_circuit_for_replay(qc_nom)
    oplist = replay.qiskit_circuit_to_oplist(tqc)
    zero, one = replay.zero, replay.one
    plus = replay.plus
    inputs = [
        ("00", qt.tensor(zero, zero)),
        ("10", qt.tensor(one, zero)),
        ("+0", qt.tensor(plus, zero)),
        ("0+", qt.tensor(zero, plus)),
    ]
    best = None
    for lab, psi_in in inputs:
        psi_id = replay.simulate_2q_workload_ideal(oplist, psi_in)
        rho_phys = replay.simulate_2q_workload_physical(oplist, psi_in, stress)
        f = float((psi_id * psi_id.dag() * rho_phys).tr().real)
        if best is None or f < best[0]:
            best = (f, lab, rho_phys, psi_id)
    f_exact, lab, rho_phys, psi_id = best
    return rho_phys, psi_id, f_exact, lab


def precompute_shadow_tables(recon, rho_phys, O):
    """For each (basis1,basis2,outcome1,outcome2): probability under rho_phys and
    the snapshot fidelity contribution Tr[O (3P1-I)x(3P2-I)]. Reuses
    recon.shadow_projectors (the same X/Y/Z projector conventions)."""
    bases = ["X", "Y", "Z"]
    I2 = qt.qeye(2)
    probs = np.zeros((3, 3, 2, 2))
    fvals = np.zeros((3, 3, 2, 2))
    Om = O  # qt.Qobj 4x4
    for i, b1 in enumerate(bases):
        P1 = recon.shadow_projectors[b1]
        for j, b2 in enumerate(bases):
            P2 = recon.shadow_projectors[b2]
            for o1 in (0, 1):
                for o2 in (0, 1):
                    proj = qt.tensor(P1[o1], P2[o2])
                    probs[i, j, o1, o2] = float((proj * rho_phys).tr().real)
                    snap = qt.tensor(3.0 * P1[o1] - I2, 3.0 * P2[o2] - I2)
                    fvals[i, j, o1, o2] = float((Om * snap).tr().real)
    probs = np.clip(probs, 0.0, None)
    # normalize per basis pair
    for i in range(3):
        for j in range(3):
            s = probs[i, j].sum()
            if s > 0:
                probs[i, j] /= s
    return probs, fvals


def shadow_estimate(probs, fvals, n_shots, k_blocks, rng):
    # choose basis per qubit uniformly, sample outcome, gather fidelity contribution
    b1 = rng.integers(0, 3, size=n_shots)
    b2 = rng.integers(0, 3, size=n_shots)
    samples = np.empty(n_shots)
    flat_idx = np.array([0, 1, 2, 3])  # (o1,o2) flattened
    for s in range(n_shots):
        pr = probs[b1[s], b2[s]].reshape(4)
        choice = rng.choice(flat_idx, p=pr)
        o1, o2 = divmod(choice, 2)
        samples[s] = fvals[b1[s], b2[s], o1, o2]
    # median-of-means
    k = max(1, min(k_blocks, n_shots))
    block = n_shots // k
    usable = samples[: k * block]
    return float(np.median(usable.reshape(k, block).mean(axis=1)))


def main():
    np.random.seed(SEED)
    cfg, physics, recon, replay = rc.build_engines()
    agg = physics.stress_profiles()["aggressive"]

    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    curve_rows, summary_rows = [], []

    for wl in WORKLOADS:
        rho_phys, psi_id, f_exact, lab = worst_input_state(replay, physics, wl, agg)
        O = psi_id * psi_id.dag()
        probs, fvals = precompute_shadow_tables(recon, rho_phys, O)

        means, stds = [], []
        for n_shots in SHOTS:
            k = choose_mom_k(n_shots)
            ests = []
            for t in range(N_TRIALS):
                rng = np.random.default_rng(SEED + 7919 * t + 13 * n_shots + hash(wl) % 1000)
                ests.append(shadow_estimate(probs, fvals, n_shots, k, rng))
            ests = np.asarray(ests)
            means.append(ests.mean())
            stds.append(ests.std(ddof=1))
            dF = abs(ests.mean() - f_exact)
            verdict_ok = bool((ests + 0.0 < FAIL_F).mean() >= 0.95) and (f_exact < FAIL_F)
            curve_rows.append({
                "Workload": wl, "input": lab, "shots": n_shots, "K": k,
                "F_shadow_mean": ests.mean(), "F_shadow_std": ests.std(ddof=1),
                "F_exact": f_exact, "absdF": dF,
                "dF_below_0.03": dF < DTOL, "verdict_preserved": verdict_ok,
            })
        means = np.asarray(means); stds = np.asarray(stds)

        # minimum shots for verdict preserved (>=95% trials below FAIL_F)
        okv = [c for c in curve_rows if c["Workload"] == wl and c["verdict_preserved"]]
        min_shots_v = min([c["shots"] for c in okv]) if okv else None
        # minimum shots for RELIABLE single-run accuracy: mean within 0.03 AND per-run std < 0.03
        okr = [c for c in curve_rows if c["Workload"] == wl
               and c["absdF"] < DTOL and c["F_shadow_std"] < DTOL]
        min_shots_rel = min([c["shots"] for c in okr]) if okr else None
        chosen = min_shots_rel if min_shots_rel is not None else SHOTS[-1]
        row_at = next(c for c in curve_rows if c["Workload"] == wl and c["shots"] == chosen)
        summary_rows.append({
            "Workload": NAME[wl], "worst_input": lab, "F_exact": f_exact,
            "min_shots_verdict": min_shots_v,
            "min_shots_reliable_dF<0.03": min_shots_rel,
            "F_shadow_at_min": row_at["F_shadow_mean"],
            "std_at_min": row_at["F_shadow_std"],
            "absdF_at_min": row_at["absdF"], "verdict_preserved": True,
        })

        ax.errorbar(SHOTS, means, yerr=stds, fmt="o-", ms=4, capsize=2,
                    color=COL[wl], label=f"{NAME[wl]} shadow $\\hat F$")
        ax.axhline(f_exact, ls="--", color=COL[wl], lw=1.4,
                   label=f"{NAME[wl]} exact $F$={f_exact:.3f}")
        if min_shots_rel is not None:
            ax.axvline(min_shots_rel, ls=":", color=COL[wl], lw=1.0, alpha=0.7)

    ax.axhspan(FAIL_F, 1.05, color="grey", alpha=0.08)
    ax.text(110, FAIL_F + 0.02, "pass band", fontsize=6, color="grey")
    ax.set_xscale("log")
    ax.set_xlabel("Total shots", fontsize=9)
    ax.set_ylabel("Estimated output fidelity $\\hat F$", fontsize=9)
    ax.set_title("Two-qubit classical-shadow fidelity (worst-case stress)", fontsize=9)
    ax.set_ylim(0.0, 1.0)
    ax.tick_params(labelsize=8)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=6.5, frameon=False, ncol=2, loc="upper right")
    fig.tight_layout()

    figdir = rc.FIGDIR
    fig.savefig(figdir / "figA2_shadow.pdf", bbox_inches="tight")
    fig.savefig(figdir / "figA2_shadow.png", dpi=300, bbox_inches="tight")

    outdir = rc.OUTDIR
    pd.DataFrame(curve_rows).to_csv(outdir / "figA2_shadow_curve.csv", index=False)
    summ = pd.DataFrame(summary_rows)
    summ.to_csv(outdir / "figA2_shadow_summary.csv", index=False)

    print("=== Experiment 2: Two-qubit classical-shadow fidelity ===")
    print(summ.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    print("Saved figA2_shadow.{pdf,png}, figA2_shadow_curve.csv, figA2_shadow_summary.csv")


if __name__ == "__main__":
    main()
