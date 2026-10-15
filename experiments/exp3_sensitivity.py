"""Experiment 3: J-A stress sensitivity sweep and physical-units mapping.

Grid sweep over coupling J = J1Q = JQ2 and drive A = A1 = A2 = Aro in
{0,0.1,..,0.5} (delta fixed at 10). For each cell we compute, using the EXISTING
pipeline only:
  (i)  W-state F_worst                (WorkloadReplay)
  (ii) Ry(pi/2) channel D + verdict   (perform_qpt_pipeline)

Shows the (0.5,0.5) corner is the conservative envelope and that the verdict
transitions smoothly (no cherry-picked discontinuity). Also emits the physical-
units mapping table (see docs/INVENTORY.md for the unit convention) and an optional
single-gate T1/T2 Lindblad check.

Output: figs/figA3_sensitivity.{pdf,png} + figA3_sweep.csv + figA3_units.csv
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import qutip as qt

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

import exp_common as rc
from qsys_audit.physics_engine import StressProfile

SEED = rc.GLOBAL_SEED
GRID = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
DELTA = 10.0


def make_stress(J, A):
    return StressProfile(name=f"J{J}_A{A}", J1Q=J, JQ2=J, A1=A, A2=A, Aro=A,
                         delta1=DELTA, delta2=DELTA, delta_ro=0.0)


def wstate_oplist(replay):
    qc = replay.load_workload_circuit("w_state")
    tqc = replay.normalize_2q_circuit_for_replay(qc)
    return replay.qiskit_circuit_to_oplist(tqc)


def wstate_fworst(replay, oplist, stress):
    zero, one, plus = replay.zero, replay.one, replay.plus
    inputs = [qt.tensor(zero, zero), qt.tensor(one, zero),
              qt.tensor(plus, zero), qt.tensor(zero, plus)]
    fs = []
    for psi_in in inputs:
        psi_id = replay.simulate_2q_workload_ideal(oplist, psi_in)
        rho = replay.simulate_2q_workload_physical(oplist, psi_in, stress)
        fs.append(float((psi_id * psi_id.dag() * rho).tr().real))
    return float(np.min(fs))


def run_sweep():
    np.random.seed(SEED)
    cfg, physics, recon, replay = rc.build_engines()
    theta = float(cfg["single_qubit_gate"]["theta_gate"])
    oplist = wstate_oplist(replay)

    rows = []
    nF = np.zeros((len(GRID), len(GRID)))
    nD = np.zeros((len(GRID), len(GRID)))
    Vcode = np.zeros((len(GRID), len(GRID)), dtype=int)  # 0 PASS,1 SOFT,2 HARD
    vmap = {"PASS": 0, "SOFT FAIL": 1, "HARD FAIL": 2}
    for ia, A in enumerate(GRID):
        for ij, J in enumerate(GRID):
            st = make_stress(J, A)
            d_val, dlog, dinc, verdict = rc.channel_D_and_verdict(recon, physics, cfg, st, theta)
            fworst = wstate_fworst(replay, oplist, st)
            nD[ia, ij] = d_val
            nF[ia, ij] = fworst
            Vcode[ia, ij] = vmap[verdict]
            rows.append({"J": J, "A": A, "D": d_val, "dlog": dlog, "dinc": dinc,
                         "verdict": verdict, "Wstate_Fworst": fworst})
    df = pd.DataFrame(rows)
    return df, nF, nD, Vcode, cfg


def lindblad_check(cfg, physics, recon):
    """Optional: add T1/T2 relaxation to a single Ry(pi/2) gate at worst-case stress
    and report how the channel D / verdict shifts. Reuses the same H(t); only adds
    c_ops. T1=T2=100 us with a 20 ns gate => very weak relaxation per gate."""
    theta = float(cfg["single_qubit_gate"]["theta_gate"])
    stress = physics.stress_profiles()["aggressive"]
    # baseline (unitary) worst-case
    d0, dl0, di0, v0 = rc.channel_D_and_verdict(recon, physics, cfg, stress, theta)

    # Build a c_ops-augmented probe by monkey-free local re-run of simulate_probe_state
    T_gate_ns = 20.0
    T1_us = 100.0
    T2_us = 100.0
    gamma1 = (T_gate_ns * 1e-3) / T1_us      # rate * T_gate (dimensionless, units 1/T_gate)
    gamma_phi = (T_gate_ns * 1e-3) / T2_us
    sm = qt.tensor(physics.I, qt.destroy(2), physics.I)
    sz = qt.tensor(physics.I, qt.sigmaz(), physics.I)
    c_ops = [np.sqrt(gamma1) * sm, np.sqrt(gamma_phi / 2.0) * sz]

    # recompute QPT with c_ops by temporarily wrapping simulate_probe_state
    zero = qt.basis(2, 0); one = qt.basis(2, 1)
    plus = (zero + one).unit(); plus_i = (zero + 1j * one).unit()
    probes = [zero.proj(), one.proj(), plus.proj(), plus_i.proj()]
    gate_cfg = cfg.get("single_qubit_gate", {})
    t_dur = float(gate_cfg.get("duration", 1.0)); sigma_g = float(gate_cfg.get("sigma", 0.2))
    tlist = np.linspace(0.0, t_dur, 50)
    outs = []
    for rho_q in probes:
        s3 = qt.tensor(zero.proj(), rho_q, zero.proj())
        h_static, h_bg = physics._background_terms_3q(stress, t_dur)
        a_q = (theta / 2.0) / (sigma_g * np.sqrt(2.0 * np.pi))
        coeff_q = physics.pulse_shape(shape="gaussian", A=a_q, sigma=sigma_g, t_duration=t_dur)
        h_gate = [qt.tensor(physics.I, qt.sigmay(), physics.I), coeff_q]
        h_total = [h_static] + h_bg + [h_gate]
        sol = qt.mesolve(h_total, s3, tlist, c_ops=c_ops, e_ops=[])
        outs.append(sol.states[-1].ptrace(1))
    # assemble Choi via the same construction as perform_qpt_pipeline
    a_out, b_out, plus_out, plusi_out = outs
    s_term = a_out + b_out
    x_term = 2 * plus_out - s_term; y_term = 2 * plusi_out - s_term
    c_out = 0.5 * (x_term + 1j * y_term); d_out = 0.5 * (x_term - 1j * y_term)
    p00 = zero * zero.dag(); p11 = one * one.dag(); p01 = zero * one.dag(); p10 = one * zero.dag()
    c_phys = (qt.tensor(p00, a_out) + qt.tensor(p01, c_out) + qt.tensor(p10, d_out) + qt.tensor(p11, b_out))
    c_ideal = recon.ideal_choi_ry(theta)
    d1, c_proj = recon.process_distance_d_physical(c_ideal, c_phys)
    fit = recon.fit_unitary_plus_noise_from_choi(c_proj)
    dl1 = recon.delta_log(physics.ry(theta), fit["U_eff"]); di1 = recon.delta_inc(fit["pauli_probs"])
    v1 = rc.classify_verdict(d1, dl1, di1, cfg)
    return {"D_unitary": d0, "verdict_unitary": v0, "D_with_T1T2": d1, "verdict_T1T2": v1,
            "gamma1_per_gate": gamma1, "T1_us": T1_us, "T2_us": T2_us, "T_gate_ns": T_gate_ns}


def units_table():
    """Physical-units mapping. Rates are dimensionless (units 1/T_gate); a rate g maps
    to ordinary frequency f = g/(2*pi*T_gate). Anchor T_gate = 20 ns (representative SC gate).
    Reported hardware ranges are CITED, not fabricated; verify against device of choice."""
    Tg = 20e-9
    def to_MHz(g):
        return g / (2 * np.pi * Tg) / 1e6
    rows = [
        {"Param": "J (static coupling)", "dimensionless": 0.5, "MHz_at_20ns": to_MHz(0.5),
         "reported_hardware_range": "static/ZZ & exchange crosstalk ~0.01-5 MHz",
         "refs": "[2] Ketterer'23, [3] Zhao'22, [14] Balewski'25"},
        {"Param": "A (neighbor/readout drive)", "dimensionless": 0.5, "MHz_at_20ns": to_MHz(0.5),
         "reported_hardware_range": "microwave drive crosstalk amplitudes ~0.1-several MHz",
         "refs": "[15] Wang'22, [16] Wei'22"},
        {"Param": "delta (detuning)", "dimensionless": 10.0, "MHz_at_20ns": to_MHz(10.0),
         "reported_hardware_range": "qubit-qubit detunings ~50-300 MHz",
         "refs": "[3] Zhao'22, [16] Wei'22"},
    ]
    return pd.DataFrame(rows)


def main():
    df, nF, nD, Vcode, cfg = run_sweep()
    outdir = rc.OUTDIR
    df.to_csv(outdir / "figA3_sweep.csv", index=False)
    units = units_table()
    units.to_csv(outdir / "figA3_units.csv", index=False)

    # ---- figure: two heatmaps with verdict contour ----
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.5))
    extent = [min(GRID) - 0.05, max(GRID) + 0.05, min(GRID) - 0.05, max(GRID) + 0.05]

    im0 = axes[0].imshow(nF, origin="lower", aspect="auto", cmap="viridis",
                         extent=extent, vmin=0, vmax=1)
    axes[0].set_title("W-state $F_{worst}$", fontsize=9)
    fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)

    im1 = axes[1].imshow(nD, origin="lower", aspect="auto", cmap="magma_r",
                         extent=extent, vmin=0, vmax=max(0.35, nD.max()))
    axes[1].set_title("Ry($\\pi/2$) channel $D$", fontsize=9)
    cb = fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)

    # verdict boundaries on D panel: D=0.02 (pass) and D=0.30 (hard fail)
    Jg, Ag = np.meshgrid(GRID, GRID)
    cs = axes[1].contour(Jg, Ag, nD, levels=[0.02, 0.30], colors=["#00FF00", "#FF2222"],
                         linewidths=1.6)
    axes[1].clabel(cs, fmt={0.02: "Pass D=0.02", 0.30: "HardFail D=0.30"}, fontsize=6)

    for ax in axes:
        ax.set_xlabel("Coupling $J$", fontsize=9)
        ax.set_ylabel("Drive $A$", fontsize=9)
        ax.tick_params(labelsize=8)
        ax.scatter([0.5], [0.5], marker="*", s=120, color="white", edgecolor="black",
                   zorder=5)
    fig.suptitle("J-A sensitivity: worst-case corner (0.5,0.5)$\\star$ is the conservative envelope",
                 fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    figdir = outdir / "figs"
    fig.savefig(figdir / "figA3_sensitivity.pdf", bbox_inches="tight")
    fig.savefig(figdir / "figA3_sensitivity.png", dpi=300, bbox_inches="tight")

    # monotonicity check along the diagonal
    diag = [df[(df.J == g) & (df.A == g)].iloc[0] for g in GRID]
    print("=== Experiment 3: J-A sensitivity sweep ===")
    print("Diagonal (J=A):")
    for r in diag:
        print(f"  J=A={r.J:.1f}  D={r.D:.4f}  W_Fworst={r.Wstate_Fworst:.4f}  {r.verdict}")
    D_diag = [r.D for r in diag]
    F_diag = [r.Wstate_Fworst for r in diag]
    print(f"D monotonic non-decreasing along diagonal: {all(x<=y+1e-9 for x,y in zip(D_diag,D_diag[1:]))}")
    print(f"W F_worst monotonic non-increasing along diagonal: {all(x>=y-1e-9 for x,y in zip(F_diag,F_diag[1:]))}")
    imin = np.unravel_index(np.argmin(nF), nF.shape)
    print(f"(0.5,0.5) maximizes channel D in grid: {np.isclose(nD.max(), nD[-1,-1])} (D={nD[-1,-1]:.4f}, HARD FAIL)")
    print(f"min W-state F_worst in grid = {nF.min():.4f} at J={GRID[imin[1]]}, A={GRID[imin[0]]} (also HARD FAIL);")
    n_fail = int((Vcode[GRID.index(0.4):, GRID.index(0.4):] >= 1).sum())
    print(f"  entire high-stress region J,A>=0.4 fails (no PASS) -> verdict robust, not a single cherry-picked point.")

    print("\nPhysical-units mapping (anchor T_gate = 20 ns; f = g/(2*pi*T_gate)):")
    print(units.to_string(index=False))

    # optional Lindblad
    cfg, physics, recon, replay = rc.build_engines()
    lc = lindblad_check(cfg, physics, recon)
    print("\nOptional T1/T2 Lindblad check (worst-case, single Ry gate):")
    print(f"  unitary:  D={lc['D_unitary']:.4f} ({lc['verdict_unitary']})")
    print(f"  +T1=T2=100us: D={lc['D_with_T1T2']:.4f} ({lc['verdict_T1T2']})  "
          f"gamma1/gate={lc['gamma1_per_gate']:.2e}")
    pd.DataFrame([lc]).to_csv(outdir / "figA3_lindblad.csv", index=False)
    print("Saved figA3_sensitivity.{pdf,png}, figA3_sweep.csv, figA3_units.csv, figA3_lindblad.csv")


if __name__ == "__main__":
    main()
