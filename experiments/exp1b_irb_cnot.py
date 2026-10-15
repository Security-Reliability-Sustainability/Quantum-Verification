"""Experiment 1b: interleaved RB on the echoed CR-CNOT, which ties the RB
measurement to the two-qubit workloads.

Mirrors the paper's workload-execution convention exactly: single-qubit gates are
IDEAL (noiseless), and the CR-CNOT is the physical, stressed
`evolve_echoed_cr_cnot`. We build the noisy CR-CNOT as a 16x16 superoperator (via
2-qubit process reconstruction over the 16 operator-basis inputs, q0 idle in |0>),
and run interleaved RB:

  * Reference sequence: m random 2-qubit LOCAL Cliffords (24x24 group, ideal) +
    ideal recovery  ->  p_ref ~ 1 (noiseless reference, by the model's convention).
  * Interleaved sequence: each random local Clifford is followed by the physical
    CR-CNOT  ->  decay p_int.
  * CR-CNOT average error: r = (d-1)/d * (1 - p_int/p_ref), d = 4.

The random local-Clifford twirl differs every step, so the decay is a GENUINE
twirled average (not one fixed channel raised to a power); the fit R^2 is
meaningful.

Tie to workloads: RB reports the CR-CNOT AVERAGE gate fidelity. Extrapolating it
independently across a workload's CX count (F_CRCNOT^n_cx) predicts a fidelity that
is OVER-OPTIMISTIC versus the audit's worst-case F_worst for the low-depth
entangling workloads (GHZ, W-state) -- the coherent worst-case penalty RB misses.

Output: figs/figA1_irb.{pdf,png} + figA1_irb_curve.csv + figA1_irb_summary.csv
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import qutip as qt
from scipy.optimize import curve_fit

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import exp_common as rc

SEED = rc.GLOBAL_SEED
LENGTHS = [1, 2, 4, 8, 16, 32]
N_SEQ = 30
REGIMES = ["moderate", "aggressive"]
COL = {"moderate": "#E69F00", "aggressive": "#0072B2"}
LABEL = {"moderate": "Moderate", "aggressive": "Worst-case"}
D = 4
FLOOR = 1.0 / D
U_CNOT = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]], dtype=complex)
WL = [("ghz", "GHZ"), ("w_state", "W-state"), ("qft", "QFT"),
      ("vqe_su2", "VQE-SU2"), ("vqe_two_local", "VQE-TwoLocal")]
WCOL = {"GHZ": "#0072B2", "W-state": "#D55E00", "QFT": "#009E73",
        "VQE-SU2": "#CC79A7", "VQE-TwoLocal": "#56B4E9"}


def _vec16(M):
    return M.reshape(16, order="F")


def _usuper(U):
    return np.kron(np.conj(U), U)


def cnot_superop(physics, cfg, stress):
    ccfg = cfg["cnot"]
    z = qt.basis(2, 0)
    S = np.zeros((16, 16), dtype=complex)
    for i in range(4):
        for j in range(4):
            E = np.zeros((4, 4), dtype=complex); E[i, j] = 1.0
            rho3 = qt.tensor(z.proj(), qt.Qobj(E, dims=[[2, 2], [2, 2]]))
            out = physics.evolve_echoed_cr_cnot(
                rho3, stress=stress, theta_zx=float(ccfg["theta_zx"]),
                rx_target_corr=float(ccfg["rx_target_corr"]),
                rz_control_corr=float(ccfg["rz_control_corr"]),
                beta_ix=float(ccfg.get("beta_ix", 0.2)), gamma_zi=float(ccfg.get("gamma_zi", 0.05)),
                sigma=float(ccfg.get("sigma", 0.25)), t_duration=float(ccfg.get("duration", 1.5)),
                nsteps=int(ccfg.get("nsteps", 120)))
            S[:, j * 4 + i] = _vec16(out.ptrace([1, 2]).full())
    return S


def local_clifford_superops():
    cliffs = [c.full() for c in rc.single_qubit_cliffords()]  # 24 ideal 1q
    sup, uni = [], []
    for a in cliffs:
        for b in cliffs:
            U = np.kron(a, b)
            uni.append(U); sup.append(_usuper(U))
    return sup, uni  # 576 local 2q Cliffords


def _fit(lengths, surv):
    def f(m, A, p, B):
        return A * p ** m + B
    try:
        po, _ = curve_fit(f, lengths, surv, p0=[surv[0] - FLOOR, 0.6, FLOOR],
                          bounds=([0, 0, 0], [1.5, 1, 1]), maxfev=40000)
    except Exception:
        po = [0.7, 0.6, FLOOR]
    yh = f(np.array(lengths, float), *po)
    ss_tot = float(np.sum((surv - surv.mean()) ** 2))
    r2 = 1 - float(np.sum((surv - yh) ** 2)) / ss_tot if ss_tot > 1e-12 else float("nan")
    return po[0], po[1], po[2], r2


def run_irb(Scnot, locsup, locU, interleave, seed):
    rng = np.random.default_rng(seed)
    E00 = np.zeros((4, 4), dtype=complex); E00[0, 0] = 1.0
    rho0 = _vec16(E00); proj00 = _vec16(E00).conj()
    nloc = len(locsup)
    surv = []
    for m in LENGTHS:
        ss = []
        for _ in range(N_SEQ):
            S = np.eye(16, dtype=complex); Un = np.eye(4, dtype=complex)
            for _k in range(m):
                li = int(rng.integers(0, nloc))
                S = locsup[li] @ S; Un = locU[li] @ Un
                if interleave:
                    S = Scnot @ S; Un = U_CNOT @ Un
            S = _usuper(np.conj(Un).T) @ S          # ideal recovery
            ss.append(float(np.real(np.dot(proj00, S @ rho0))))
        surv.append(float(np.mean(ss)))
    return np.asarray(surv)


def main():
    np.random.seed(SEED)
    cfg, physics, recon, replay = rc.build_engines()
    locsup, locU = local_clifford_superops()

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(8.8, 3.5))
    summary, curve_rows = [], []
    rcnot_byreg = {}

    for reg in REGIMES:
        st = physics.stress_profiles()[reg]
        Scnot = cnot_superop(physics, cfg, st)
        s_ref = run_irb(Scnot, locsup, locU, interleave=False, seed=SEED)
        s_int = run_irb(Scnot, locsup, locU, interleave=True, seed=SEED + 1)
        _, p_ref, _, r2_ref = _fit(LENGTHS, s_ref)
        Aint, p_int, Bint, r2_int = _fit(LENGTHS, s_int)
        r_cnot = (1 - p_int / p_ref) * (D - 1) / D
        f_cnot = 1 - r_cnot
        rcnot_byreg[reg] = (r_cnot, f_cnot)
        for m, sr, si in zip(LENGTHS, s_ref, s_int):
            curve_rows.append({"regime": reg, "m": m, "survival_ref": sr, "survival_int": si})
        summary.append({"regime": reg, "p_ref": p_ref, "p_int": p_int, "fit_R2_int": r2_int,
                        "r_CRCNOT_irb": r_cnot, "F_CRCNOT_irb": f_cnot})

        def f(m, A, p, B):
            return A * p ** m + B
        mm = np.array(LENGTHS, float); mf = np.linspace(min(mm), max(mm), 200)
        axL.plot(mm, s_int, "o", ms=4, color=COL[reg])
        axL.plot(mf, f(mf, Aint, p_int, Bint), "-", color=COL[reg], lw=2,
                 label=f"{LABEL[reg]} interleaved: $r_{{CR}}$={r_cnot:.3f}, $R^2$={r2_int:.3f}")
        axL.plot(mm, s_ref, "s:", ms=3, color=COL[reg], alpha=0.55,
                 label=f"{LABEL[reg]} reference ($p$={p_ref:.3f})")

    axL.axhline(FLOOR, ls=":", color="grey", lw=1)
    axL.text(1.05, FLOOR + 0.01, "floor $1/d$", fontsize=6, color="grey")
    axL.set_xscale("log", base=2)
    axL.set_xlabel("Sequence length $m$ (Cliffords)", fontsize=9)
    axL.set_ylabel("Survival probability", fontsize=9)
    axL.set_title("(a) Interleaved RB of the echoed CR-CNOT", fontsize=9)
    axL.set_ylim(0.2, 1.02); axL.tick_params(labelsize=8); axL.grid(alpha=0.25)
    axL.legend(fontsize=6.3, frameon=False)

    # ---- panel (b): RB per-CNOT extrapolation vs actual workload fidelity (worst-case) ----
    agg = physics.stress_profiles()["aggressive"]
    f_cnot_agg = rcnot_byreg["aggressive"][1]
    names, ncx, pred, fexp, fworst = [], [], [], [], []
    for key, disp in WL:
        m = replay.benchmark_workload(key, agg)
        names.append(disp); ncx.append(m.n_cx)
        pred.append(f_cnot_agg ** m.n_cx)
        fexp.append(m.f_expected); fworst.append(m.f_worst)
        curve_rows.append({"regime": "aggressive", "m": -1,
                           "workload": disp, "n_cx": m.n_cx,
                           "RB_pred_fid": f_cnot_agg ** m.n_cx,
                           "F_expected": m.f_expected, "F_worst": m.f_worst})
    x = np.arange(len(names)); w = 0.27
    axR.bar(x - w, pred, w, color="black", alpha=0.65, label="RB-extrapolated $F_{CR}^{n_{cx}}$")
    axR.bar(x, fexp, w, color="#999999", label="actual $F_{expected}$")
    axR.bar(x + w, fworst, w, color=[WCOL[n] for n in names], edgecolor="black", lw=0.4,
            label="actual $F_{worst}$ (audit)")
    axR.set_xticks(x); axR.set_xticklabels([f"{n}\n($n_{{cx}}$={c})" for n, c in zip(names, ncx)],
                                           fontsize=6.5)
    axR.set_ylabel("Fidelity", fontsize=9)
    axR.set_title("(b) RB per-CNOT extrapolation vs workload audit", fontsize=9)
    axR.set_ylim(0, 0.7); axR.tick_params(labelsize=8); axR.grid(axis="y", alpha=0.25)
    axR.legend(fontsize=6.3, frameon=False)
    fig.tight_layout()

    figdir = rc.FIGDIR
    fig.savefig(figdir / "figA1_irb.pdf", bbox_inches="tight")
    fig.savefig(figdir / "figA1_irb.png", dpi=300, bbox_inches="tight")

    outdir = rc.OUTDIR
    pd.DataFrame(curve_rows).to_csv(outdir / "figA1_irb_curve.csv", index=False)
    summ = pd.DataFrame(summary)
    summ.to_csv(outdir / "figA1_irb_summary.csv", index=False)

    print("=== Experiment 1b: interleaved RB on echoed CR-CNOT ===")
    print(summ.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print("\nWorst-case: RB per-CNOT extrapolation vs actual workload fidelity:")
    for n, c, pr, fe, fw in zip(names, ncx, pred, fexp, fworst):
        print(f"  {n:14s} n_cx={c}  RB_pred={pr:.3f}  F_expected={fe:.3f}  F_worst={fw:.3f}")
    print("Saved figA1_irb.{pdf,png}, figA1_irb_curve.csv, figA1_irb_summary.csv")


if __name__ == "__main__":
    main()
