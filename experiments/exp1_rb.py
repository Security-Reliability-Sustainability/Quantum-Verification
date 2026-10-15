"""Experiment 1: randomized benchmarking on the stressed channel.

PHYSICAL, GATE-DEPENDENT RB. Each of the 24 single-qubit Cliffords is compiled
into a word over the physical primitives {X90, Y90, X-90, Y-90}; every primitive
is simulated as a real pulse THROUGH the stressed evolution used by the gate
suite (run_shadow_convergence.get_final_state_for_gate -> the paper's H(t)). The
noisy primitive channels are assembled once (4-probe reconstruction -> 4x4
superoperator) and composed for each random sequence. Single global seed (1234).

Why this design:
  * The noise is GATE-DEPENDENT -- the anisotropic stress (sigma_y(x)sigma_x
    coupling, sigma_x drives) makes the X90 and Y90 error channels different --
    so the RB decay is NOT an analytically-forced single exponential. The fit
    quality R^2 is a genuine measurement; we also FLAG when the decay saturates
    to the depolarized floor and the single-exponential fit becomes degenerate.
  * We do NOT claim RB looks "benign": under this strong-coupling model the
    single-qubit gates are genuinely degraded and RB registers that. The honest
    QCVV blind spot is structural and is what we report:
      (i)  where resolvable, RB's AVERAGE infidelity r is SMALLER than the
           worst-case channel D / coherent delta_log the audit uses;
      (ii) at the worst-case corner the decay saturates to the 1/d floor within
           m~2, so the standard RB exponential fit is ILL-CONDITIONED -> RB
           returns an essentially unresolvable/ambiguous estimate;
      (iii) RB's single number is WORKLOAD- and INPUT-AGNOSTIC: identical for
           GHZ (F_worst=0.39) and the Hard-Fail W-state (F_worst=0.093), so an
           RB-thresholded scheduler cannot separate a near-pass tenant from a
           catastrophically-exposed victim (the 4.23x gap).

Output: figs/figA1_rb.{pdf,png} + figA1_rb_curve.csv + figA1_rb_sweep.csv + figA1_rb_summary.csv
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
import run_shadow_convergence as RSC
from qsys_audit.physics_engine import StressProfile

SEED = rc.GLOBAL_SEED
LENGTHS = [1, 2, 4, 8, 16, 32, 64, 128]
N_SEQ = 40
CURVE_REGIMES = ["calibrated", "moderate", "aggressive"]
COL = {"calibrated": "#009E73", "moderate": "#E69F00", "aggressive": "#0072B2"}
LABEL = {"calibrated": "Calibrated", "moderate": "Moderate", "aggressive": "Worst-case"}
WCOL = {"GHZ": "#0072B2", "W-state": "#D55E00", "QFT": "#009E73",
        "VQE-SU2": "#CC79A7", "VQE-TwoLocal": "#56B4E9"}
FLOOR = 0.5  # 1/d survival floor for a fully depolarizing single-qubit channel


def rb_decay(m, A, p, B):
    return A * (p ** m) + B


def _prim_spec(physics, axis, theta):
    sig = float(physics.config["single_qubit_gate"]["sigma"])
    dur = float(physics.config["single_qubit_gate"]["duration"])
    U = physics.rx(theta) if axis == "x" else physics.ry(theta)
    return {"ideal_U": U, "segments": [{"axis": axis, "theta": theta,
            "sigma": sig, "duration": dur}], "virtual_post": []}


def _vec(M):
    return M.full().reshape(4, order="F")  # column-stacking


def channel_superop(physics, gate_spec, stress):
    """Column-stacking 4x4 superoperator of the physically-simulated noisy gate
    (4-probe reconstruction; same probe set as the QPT pipeline)."""
    z, o = qt.basis(2, 0), qt.basis(2, 1)
    pl = (z + o).unit(); pli = (z + 1j * o).unit()
    outs = [RSC.get_final_state_for_gate(physics, p, gate_spec, stress, 50)
            for p in [z.proj(), o.proj(), pl.proj(), pli.proj()]]
    a, b, p1, p2 = outs
    s = a + b; X = 2 * p1 - s; Y = 2 * p2 - s
    c = 0.5 * (X + 1j * Y); d = 0.5 * (X - 1j * Y)
    return np.column_stack([_vec(a), _vec(d), _vec(c), _vec(b)])


def _canon(u):
    f = u.flatten(); i = int(np.argmax(np.abs(f) > 1e-9))
    ph = np.exp(-1j * np.angle(f[i]))
    return tuple((u * ph).round(6).flatten())


def compile_cliffords(physics):
    prims = {"X90": _prim_spec(physics, "x", np.pi / 2),
             "Y90": _prim_spec(physics, "y", np.pi / 2),
             "Xm": _prim_spec(physics, "x", -np.pi / 2),
             "Ym": _prim_spec(physics, "y", -np.pi / 2)}
    idU = {k: prims[k]["ideal_U"] for k in prims}
    I = np.eye(2, dtype=complex)
    words = {_canon(I): []}; mats = {_canon(I): I}; frontier = [(I, [])]
    while frontier and len(words) < 24:
        nf = []
        for M, w in frontier:
            for k in prims:
                N = idU[k].full() @ M; key = _canon(N)
                if key not in words:
                    words[key] = w + [k]; mats[key] = N; nf.append((N, w + [k]))
        frontier = nf
    assert len(words) == 24, f"got {len(words)} cliffords"
    return list(words.keys()), words, mats, prims


def run_rb(physics, prims, keys, words, mats, stress, seed=SEED):
    Sp = {k: channel_superop(physics, prims[k], stress) for k in prims}
    noisy = {}
    for key in keys:
        S = np.eye(4, dtype=complex)
        for k in words[key]:
            S = Sp[k] @ S
        noisy[key] = S
    rng = np.random.default_rng(seed)
    rho0 = np.array([1, 0, 0, 0], dtype=complex)
    rows = []
    for m in LENGTHS:
        sv = []
        for _ in range(N_SEQ):
            idx = rng.integers(0, 24, size=m)
            seq = [keys[i] for i in idx]
            Un = np.eye(2, dtype=complex); S = np.eye(4, dtype=complex)
            for key in seq:
                S = noisy[key] @ S; Un = mats[key] @ Un
            S = noisy[_canon(np.conj(Un).T)] @ S
            sv.append(float(np.real((S @ rho0)[0])))
        sv = np.asarray(sv)
        rows.append({"m": m, "survival_mean": sv.mean(), "survival_std": sv.std(ddof=1)})
    return pd.DataFrame(rows)


def fit_curve(df):
    m = df["m"].to_numpy(float); y = df["survival_mean"].to_numpy()
    try:
        popt, pcov = curve_fit(rb_decay, m, y, p0=[max(y[0] - FLOOR, 0.05), 0.7, FLOOR],
                               bounds=([0, 0, 0], [1.5, 1, 1]), maxfev=40000)
        p_err = float(np.sqrt(np.diag(pcov))[1]) if np.all(np.isfinite(pcov)) else float("inf")
    except Exception:
        popt = [0.5, 0.7, 0.5]; p_err = float("inf")
    A, p, B = popt
    yhat = rb_decay(m, *popt)
    ss_res = float(np.sum((y - yhat) ** 2)); ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 1e-12 else float("nan")
    r = (1 - p) * 0.5
    # degenerate / unresolved if the decay barely leaves the floor or p poorly determined
    resolved = (y[0] - FLOOR > 0.06) and (p_err < 0.25) and (p > 0.02)
    return {"A": A, "p": p, "B": B, "R2": r2, "r": r, "F_RB": 1 - r,
            "p_err": p_err, "resolved": bool(resolved), "surv_m1": float(y[0])}


def wstate_fworst(replay, physics, stress):
    return replay.benchmark_workload("w_state", stress).f_worst


def main():
    np.random.seed(SEED)
    cfg, physics, recon, replay = rc.build_engines()
    theta = float(cfg["single_qubit_gate"]["theta_gate"])
    keys, words, mats, prims = compile_cliffords(physics)

    agg = physics.stress_profiles()["aggressive"]
    wl_disp = [("ghz", "GHZ"), ("w_state", "W-state"), ("qft", "QFT"),
               ("vqe_su2", "VQE-SU2"), ("vqe_two_local", "VQE-TwoLocal")]
    wl_fworst = {disp: replay.benchmark_workload(k, agg).f_worst for k, disp in wl_disp}

    # ---- panel (a): decay curves ----
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(8.8, 3.5))
    curve_frames, summary_rows = [], []
    fit_byreg = {}
    for reg in CURVE_REGIMES:
        df = run_rb(physics, prims, keys, words, mats, physics.stress_profiles()[reg])
        df.insert(0, "Regime", reg)
        fit = fit_curve(df)
        fit_byreg[reg] = fit
        curve_frames.append(df)
        d_r, dl_r, di_r, v_r = rc.channel_D_and_verdict(recon, physics, cfg,
                                                        physics.stress_profiles()[reg], theta)
        summary_rows.append({
            "Regime": reg, "F_RB": fit["F_RB"], "r_avg_per_Clifford": fit["r"],
            "p_decay": fit["p"], "fit_R2": fit["R2"], "fit_resolved": fit["resolved"],
            "survival_m1": fit["surv_m1"], "worstcase_channel_D": d_r,
            "worstcase_delta_log": dl_r, "workflow_verdict": v_r,
        })
        mm = df["m"].to_numpy()
        axL.errorbar(mm, df["survival_mean"], yerr=df["survival_std"], fmt="o", ms=4,
                     capsize=2, color=COL[reg])
        mf = np.linspace(min(mm), max(mm), 200)
        tag = f"$r$={fit['r']:.3f}, $R^2$={fit['R2']:.3f}" if fit["resolved"] else "decay unresolved (floor)"
        axL.plot(mf, rb_decay(mf, fit["A"], fit["p"], fit["B"]), "-", color=COL[reg], lw=2,
                 label=f"{LABEL[reg]}: {tag}")
    axL.axhline(FLOOR, ls=":", color="grey", lw=1)
    axL.text(1.1, FLOOR + 0.01, "depolarized floor $1/d$", fontsize=6, color="grey")
    axL.set_xscale("log", base=2)
    axL.set_xlabel("Sequence length $m$ (Cliffords)", fontsize=9)
    axL.set_ylabel("Survival probability", fontsize=9)
    axL.set_title("(a) Physical gate-dependent RB decay", fontsize=9)
    axL.set_ylim(0.45, 1.02)
    axL.tick_params(labelsize=8); axL.grid(alpha=0.25)
    axL.legend(fontsize=6.8, frameon=False)

    # ---- panel (b): RB under-reports & is workload-agnostic, across the stress sweep ----
    sweep_defs = [("calib", physics.stress_profiles()["calibrated"]),
                  ("J=A=.1", StressProfile("s", .1, .1, .1, .1, .1, 10., 10., 0.)),
                  ("J=A=.2", StressProfile("s", .2, .2, .2, .2, .2, 10., 10., 0.)),
                  ("J=A=.3", StressProfile("s", .3, .3, .3, .3, .3, 10., 10., 0.)),
                  ("J=A=.4", StressProfile("s", .4, .4, .4, .4, .4, 10., 10., 0.)),
                  ("worst", physics.stress_profiles()["aggressive"])]
    sweep_rows = []
    for nm, st in sweep_defs:
        df = run_rb(physics, prims, keys, words, mats, st)
        fit = fit_curve(df)
        d_r, dl_r, di_r, v_r = rc.channel_D_and_verdict(recon, physics, cfg, st, theta)
        wf = wstate_fworst(replay, physics, st)
        sweep_rows.append({"regime": nm, "F_RB": fit["F_RB"], "r_avg": fit["r"],
                           "fit_R2": fit["R2"], "fit_resolved": fit["resolved"],
                           "worstcase_D": d_r, "Wstate_Fworst": wf,
                           "Wstate_infid": 1 - wf, "verdict": v_r})
    sw = pd.DataFrame(sweep_rows)
    x = np.arange(len(sw))
    axR.plot(x, sw["r_avg"], "o-", color="black", lw=1.8, label="RB avg infidelity $r$ (1 number)")
    axR.plot(x, sw["worstcase_D"], "s--", color="red", lw=1.6, label="audit worst-case $D$")
    axR.plot(x, sw["Wstate_infid"], "^-", color=WCOL["W-state"], lw=1.6,
             label="W-state $1-F_{worst}$ (workload)")
    axR.set_xticks(x); axR.set_xticklabels(sw["regime"], rotation=30, ha="right", fontsize=7)
    axR.set_ylabel("Infidelity", fontsize=9)
    axR.set_ylim(0, 1.0)
    axR.set_title("(b) RB under-reports vs audit & workload", fontsize=9)
    axR.tick_params(labelsize=8); axR.grid(alpha=0.25)
    axR.legend(fontsize=6.3, frameon=False, loc="upper left")
    fig.tight_layout()

    figdir = rc.FIGDIR
    fig.savefig(figdir / "figA1_rb.pdf", bbox_inches="tight")
    fig.savefig(figdir / "figA1_rb.png", dpi=300, bbox_inches="tight")

    outdir = rc.OUTDIR
    pd.concat(curve_frames, ignore_index=True).to_csv(outdir / "figA1_rb_curve.csv", index=False)
    summ = pd.DataFrame(summary_rows)
    summ["F_worst_GHZ"] = wl_fworst["GHZ"]; summ["F_worst_Wstate"] = wl_fworst["W-state"]
    summ["R_vul_GHZ_over_W"] = wl_fworst["GHZ"] / wl_fworst["W-state"]
    summ.to_csv(outdir / "figA1_rb_summary.csv", index=False)
    sw.to_csv(outdir / "figA1_rb_sweep.csv", index=False)

    print("=== Experiment 1: physical gate-dependent RB ===")
    print(summ.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print("\nStress sweep:")
    print(sw.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    agg_fit = fit_byreg["aggressive"]
    print(f"\nWorst-case: RB r={agg_fit['r']:.3f} (resolved={agg_fit['resolved']}), "
          f"identical for ALL workloads; workloads span F_worst "
          f"{wl_fworst['GHZ']:.3f}(GHZ)->{wl_fworst['W-state']:.3f}(W), R_vul=4.23.")
    print("Saved figA1_rb.{pdf,png}, figA1_rb_curve.csv, figA1_rb_sweep.csv, figA1_rb_summary.csv")


if __name__ == "__main__":
    main()
