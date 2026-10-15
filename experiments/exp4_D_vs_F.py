"""Experiment 4: channel-D versus workload-fidelity correlation.

Answers whether the channel-level D threshold connects to workload output
fidelity. Combines:
  (a) the W-state J-A sweep (continuous D vs F_worst) from exp3 (figA3_sweep.csv), and
  (b) all 5 workloads x 3 regimes: regime channel D (Ry audit) vs workload F_worst.

Shades Pass (D<=0.02) and Hard-Fail (D>=0.30) bands, annotates the corresponding
F_worst bands, reports Pearson/Spearman correlation, and overlays the coherent-
accumulation trend 1 - F_worst = O(n_cx^2 * eps^2) discussed in the paper.

Output: figs/figA4_D_vs_F.{pdf,png} + figA4_points.csv
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import qutip as qt
from scipy.stats import pearsonr, spearmanr

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import exp_common as rc

SEED = rc.GLOBAL_SEED
WORKLOADS = [("ghz", "GHZ"), ("w_state", "W-state"), ("qft", "QFT"),
             ("vqe_su2", "VQE-SU2"), ("vqe_two_local", "VQE-TwoLocal")]
REGIMES = ["calibrated", "moderate", "aggressive"]
WCOL = {"GHZ": "#0072B2", "W-state": "#D55E00", "QFT": "#009E73",
        "VQE-SU2": "#CC79A7", "VQE-TwoLocal": "#56B4E9"}


def main():
    np.random.seed(SEED)
    cfg, physics, recon, replay = rc.build_engines()
    theta = float(cfg["single_qubit_gate"]["theta_gate"])

    # channel D per regime (single-qubit Ry audit)
    Dreg = {}
    for reg in REGIMES:
        d, dl, di, v = rc.channel_D_and_verdict(recon, physics, cfg,
                                                physics.stress_profiles()[reg], theta)
        Dreg[reg] = (d, v)

    # workload F_worst per (workload, regime)
    rows = []
    for key, name in WORKLOADS:
        for reg in REGIMES:
            m = replay.benchmark_workload(key, physics.stress_profiles()[reg])
            rows.append({"workload": name, "regime": reg, "n_cx": m.n_cx,
                         "channel_D": Dreg[reg][0], "regime_verdict": Dreg[reg][1],
                         "F_worst": m.f_worst})
    df = pd.DataFrame(rows)

    # W-state continuous sweep (reuse exp3 output)
    sweep_path = rc.OUTDIR / "figA3_sweep.csv"
    sweep = pd.read_csv(sweep_path) if sweep_path.exists() else None

    outdir = rc.OUTDIR
    df.to_csv(outdir / "figA4_points.csv", index=False)

    # correlations
    pr, _ = pearsonr(df["channel_D"], df["F_worst"])
    sr, _ = spearmanr(df["channel_D"], df["F_worst"])
    pr_sw = sr_sw = float("nan")
    if sweep is not None:
        pr_sw, _ = pearsonr(sweep["D"], sweep["Wstate_Fworst"])
        sr_sw, _ = spearmanr(sweep["D"], sweep["Wstate_Fworst"])

    # coherent-accumulation collapse: 1 - F_worst vs n_cx^2 * D
    df["infid"] = 1.0 - df["F_worst"]
    df["ncx2_D"] = (df["n_cx"] ** 2) * df["channel_D"]
    sub = df[df["channel_D"] > 0.05]
    pr_coh, _ = pearsonr(sub["ncx2_D"], sub["infid"]) if len(sub) > 2 else (float("nan"), 0)

    # ---- figure ----
    fig, ax = plt.subplots(figsize=(5.4, 3.6))
    if sweep is not None:
        ss = sweep.sort_values("D")
        ax.plot(ss["D"], ss["Wstate_Fworst"], "-", color="grey", lw=1.3, alpha=0.8,
                label="W-state J–A sweep")
    for name in [n for _, n in WORKLOADS]:
        d = df[df.workload == name]
        ax.scatter(d["channel_D"], d["F_worst"], s=48, color=WCOL[name],
                   edgecolor="black", lw=0.4, zorder=4, label=name)

    ax.axvspan(0.0, 0.02, color="#2ca02c", alpha=0.12)
    ax.axvspan(0.30, max(0.42, df["channel_D"].max() + 0.02), color="#d62728", alpha=0.10)
    ax.text(0.001, 0.05, "Pass\nD≤0.02", fontsize=6, color="#2ca02c")
    ax.text(0.305, 0.78, "Hard Fail\nD≥0.30", fontsize=6, color="#d62728")

    ax.set_xlabel("Channel process infidelity $D$ (Ry audit)", fontsize=9)
    ax.set_ylabel("Workload $F_{worst}$", fontsize=9)
    ax.set_title(f"D ↔ workload fidelity (Pearson r={pr:.2f}, sweep r={pr_sw:.2f})", fontsize=9)
    ax.set_xlim(-0.01, max(0.42, df["channel_D"].max() + 0.02))
    ax.set_ylim(0.0, 1.02)
    ax.tick_params(labelsize=8)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=6.5, frameon=False, ncol=2, loc="upper right")
    fig.tight_layout()

    figdir = outdir / "figs"
    fig.savefig(figdir / "figA4_D_vs_F.pdf", bbox_inches="tight")
    fig.savefig(figdir / "figA4_D_vs_F.png", dpi=300, bbox_inches="tight")

    print("=== Experiment 4: channel-D vs workload-F correlation ===")
    print(df.pivot(index="workload", columns="regime", values="F_worst").to_string(
        float_format=lambda x: f"{x:.4f}"))
    print(f"\nchannel D by regime: " + ", ".join(f"{r}={Dreg[r][0]:.4f}({Dreg[r][1]})" for r in REGIMES))
    print(f"Pearson r (pooled 15 pts) = {pr:.3f}; Spearman = {sr:.3f}")
    print(f"Pearson r (W-state sweep) = {pr_sw:.3f}; Spearman = {sr_sw:.3f}")
    print(f"Coherent-accumulation collapse: corr(1-F_worst , n_cx^2*D) = {pr_coh:.3f}")
    print("F_worst bands: D<=0.02 -> F_worst ~ %.3f (pass); D>=0.30 -> F_worst in [%.3f, %.3f] (hard fail)" % (
        df[df.channel_D <= 0.02].F_worst.min() if (df.channel_D<=0.02).any() else float('nan'),
        df[df.channel_D >= 0.30].F_worst.min(), df[df.channel_D >= 0.30].F_worst.max()))
    print("Saved figA4_D_vs_F.{pdf,png}, figA4_points.csv")


if __name__ == "__main__":
    main()
