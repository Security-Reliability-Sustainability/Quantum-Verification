#!/usr/bin/env python3
"""
Regenerate the paper figures (A1, A1-IRB, A2, A3, A4) from the experiment
CSVs, with uniform publication styling and the terminology used in the paper.

Usage:
    python experiments/make_paper_figures.py
Reads CSVs from experiments/ and writes PDFs+PNGs to experiments/figures_paper/.
Override with the environment variables DATA_DIR and OUT_DIR.
"""
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

_HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.environ.get("DATA_DIR", _HERE)
OUT_DIR = os.environ.get("OUT_DIR", os.path.join(_HERE, "figures_paper"))
os.makedirs(OUT_DIR, exist_ok=True)

mpl.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.titlesize": 10,
    "axes.labelsize": 10,
    "legend.fontsize": 8.5,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linewidth": 0.5,
    "axes.axisbelow": True,
    "axes.unicode_minus": False,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})
# Okabe-Ito colourblind-safe palette
C = dict(blue="#0072B2", orange="#E69F00", green="#009E73", verm="#D55E00",
         purple="#CC79A7", sky="#56B4E9", yellow="#F0E442", black="#000000",
         grey="#7F7F7F")

def p(name):
    return os.path.join(DATA_DIR, name)

def save(fig, stem):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT_DIR, f"{stem}.{ext}"))
    plt.close(fig)
    print("wrote", stem)

# ===========================================================================
# Fig A1 : single-qubit gate-dependent RB on the stressed model
# ===========================================================================
def fig_a1():
    cur = pd.read_csv(p("figA1_rb_curve.csv"))
    sw = pd.read_csv(p("figA1_rb_sweep.csv"))

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(9.2, 3.5))

    style = {"calibrated": (C["green"], "Calibrated  $F_{RB}=0.9999$"),
             "moderate":   (C["orange"], "Moderate stress  $F_{RB}=0.786$"),
             "aggressive": (C["verm"], "Worst-Case stress  $F_{RB}=0.780$")}
    for reg, (col, lab) in style.items():
        d = cur[cur.Regime == reg]
        axL.errorbar(d.m, d.survival_mean, yerr=d.survival_std, marker="o", ms=4,
                     lw=1.6, color=col, capsize=2, label=lab)
    axL.axhline(0.5, ls=":", color=C["grey"], lw=1.0)
    axL.text(120, 0.515, "$1/d$ floor ($d{=}2$)", color=C["grey"], ha="right", fontsize=8)
    axL.set_xscale("log", base=2)
    axL.set_xlabel("Sequence length $m$ (Cliffords)")
    axL.set_ylabel("Survival probability")
    axL.set_title("(a) Randomized benchmarking on stressed $H(t)$")
    axL.legend(loc="upper right", framealpha=0.95)
    axL.set_ylim(0.46, 1.02)

    x = np.arange(len(sw))
    axR.plot(x, sw.Wstate_infid, marker="s", lw=1.8, color=C["verm"],
             label=r"W-state  $1-F_{worst}$")
    axR.plot(x, sw.worstcase_D, marker="^", lw=1.8, color=C["blue"],
             label=r"process infidelity  $D$")
    axR.plot(x, sw.r_avg, marker="o", lw=1.8, color=C["black"],
             label=r"RB error rate  $r$")
    axR.axhline(0.30, ls="--", color=C["grey"], lw=1.0)
    axR.text(0.05, 0.315, "Hard Fail ($D=0.30$)", color=C["grey"], fontsize=8)
    axR.set_xticks(x)
    axR.set_xticklabels(["Calib", "0.1", "0.2", "0.3", "0.4", "Worst"])
    axR.set_xlabel("Stress level ($J=A$)")
    axR.set_ylabel("Value")
    axR.set_title("(b) RB error $r$ vs $D$ and workload $1-F_{worst}$")
    axR.legend(loc="center left", framealpha=0.95)
    axR.set_ylim(-0.03, 1.0)

    fig.tight_layout()
    save(fig, "figA1_rb")

# ===========================================================================
# Fig A1-IRB : interleaved RB of the echoed CR-CNOT + per-gate extrapolation
# ===========================================================================
def fig_a1_irb():
    cur = pd.read_csv(p("figA1_irb_curve.csv"))

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(9.2, 3.5))

    decay = cur[cur.m > 0]
    for reg, col, lab in [("moderate", C["orange"], "Moderate  $F_{CR}=0.486$ (interleaved)"),
                          ("aggressive", C["blue"], "Worst-Case  $F_{CR}=0.564$ (interleaved)")]:
        d = decay[decay.regime == reg]
        axL.plot(d.m, d.survival_int, marker="o", ms=4, lw=1.7, color=col, label=lab)
        axL.plot(d.m, d.survival_ref, marker="s", ms=3, lw=1.0, ls=":", color=col,
                 alpha=0.8, label=f"{('Moderate' if reg=='moderate' else 'Worst-Case')} reference ($p=1$)")
    axL.axhline(0.25, ls=":", color=C["grey"], lw=1.0)
    axL.text(28, 0.265, "$1/d$ floor ($d{=}4$)", color=C["grey"], ha="right", fontsize=8)
    axL.set_xscale("log", base=2)
    axL.set_xlabel("Sequence length $m$ (Cliffords)")
    axL.set_ylabel("Survival probability")
    axL.set_title("(a) Interleaved RB of the echoed CR-CNOT")
    axL.legend(loc="upper right", framealpha=0.95, fontsize=7.5)
    axL.set_ylim(0.2, 1.04)

    wl = cur[cur.m == -1].copy()
    order = ["GHZ", "W-state", "QFT", "VQE-SU2", "VQE-TwoLocal"]
    wl["workload"] = pd.Categorical(wl.workload, categories=order, ordered=True)
    wl = wl.sort_values("workload")
    xx = np.arange(len(wl)); w = 0.27
    axR.bar(xx - w, wl.RB_pred_fid, w, color="#444444", label=r"RB prediction  $F_{CR}^{\,n_{cx}}$")
    axR.bar(xx, wl.F_expected, w, color=C["sky"], label=r"expected  $F_{exp}$")
    axR.bar(xx + w, wl.F_worst, w, color=C["verm"], label=r"worst-case  $F_{worst}$")
    axR.set_xticks(xx)
    axR.set_xticklabels([f"{r.workload}\n($n_{{cx}}$={int(r.n_cx)})" for _, r in wl.iterrows()],
                        fontsize=7.5)
    axR.set_ylabel("Output-state fidelity")
    axR.set_title("(b) Per-gate RB prediction vs workload audit")
    axR.legend(loc="upper right", framealpha=0.95, fontsize=7.5)
    axR.set_ylim(0, 0.62)
    fig.tight_layout()
    save(fig, "figA1_irb")

# ===========================================================================
# Fig A2 : two-qubit shadow-assisted output-state fidelity (Path B)
# ===========================================================================
def fig_a2():
    cur = pd.read_csv(p("figA2_shadow_curve.csv"))
    fig, ax = plt.subplots(figsize=(5.4, 3.6))
    for key, col, lab in [("w_state", C["blue"], "W-state  exact $F_{worst}=0.093$"),
                          ("ghz", C["green"], "GHZ  exact $F_{worst}=0.393$")]:
        d = cur[cur.Workload == key]
        ax.errorbar(d.shots, d.F_shadow_mean, yerr=d.F_shadow_std, marker="o", ms=4,
                    lw=1.5, color=col, capsize=2)
        ax.axhline(d.F_exact.iloc[0], ls="--", lw=1.2, color=col, alpha=0.8)
    ax.axvline(2000, ls=":", color=C["grey"], lw=1.2)
    ax.text(2150, 0.62, "2000 shots\n($|\\Delta F|<0.03$)", color=C["grey"], fontsize=7.5)
    ax.axhspan(0.5, 1.0, color="0.92", zorder=0)
    ax.text(115, 0.93, "utility reference band ($F>0.5$)", color=C["grey"], fontsize=7.5)
    ax.set_xscale("log")
    ax.set_xlabel("Total shots")
    ax.set_ylabel("Estimated output-state fidelity  $\\hat{F}$")
    ax.set_title("Shadow-assisted output-state fidelity (Path B, Worst-Case stress)")
    ax.set_ylim(0, 1.0)
    custom = [Line2D([0],[0], color=C["blue"], marker="o", lw=1.5),
              Line2D([0],[0], color=C["green"], marker="o", lw=1.5),
              Line2D([0],[0], color="k", ls="--", lw=1.2)]
    ax.legend(custom, ["W-state  exact $F_{worst}=0.093$", "GHZ  exact $F_{worst}=0.393$",
                       "exact value (dashed)"], loc="center right", framealpha=0.95)
    fig.tight_layout()
    save(fig, "figA2_shadow")

# ===========================================================================
# Fig A3 : (J, A) sensitivity heatmaps  (W-state F_worst ; process infidelity D)
# ===========================================================================
def fig_a3():
    sw = pd.read_csv(p("figA3_sweep.csv"))
    Js = sorted(sw.J.unique()); As = sorted(sw.A.unique())
    Dm = sw.pivot(index="A", columns="J", values="D").reindex(index=As, columns=Js)
    Fm = sw.pivot(index="A", columns="J", values="Wstate_Fworst").reindex(index=As, columns=Js)
    ext = [min(Js)-0.05, max(Js)+0.05, min(As)-0.05, max(As)+0.05]

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(9.2, 3.7))
    im1 = axL.imshow(Fm.values, origin="lower", extent=ext, aspect="auto",
                     cmap="viridis", vmin=0, vmax=1)
    axL.set_title("(a) W-state worst-case fidelity $F_{worst}$")
    fig.colorbar(im1, ax=axL, fraction=0.046, pad=0.04)
    im2 = axR.imshow(Dm.values, origin="lower", extent=ext, aspect="auto",
                     cmap="magma_r", vmin=0, vmax=0.35)
    axR.set_title("(b) Process infidelity $D$")
    fig.colorbar(im2, ax=axR, fraction=0.046, pad=0.04)
    for ax in (axL, axR):
        ax.set_xlabel("Coupling $J$"); ax.set_ylabel("Drive $A$")
        ax.scatter([0.5], [0.5], marker="*", s=180, c="white",
                   edgecolors="k", linewidths=0.8, zorder=5)
    axR.contour(Dm.columns, Dm.index, Dm.values, levels=[0.30],
                colors=["#00FF7F"], linewidths=1.5)
    fig.suptitle("Joint $(J,A)$ stress sweep: Worst-Case corner $(0.5,0.5)\\star$ maximizes $D$; "
                 "region $J,A\\geq0.4$ fails broadly", fontsize=9.5)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save(fig, "figA3_sensitivity")

# ===========================================================================
# Fig A4 : process infidelity D  vs  worst-case output-state fidelity F_worst
# ===========================================================================
def fig_a4():
    pts = pd.read_csv(p("figA4_points.csv"))
    sw = pd.read_csv(p("figA3_sweep.csv"))
    s = sw.sort_values("D")
    r_pearson = np.corrcoef(sw.D, sw.Wstate_Fworst)[0, 1]
    from scipy.stats import spearmanr
    r_spear = spearmanr(sw.D, sw.Wstate_Fworst).statistic
    print(f"A4 sweep Pearson={r_pearson:.3f} Spearman={r_spear:.3f}")

    fig, ax = plt.subplots(figsize=(6.0, 3.9))
    ax.axvspan(0, 0.02, color=C["green"], alpha=0.12)
    ax.axvspan(0.30, 0.42, color=C["verm"], alpha=0.10)
    ax.text(0.005, 0.05, "Pass\n($D\\leq0.02$)", color=C["green"], fontsize=7.5)
    ax.text(0.305, 0.86, "Hard Fail\n($D\\geq0.30$)", color=C["verm"], fontsize=7.5)
    ax.plot(s.D, s.Wstate_Fworst, color=C["grey"], lw=1.2, alpha=0.7, zorder=1,
            label="W-state $(J,A)$ sweep")
    wcol = {"GHZ": C["blue"], "W-state": C["verm"], "QFT": C["green"],
            "VQE-SU2": C["purple"], "VQE-TwoLocal": C["sky"]}
    for wl, col in wcol.items():
        d = pts[pts.workload == wl]
        ax.scatter(d.channel_D, d.F_worst, s=46, color=col, edgecolors="k",
                   linewidths=0.5, zorder=3, label=wl)
    ax.set_xlabel("Process infidelity $D$ ($R_y$ audit)")
    ax.set_ylabel("Worst-case fidelity $F_{worst}$")
    ax.set_title("Process infidelity $D$ vs worst-case fidelity $F_{worst}$\n"
                 f"(Pearson r = {r_pearson:.2f}, Spearman rho = {r_spear:.2f})",
                 fontsize=9)
    ax.set_xlim(-0.02, 0.42); ax.set_ylim(0, 1.03)
    ax.legend(loc="upper right", ncol=2, framealpha=0.95, fontsize=7.5)
    fig.tight_layout()
    save(fig, "figA4_D_vs_F")

if __name__ == "__main__":
    fig_a1(); fig_a1_irb(); fig_a2(); fig_a3(); fig_a4()
    print("All figures written to", OUT_DIR)
