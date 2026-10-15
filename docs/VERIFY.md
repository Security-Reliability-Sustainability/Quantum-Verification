# VERIFY — reproduction of Section-3 ground-truth numbers

Re-ran the **existing** pipeline (no physics re-implementation) with global seed 1234.
All published numbers reproduce within rounding / Monte-Carlo noise. ✅ **Gate passed — new experiments are trustworthy.**

## Ry(π/2) single-qubit channel (`perform_qpt_pipeline` → `process_distance_d_physical`)

| Regime | Published | Reproduced | Verdict |
|---|---|---|---|
| calibrated | D ≈ 1e−4, PASS | D = 0.0001, δ_log = 0.0098, δ_inc = 0.0000 | PASS ✅ |
| moderate | D = 0.269, δ_log = 0.11, δ_inc = 0.26, SOFT FAIL | D = 0.2687, δ_log = 0.1065, δ_inc = 0.2605 | SOFT FAIL ✅ |
| worst-case (aggressive) | D = 0.338, δ_log = 0.39, HARD FAIL | D = 0.3377, δ_log = 0.3907, δ_inc = 0.2618 | HARD FAIL ✅ |

## Workload worst-case fidelity, aggressive regime (`benchmark_workload`)

| Workload | Published F_worst | Reproduced F_worst |
|---|---|---|
| GHZ | 0.393 | 0.3930 ✅ |
| W-state | 0.093 | 0.0928 ✅ |
| QFT | 0.125 | 0.1251 ✅ |
| VQE-SU2 | 0.174 | 0.1737 ✅ |
| VQE-TwoLocal | 0.102 | 0.1016 ✅ |

Vulnerability ratio R_vul = F_worst(GHZ)/F_worst(W) = 0.3930/0.0928 = **4.234** (published 4.23 ✅).

## T-gate multi-basis masking (`plot_tgate_masking` path)
Worst-basis concealed gap `Gap_min` (= max over MUB of 1−F; the {|0⟩,|1⟩,|+i⟩,|−i⟩} bases
all show this gap while {|+⟩,|−⟩} conceal it at gap 0):

| Regime | Published | Reproduced |
|---|---|---|
| moderate | 0.354 (~35%) | 0.3540 ✅ |
| worst-case | 0.416 (~42%) | 0.4159 ✅ |

## Single-qubit shadow vs exact (`run_shadow_convergence` path)
Ry(π/2), aggressive, 5000 shots, K=10, seed 1234:
- D_exact = 0.3377, D_shadow = 0.3309, **|ΔD| = 0.0068 < 0.03** ✅
- Verdicts agree: HARD FAIL / HARD FAIL ✅

## Conclusion
Every Section-3 ground-truth number reproduces. The build (qutip 5.2.3 from source, qiskit 2.4.2,
mqt.bench 2.2.2) is numerically consistent with the submitted paper. New experiments proceed.
