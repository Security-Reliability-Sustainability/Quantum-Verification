# RESULTS - experiment numbers

Seed 1234. All scripts import the existing pipeline (no physics re-implementation).
Section-3 ground-truth reproduced first (see `VERIFY.md`). Figures in `figs/`, raw numbers in the
adjacent CSVs. Nothing below **contradicts** the submitted paper (see "Consistency" at the end).

---

## Exp 1 - Randomized benchmarking on the stressed channel (Fig A1)
`exp1_rb.py` → `figA1_rb.{pdf,png}`, `figA1_rb_curve.csv`, `figA1_rb_sweep.csv`, `figA1_rb_summary.csv`

**Physical, gate-dependent RB** (single seed 1234, N_seq = 40, m ∈ {1…128}). Each of the 24
single-qubit Cliffords is compiled into a word over the physical primitives {X90, Y90, X-90, Y-90};
every primitive is simulated as a real pulse through the paper's stressed H(t) and assembled into a
4×4 superoperator. Because the anisotropic stress makes X90 ≠ Y90 error, the noise is genuinely
**gate-dependent** — so the decay is **not** an analytically-forced exponential and R² (plus a
resolved/unresolved flag) is a real measurement, not a modeling artifact.

| Regime | F_RB | RB avg r/Clifford | fit R² | fit resolved? | survival(m=1) | audit worst-case D | verdict |
|---|---|---|---|---|---|---|---|
| Calibrated | 0.9999 | 0.0001 | 0.994 | yes | 1.000 | 0.0001 | PASS |
| Moderate | 0.787 | 0.214 | 0.998 | **no (floor)** | 0.559 | 0.269 | SOFT FAIL |
| Worst-case | 0.780 | 0.220 | 0.989 | **no (floor)** | 0.530 | 0.338 | HARD FAIL |

Stress sweep (RB stays low while audit D and workload infidelity climb):

| Regime | RB r | audit D | W-state 1−F_worst | verdict |
|---|---|---|---|---|
| calib | 0.000 | 0.000 | 0.000 | PASS |
| J=A=0.1 | 0.040 | 0.021 | 0.093 | SOFT FAIL |
| J=A=0.2 | 0.129 | 0.077 | 0.324 | SOFT FAIL |
| J=A=0.3 | 0.165 | 0.159 | 0.593 | SOFT FAIL |
| J=A=0.4 | 0.160 | 0.250 | 0.802 | SOFT FAIL |
| Worst-case | 0.220 | 0.338 | 0.907 | HARD FAIL |

**Honest reading (not benign — three structural blind spots):** under this strong-coupling model the
single-qubit gates are genuinely degraded, so RB *does* register error — we do **not** claim it looks
clean. The QCVV gap is structural: **(i) under-reporting** — RB's average r (≤ 0.22) is far below the
audit worst-case D (0.34) and the coherent δ_log (0.39), and dramatically below the workload
infidelity it should predict (W-state 1−F_worst = 0.91); at J=A=0.4 RB reads r = 0.16 (a merely
"mediocre" gate) while the W-state is already 0.80 corrupted. **(ii) ambiguity** — at moderate/worst
stress the survival saturates to the depolarized 1/d floor by m ≈ 2, so the single-exponential RB fit
is ill-conditioned (flagged `resolved = False`): RB cannot return a well-defined decay constant.
**(iii) workload-agnostic** — RB's one number is identical for GHZ (F_worst 0.39) and the Hard-Fail
W-state (0.093), so an RB-thresholded scheduler issues the same verdict to a near-pass tenant and a
catastrophically-exposed victim, blind to the 4.23× gap.

**Phrasing:** *"Gate-dependent RB on the stressed model returns a single average infidelity that both
under-reports and cannot localize the audited risk: at worst-case stress its decay saturates to the
1/d floor by m ≈ 2 (no well-defined r; Fig. A1a), and across the J–A sweep RB's r stays ≤ 0.22 while
the same stress drives the W-state to 1−F_worst = 0.91 and the audit D to 0.34 (Fig. A1b). RB also
reports one workload-agnostic number — identical for GHZ (F_worst 0.39) and the Hard-Fail W-state
(0.093) — so it cannot flag the 4.23× victim-specific exposure the cross-layer audit catches."*

---

## Exp 1b - Interleaved RB of the echoed CR-CNOT (Fig A1-IRB)
`exp1b_irb_cnot.py` → `figA1_irb.{pdf,png}`, `figA1_irb_curve.csv`, `figA1_irb_summary.csv`

Mirrors the workload-execution convention (single-qubit gates ideal, CR-CNOT physical). Noisy
CR-CNOT built as a 16×16 superoperator from 2-qubit process reconstruction; interleaved RB with a
random 24×24 local-Clifford twirl each step (genuine twirled decay, reference noiseless ⇒ p_ref = 1).

| Regime | p_int | fit R² (int) | r_CR-CNOT (IRB) | F_CR-CNOT (IRB) |
|---|---|---|---|---|
| Moderate | 0.315 | 0.999 | 0.514 | 0.486 |
| Worst-case | 0.419 | 1.000 | 0.436 | 0.564 |

The CR-CNOT is the dominant failure source (average gate fidelity ≈ 0.49–0.56 — a clearly broken
entangler, consistent with HARD FAIL). (Mild non-monotonicity — worst-case F_CR slightly above
moderate — comes from partial echo refocusing under the active drive, the same effect seen in
Exp 3 where W-state F_worst bottoms at A = 0.4.)

**RB per-CNOT extrapolation vs the workload audit (worst-case, F_CR = 0.564):**

| Workload | n_cx | RB pred = F_CRⁿ | actual F_expected | actual F_worst |
|---|---|---|---|---|
| GHZ | 1 | 0.564 | 0.408 | 0.393 |
| W-state | 2 | 0.318 | 0.304 | **0.093** |
| QFT | 5 | 0.057 | 0.256 | 0.125 |
| VQE-SU2 | 3 | 0.179 | 0.321 | 0.174 |
| VQE-TwoLocal | 3 | 0.179 | 0.309 | 0.102 |

The naive RB per-gate extrapolation mis-estimates the audited worst case in **both** directions:
it **over-predicts** fidelity for the low-depth entanglers (W-state 0.318 vs actual 0.093 — a 3.4×
optimism gap; GHZ 0.564 vs 0.393) and **under-predicts** for deep circuits (QFT 0.057 vs 0.125). So a
per-gate RB number — even a correctly measured one — is not a reliable proxy for workload worst-case
correctness; the cross-layer audit is required.

**Phrasing:** *"Interleaved RB measures the echoed CR-CNOT at average gate fidelity F_CR ≈ 0.56
(worst-case; clean twirled decay, R² = 1.0, Fig. A1-IRB-a), confirming the entangler as the dominant
failure mode. But extrapolating this single per-gate number across circuit depth mis-predicts the
audited worst case in both directions — over-optimistic for the W-state (0.32 predicted vs 0.093
measured, 3.4×) and over-pessimistic for QFT — so RB per-gate fidelity cannot substitute for the
cross-layer worst-case workload audit (Fig. A1-IRB-b)."*

---

## Exp 2 - Two-qubit classical-shadow fidelity (Fig A2)
`exp2_shadow2q.py` → `figA2_shadow.{pdf,png}`, `figA2_shadow_curve.csv`, `figA2_shadow_summary.csv`

Path-B extension: executed 2-qubit output state from the existing workload path; output fidelity
estimated by randomized single-qubit Pauli-basis shadows (same conventions as `reconstruction.py`,
median-of-means). 30 trials/point.

| Workload | exact F_worst | min shots: verdict preserved | min shots: reliable \|ΔF\|<0.03 (per-run std<0.03) | F̂ at that point | \|ΔF\| |
|---|---|---|---|---|---|
| **W-state** | 0.093 | 100 | **2000** | 0.0924 | 0.0005 |
| GHZ | 0.393 | 200 | 1000 | 0.3875 | 0.0055 |

The shadow estimator is **unbiased even at 100 shots** (mean \|ΔF\| < 0.015 throughout); the binding
budget is per-run variance.

**Phrasing:** *"Two-qubit Pauli-shadow estimation recovers the W-state worst-case fidelity
(exact F = 0.093) to within \|ΔF\| < 0.03 using ≈2000 shots, and preserves the HARD-FAIL verdict
(F ≪ 0.5) from as few as 100 shots (Fig. A2) — closing the single-qubit-only gap in Path B."*

---

## Exp 3 - J-A sensitivity sweep and physical-units mapping (Fig A3)
`exp3_sensitivity.py` → `figA3_sensitivity.{pdf,png}`, `figA3_sweep.csv`, `figA3_units.csv`, `figA3_lindblad.csv`

Grid J = J1Q = JQ2 and A = A1 = A2 = Aro ∈ {0,…,0.5}, δ = 10. Diagonal (J = A):

| J=A | channel D | W-state F_worst | verdict |
|---|---|---|---|
| 0.0 | 0.0001 | 1.000 | PASS |
| 0.1 | 0.021 | 0.907 | SOFT FAIL |
| 0.2 | 0.077 | 0.676 | SOFT FAIL |
| 0.3 | 0.159 | 0.407 | SOFT FAIL |
| 0.4 | 0.250 | 0.198 | SOFT FAIL |
| 0.5 | 0.338 | 0.093 | HARD FAIL |

- **Smooth & monotonic** (D non-decreasing, F_worst non-increasing) — no discontinuity ⇒ the
  verdict is not an artifact of one hand-picked corner.
- The (0.5,0.5) corner **maximizes channel D** (the audit metric driving the verdict). The W-state
  F_worst minimum (0.050) lies adjacent at (J=0.5, A=0.4); the **entire J,A ≥ 0.4 region fails**, so
  the corner is a conservative envelope, not a singular cherry-pick.
- **T1/T2 Lindblad check** (T1=T2=100 µs, 20 ns gate): D shifts 0.3377 → 0.3378, verdict unchanged
  (HARD FAIL) ⇒ the flagged error is **coherent crosstalk, not relaxation**.

**Physical-units mapping** (rates are dimensionless in units 1/T_gate; f = g/(2π·T_gate); anchor
T_gate = 20 ns — *scales linearly with T_gate, please sanity-check the anchor*):

| Param | dimensionless | ≈ MHz @ 20 ns | reported hardware range (cited, not fabricated) | refs |
|---|---|---|---|---|
| J (static coupling) | 0.5 | ≈ 4.0 MHz | static/ZZ & exchange crosstalk, ~0.01–5 MHz | [2],[3],[14] |
| A (neighbor/RO drive) | 0.5 | ≈ 4.0 MHz | microwave drive-crosstalk amplitudes, ~0.1–several MHz | [15],[16] |
| δ (detuning) | 10.0 | ≈ 80 MHz | qubit–qubit detunings, ~50–300 MHz | [3],[16] |

**Phrasing:** *"Across the full J–A grid the verdict transitions smoothly (Fig. A3); (0.5,0.5)
maximizes channel D and the entire J,A ≥ 0.4 region is HARD/SOFT-FAIL, so our worst-case corner is a
conservative envelope rather than a tuned point. At a 20 ns gate anchor, J = A = 0.5 ≈ 4 MHz and
δ = 10 ≈ 80 MHz — the high-but-bounded end of reported transmon crosstalk [2,3,14,15,16]. Adding
T1=T2=100 µs leaves the verdict unchanged, confirming the dominant error is coherent."*

> ⚠️ Flags for author (verified against the PDF):
> - **20 ns is an external anchor, NOT in the paper or code.** The paper has zero GHz/MHz/kHz and
>   states the stress values as bare numbers (Eq. (1)–(3); worst-case set in §4 "Step 2A", cited
>   [2,3,14]); the code carries only dimensionless durations (single-qubit = 1.0, CR-CNOT = 1.5).
>   Present 20 ns as a representative assumption; every MHz value scales linearly with it.
> - **Rates are angular** (they sit in H·exp(−iHt), ħ=1), so f = g/(2π·T_gate); the arithmetic
>   J=0.5→3.98 MHz, δ=10→79.6 MHz is exact at 20 ns. δ is specifically a *drive detuning* (Eq. (2)).
> - **Prefactor:** Eq. (3) has a J/2 prefactor — the labeled J=0.5 → 3.98 MHz, but the physical
>   σ_yσ_x coefficient (J/2=0.25) → 1.99 MHz. Quote one convention consistently.
> - Hardware ranges are paraphrased from [2,3,14,15,16] as order-of-magnitude bounds — verify exact
>   figures against the papers before the PDF.

---

## Exp 4 - Channel-D versus workload-fidelity correlation (Fig A4)
`exp4_D_vs_F.py` → `figA4_D_vs_F.{pdf,png}`, `figA4_points.csv`

- **Pearson r = −0.97** (pooled, 5 workloads × 3 regimes); **r = −0.93** on the continuous W-state
  J–A sweep (Spearman −0.93). Strong monotonic D ↔ F_worst link.
- Threshold bands map cleanly: **D ≤ 0.02 ⇒ F_worst ≈ 1.0 (Pass)**; **D ≥ 0.30 ⇒ F_worst ∈
  [0.09, 0.39] (Hard Fail)**. At fixed regime D, the residual vertical spread is the workload-
  structure (n_cx) dependence — exactly the 4.23× gap a single channel number cannot resolve.
- Coherent-accumulation overlay: corr(1−F_worst, n_cx²·D) = 0.58 (positive, consistent with the
  paper's `1−F_worst = O(n_cx²·ε²)` discussion; partial because workloads differ beyond CX count).

**Phrasing:** *"Channel infidelity D and workload F_worst are strongly anti-correlated (Pearson
r = −0.97; r = −0.93 on the continuous sweep): D ≤ 0.02 maps to F_worst ≈ 1.0 and D ≥ 0.30 to
F_worst ≤ 0.39, so the channel-level cutoffs do bound workload output fidelity (Fig. A4). The
residual spread at fixed D is the n_cx-dependent workload structure the audit is designed to
expose."*

---

## Positioning
Prior crosstalk work constructs or analyzes attacks and information leakage. This work supplies the
cross-layer, verdict-emitting verification primitive: a pre-dispatch audit rather than an attack.

---

## Consistency with the paper
All Section-3 anchors reproduced exactly (see `docs/VERIFY.md`): Ry D = 0.0001/0.269/0.338; F_worst GHZ 0.393,
W 0.093, QFT 0.125, VQE-SU2 0.174, VQE-TwoLocal 0.102; R_vul 4.23; T-gate gap 0.354/0.416; shadow
\|ΔD\| = 0.007. **No new number contradicts the paper.** New facts that merely *add* detail:
- RB (new measurement): avg r ≈ 0.22 / F_RB ≈ 0.78 at worst-case, with the decay saturating to the
  1/d floor (fit unresolved). Consistent with HARD FAIL; it adds the QCVV-gap evidence, no conflict.
- The W-state F_worst grid minimum is 0.050 at (J=0.5, A=0.4), slightly below the (0.5,0.5)
  value 0.093 — does not contradict anything; the paper's worst-case **corner** is defined by the
  channel metric D, which *is* maximal at (0.5,0.5). Worth a one-line footnote if cited.
