# INVENTORY — Paper concept → code location

This maps each concept in the paper to the canonical function in the existing repo.
The experiments in `experiments/` **import** these; they do not re-implement the physics.

Repo root: `qsys_audit/` (src-layout under `src/qsys_audit/`).
Global seed for all experiment scripts: **1234** (see `docs/ENV.txt`).

| Paper concept | Code location | Notes |
|---|---|---|
| 3-qubit time-dependent Hamiltonian `H(t)`, chain q0–Q–q2 | `src/qsys_audit/physics_engine.py :: PhysicsEngine._background_terms_3q` | Builds static + drive background. Gate term added per-call. |
| YX nearest-neighbor coupling `J1Q·σy⊗σx⊗I`, `JQ2·I⊗σx⊗σy` | `_background_terms_3q` (`h_couple_1q`, `h_couple_q2`) | Coupling is σy⊗σx (q1Q) and σx⊗σy (Qq2). |
| Neighbor stress drives `A1·σx(q0)`, `A2·σx(q2)`, readout `A_ro·σx(Q)` | `_background_terms_3q` (`h_drive1`, `h_drive2`, `h_ro`) | `A1,A2` use **chirp** envelope (`chirp_rate=0.5`), `A_ro` uses **cos** envelope at `delta_ro`. |
| Single-qubit drive (gate) term `A f(t) σy(Q)` | `PhysicsEngine.simulate_probe_state` (`h_gate_q`) | Gaussian envelope, area-calibrated: `a_q=(θ/2)/(σ√(2π))`. |
| QuTiP `mesolve`, unitary (`c_ops=[]`) | `simulate_probe_state` (probe), `evolve_echoed_cr_cnot` (CNOT), `scripts/run_shadow_convergence.py::get_final_state_for_gate` (gate suite) | All runs use `c_ops=[]` → unitary evolution. |
| Propagator / final state & partial trace to QUV channel `E_phys(ρ)=Tr_spec[U(ρ⊗ρ_spec)U†]` | `simulate_probe_state` → `sol.states[-1].ptrace(1)` | Spectators initialized to `|0⟩⟨0|` on q0 and q2; Q is index 1. |
| One-qubit QPT → 4×4 Choi | `reconstruction.py :: ReconstructionEngine.perform_qpt_pipeline` | 4-probe set {\|0⟩,\|1⟩,\|+⟩,\|+i⟩}; standard Choi assembly. |
| Ideal Choi for Ry(θ) | `ReconstructionEngine.ideal_choi_ry` | (gate-suite uses `run_shadow_convergence.ideal_choi_from_unitary`). |
| Process infidelity `D = 1 − (1/d²)Tr[C_ideal† C_phys]`, d=2 | `ReconstructionEngine.process_distance_d` and `.process_distance_d_physical` | `_physical` first CPTP-projects `C_phys` (PSD + TP alternating projection, `project_choi_to_cptp_1q`). This is the one used by the CLI/scripts. |
| Polar decomposition → `U_eff`, residual `N` | `fit_unitary_plus_noise_from_choi` (`scipy.linalg.polar` on dominant Kraus `k0`) | Returns `U_eff`, Pauli-channel probs of residual `N`. |
| `δ_log` (Fubini–Study logical deviation) | `ReconstructionEngine.delta_log` | `√(max(0,1−|Tr(U_ideal† U_eff)|²/4))`. |
| `δ_inc = 1 − p_I` (incoherent residual weight) | `ReconstructionEngine.delta_inc` | `p_I` = identity component of residual Pauli channel. |
| Verdict classifier (Pass / Soft / Hard Fail) | `src/qsys_audit/cli.py :: classify_verdict` (duplicated in `run_shadow_convergence.py`) | Thresholds from `config.yaml::verification_thresholds`. |
| Single-qubit classical-shadow estimator (shots, K MoM blocks, CPTP) | `reconstruction.py :: single_qubit_pauli_shadow_samples`, `estimate_rho_from_shadows`, `median_of_means`, `project_density_matrix` | Per shot: random Pauli basis ∈ {X,Y,Z}, sample bit, `ρ̂ = 3·P_b − I`. MoM over K blocks, then PSD/trace projection. |
| Two-qubit workload execution path (native {Rx,Ry,Rz,CNOT}, echoed CR-CNOT, probe inputs) | `workload_replay.py :: WorkloadReplay.simulate_2q_workload_physical` / `.simulate_2q_workload_ideal` / `.benchmark_workload` | Probe inputs {\|00⟩,\|10⟩,\|+0⟩,\|0+⟩}. CNOT via `evolve_echoed_cr_cnot`. |
| `F_exp / F_worst / F_spread` | `WorkloadReplay.benchmark_workload` → `WorkloadMetrics` | mean / min / (max−min) over the 4 probe inputs. |
| Echoed CR-CNOT pulse (θ_ZX, Rx/Rz corrections, β_IX, γ_ZI) | `PhysicsEngine.evolve_echoed_cr_cnot` | Two Gaussian CR halves (+/−), π echo on Q, then Rx(−π/2) on q2 + Rz(−1.6408) on Q. |
| Three stress regimes | `config.yaml::stress_profiles` (`calibrated`, `moderate`, `aggressive`) | See spec match below. |
| Thresholds | `config.yaml::verification_thresholds` | See spec match below. |

## Stress regime spec match (config.yaml)
- **calibrated**: `J1Q=JQ2=0`, `A1=A2=Aro=0`, `delta=0` ✅ (handoff "Calibrated baseline")
- **moderate**: `J1Q=JQ2=0.5`, `delta1=delta2=10`, `A1=A2=Aro=0` ✅ (handoff "Moderate / passive crosstalk")
- **aggressive**: `J1Q=JQ2=0.5`, `A1=A2=Aro=0.5`, `delta1=delta2=10` ✅ (handoff "Worst-case / bounded active stress")

> Note: the handoff calls the third regime "worst-case"; the code alias is `aggressive`. They are the same profile.

## Threshold spec match (config.yaml)
- Pass: `D ≤ 0.02` AND `δ_log ≤ 0.05` AND `δ_inc ≤ 0.05` ✅
- Hard Fail: `D ≥ 0.30` OR `δ_log ≥ 0.90` OR `δ_inc ≥ 0.50` ✅
- Soft Fail: otherwise ✅

## Unit convention (IMPORTANT — used for the Exp-3 physical-units mapping)
From the docstring of `PhysicsEngine.simulate_probe_state` (lines ~183–196):

> The simulator works in the **rotating frame** w.r.t. the bare qubit frequencies, so there is
> **no** explicit `ω·σz/2` drift term — a frame choice, not a missing effect. All Hamiltonian
> rates — static couplings `J`, drive amplitudes `A`, detunings `δ` — are **dimensionless,
> expressed in natural units of `1/T_gate`**, where `T_gate = config.single_qubit_gate.duration`.

So `J=0.5`, `A=0.5`, `δ=10` are **rates in units of 1/T_gate**. `config.single_qubit_gate.duration = 1.0`
(dimensionless). To convert to physical frequency the user must fix a concrete `T_gate`.

**Conversion rule used in Exp 3:** a dimensionless rate `g` maps to angular frequency `g/T_gate`
(rad/time), i.e. ordinary frequency `f = g / (2π·T_gate)`. Anchoring the single-qubit gate to a
representative superconducting value `T_gate = 20 ns` gives `1/T_gate = 50 MHz` (as a rate), and
`f = g/(2π·T_gate)`:
- `J = 0.5` → `0.5/(2π·20ns) ≈ 4.0 MHz` coupling.
- `A = 0.5` → `≈ 4.0 MHz` drive amplitude.
- `δ = 10` → `10/(2π·20ns) ≈ 80 MHz` detuning.

These are flagged in `RESULTS.md` for the author to sanity-check against the real device anchor
they prefer (the conversion scales linearly with the chosen `T_gate`).

### Paper §/Eq. cross-references and provenance of the 20 ns anchor (verified against the PDF)
- **Eq. (1):** `H(t) = H_Q(t) + H_rest(t)` (total = ideal victim + bounded stress).
- **Eq. (2):** single-qubit drive term; defines `ω_d^(i) = ω_i + δ_i`, i.e. **δ is a drive detuning**
  (consistent with the code, which uses δ as the stress-drive envelope frequency).
- **Eq. (3):** YX nearest-neighbor coupling `H_YX = Σ_<j,k> (J_jk/2)(σ_y^j σ_x^k + σ_x^j σ_y^k)`.
- **Section 4 methodology, "Step 2A"** sets the worst-case regime `J=0.5, A=0.5, δ=10.0`
  ("solve the time-dependent Schrödinger equation"), cited to refs [2,3,14] (and [5,15,16]).
- **The paper contains NO physical units** (0 occurrences of GHz/MHz/kHz); all stress values are
  bare dimensionless numbers. The code likewise hard-codes no physical time — only the
  dimensionless `single_qubit_gate.duration = 1.0` and `cnot.duration = 1.5`.
- **`T_gate = 20 ns` is therefore an EXTERNAL anchor introduced only for the physical-units line**,
  not present in the paper or code. Chosen as a representative superconducting single-qubit gate
  time; CR-CNOT = 1.5·T_gate = 30 ns at this anchor. All MHz values scale linearly with T_gate.
- **Prefactor caveat:** Eq. (3) carries a `J/2` prefactor (Eq. (2) an `Ω/2`). The numbers above
  convert the *labeled parameter* `J=0.5 → 3.98 MHz`. If instead quoting the *physical σ_yσ_x
  coefficient in H* (= J/2 = 0.25), that is `≈ 1.99 MHz`. Be consistent about which is cited.
