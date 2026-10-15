# QVSh / qsys_audit

Offline cross-layer verification of concurrent crosstalk in shared quantum
processors. The package simulates a three-qubit line at the pulse level under
bounded concurrent stress, reconstructs the single-qubit channel that a tenant
actually receives, replays two-qubit algorithmic workloads through that stressed
hardware model, and reports both a channel verdict and workload-level fidelity
risk.

Two analysis paths share one physical model:

* **Path A (channel).** Quantum process tomography from a tomographically
  complete probe set gives a Choi matrix. After CPTP projection it yields
  process infidelity `D`, a coherent deviation `delta_log` (the sine of the
  Fubini-Study angle between the fitted `U_eff` and the target), and an
  incoherent weight `delta_inc = 1 - p_I` from a Pauli twirl of the residual.
  Those three numbers drive a PASS / SOFT FAIL / HARD FAIL verdict.
* **Path B (workload).** Compiled two-qubit MQT Bench circuits are executed on
  the same stressed model over four probe inputs, giving expected, worst-case
  and spread output-state fidelities `F_exp`, `F_worst`, `F_spread`.

## Install

Python 3.10 or newer.

```bash
git clone https://github.com/syedshubha/QVSh.git
cd QVSh
python -m venv .venv && source .venv/bin/activate
pip install -e .
```

To reproduce the published numbers bit for bit, install the exact versions that
produced them first:

```bash
pip install -r requirements-lock.txt
pip install -e . --no-deps
```

`requirements.txt` holds the loose lower bounds; `requirements-lock.txt` holds
the exact versions recorded in `docs/ENV.txt`. `qutip-qip` is listed as the
optional `qip` extra because no module in the audited pipeline imports it.

## Quick start

```bash
qsys-audit --workload w_state --stress aggressive
qsys-audit --workload ghz --stress moderate --json
```

Equivalent without installing:

```bash
python audit_cli.py --workload w_state --stress aggressive
```

Workload aliases: `ghz`, `graph_state`, `w_state`, `qft`, `qaoa`, `vqe_su2`,
`vqe_two_local`. Stress profiles: `calibrated`, `moderate`, `aggressive`.

Output fields are `F_exp`, `F_worst`, `F_spread`, the circuit counts, the three
channel metrics `D`, `delta_log`, `delta_inc`, and the verdict.

## Repository layout

```
src/qsys_audit/
  physics_engine.py       Hamiltonian builders, QuTiP mesolve, stress profiles,
                          Gaussian single-qubit drive, echoed CR-CNOT
  reconstruction.py       QPT, CPTP projection, Kraus, polar fit, delta_log,
                          delta_inc, single-qubit classical shadows
  workload_replay.py      circuit compilation and replay, WorkloadMetrics
  cli.py                  the qsys-audit command line interface
  paths.py                config.yaml resolution
  data/config.yaml        copy of config.yaml shipped inside the wheel
experiments/              the five VeriXTalk experiments, their result CSVs and
                          figures (see experiments/README.md)
scripts/                  standalone analysis and plotting scripts
docs/                     environment record, concept-to-code inventory,
                          results, reproduction check, workflow figure
tests/                    packaging checks plus pipeline smoke tests
config.yaml               hardware limits, detunings, stress profiles,
                          calibration constants, verdict thresholds
```

## Reproducing the experiments

```bash
python experiments/run_all.py          # or: qsys-audit-experiments
python experiments/make_paper_figures.py
```

Every script is also runnable on its own and writes its CSVs next to itself and
its figures into `experiments/figs/`. Set `QVSH_OUTDIR` and `QVSH_FIGDIR` to
redirect either. The global seed is 1234 and per-trial seeds are derived
deterministically from it.

| Script | What it measures | Result CSVs | Figure |
|---|---|---|---|
| `exp1_rb.py` | gate-dependent randomized benchmarking on the stressed `H(t)`, plus a stress sweep | `figA1_rb_curve`, `figA1_rb_sweep`, `figA1_rb_summary` | `figs/figA1_rb.pdf` |
| `exp1b_irb_cnot.py` | interleaved RB of the echoed CR-CNOT and per-gate extrapolation | `figA1_irb_curve`, `figA1_irb_summary` | `figs/figA1_irb.pdf` |
| `exp2_shadow2q.py` | two-qubit shadow-assisted workload fidelity (Path B) | `figA2_shadow_curve`, `figA2_shadow_summary` | `figs/figA2_shadow.pdf` |
| `exp3_sensitivity.py` | joint `(J, A)` sensitivity sweep, Lindblad cross-check, unit mapping | `figA3_sweep`, `figA3_units`, `figA3_lindblad` | `figs/figA3_sensitivity.pdf` |
| `exp4_D_vs_F.py` | correlation of process infidelity `D` with worst-case fidelity `F_worst` | `figA4_points` | `figs/figA4_D_vs_F.pdf` |

`exp4_D_vs_F.py` reads `figA3_sweep.csv`, so run `exp3_sensitivity.py` first (or
just use `run_all.py`, which orders them correctly). `exp3_sensitivity.py` is the
longest run because it sweeps a two-dimensional `(J, A)` grid.
`make_paper_figures.py` restyles the same CSVs into `experiments/figures_paper/`
with uniform publication formatting; it recomputes nothing.

The experiments import the pipeline rather than reimplementing it: channels come
from `perform_qpt_pipeline`, workload fidelities from `benchmark_workload`, and
the verdict logic from the thresholds in `config.yaml`.

## Standalone scripts

| Script | Purpose |
|---|---|
| `scripts/run_shadow_convergence.py` | shadow convergence study for the single-qubit gate suite |
| `scripts/run_baseline_comparison.py` | static isolated-entangler baseline comparison |
| `scripts/plot_tgate_masking.py` | T-gate multi-basis diagnostic figure and table |
| `scripts/plot_workload_spectrum.py` | two-qubit workload fidelity spectrum |
| `scripts/compile_fig1.sh` | compile `docs/figure1_workflow.tex` |

## Configuration

`config.yaml` is the single place where the physics is set: hardware limits,
detunings, the three stress profiles, the single-qubit gate envelope, the echoed
CR-CNOT parameters, and the verdict thresholds. The CLI resolves it in this
order: `--config`, then `$QSYS_AUDIT_CONFIG`, then `./config.yaml`, then the
repository root, then the copy inside the installed package.

The three stress profiles and the verdict thresholds are:

* `calibrated`: all couplings and drive amplitudes zero.
* `moderate`: static coupling on, neighbor and readout drives off.
* `aggressive`: static coupling and all drive amplitudes on.
* PASS when `D <= 0.02` and `delta_log <= 0.05` and `delta_inc <= 0.05`.
* HARD FAIL when `D >= 0.30` or `delta_log >= 0.90` or `delta_inc >= 0.50`.
* SOFT FAIL otherwise.

## Documentation

* `docs/INVENTORY.md` maps every concept to the function that implements it.
* `docs/VERIFY.md` records the reproduction check against the published numbers.
* `docs/RESULTS.md` records the experiment outputs and how they were read.
* `docs/ENV.txt` records the exact environment, including the from-source QuTiP
  build, and the global seed.

## Tests

```bash
pip install -e ".[dev]"
pytest
```

The packaging tests run without the scientific stack. The pipeline smoke tests
skip themselves if QuTiP, Qiskit or MQT Bench is missing.

## Citation

See `CITATION.cff`.

## License

MIT. See `LICENSE`.
