# Experiments

Five experiments, each importing the `qsys_audit` pipeline rather than
reimplementing any physics. Global seed 1234; per-trial seeds are derived
deterministically from it and recorded in each CSV and script header.

Run everything in dependency order:

```bash
python experiments/run_all.py
```

Or one at a time:

```bash
python experiments/exp1_rb.py
python experiments/exp1b_irb_cnot.py
python experiments/exp2_shadow2q.py
python experiments/exp3_sensitivity.py     # must precede exp4
python experiments/exp4_D_vs_F.py
```

## Files

* `exp_common.py` shared helpers: config loading, engine construction, the
  verdict classifier, the stressed-channel Kraus helper, and the 24-element
  single-qubit Clifford group. Also defines the path constants `ROOT`, `HERE`,
  `OUTDIR` and `FIGDIR`.
* `exp1_rb.py`, `exp1b_irb_cnot.py`, `exp2_shadow2q.py`, `exp3_sensitivity.py`,
  `exp4_D_vs_F.py` the experiments themselves.
* `run_all.py` runs all five in order.
* `make_paper_figures.py` restyles the CSVs into `figures_paper/` with uniform
  publication formatting. It recomputes nothing, so it is safe to rerun.
* `*.csv` the committed results, so the figures can be rebuilt without rerunning
  the simulations.
* `figs/` raw figures written by the experiments.
* `figures_paper/` restyled figures written by `make_paper_figures.py`.

## Output locations

CSVs go to `OUTDIR` (this directory by default) and figures to `FIGDIR`
(`figs/` by default). Override with the environment variables `QVSH_OUTDIR` and
`QVSH_FIGDIR`. `make_paper_figures.py` reads `DATA_DIR` and writes `OUT_DIR`.
