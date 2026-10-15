# Convenience targets. Run `make help` for the list.
PY ?= python3
VENV ?= .venv

.PHONY: help install install-dev install-locked audit experiments figures test clean-pyc

help:
	@grep -E '^[a-zA-Z_-]+:.*?## ' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "  %-16s %s\n", $$1, $$2}'

install: ## install the package (editable) with its runtime dependencies
	$(PY) -m pip install -e .

install-dev: ## install with the test extras
	$(PY) -m pip install -e ".[dev]"

install-locked: ## install the exact versions used to produce the published results
	$(PY) -m pip install -r requirements-lock.txt && $(PY) -m pip install -e . --no-deps

audit: ## one example audit (W state, aggressive stress)
	$(PY) -m qsys_audit.cli --workload w_state --stress aggressive

experiments: ## run all five VeriXTalk experiments (regenerates CSVs and figs/)
	$(PY) experiments/run_all.py

figures: ## rebuild the styled paper figures from the existing CSVs
	$(PY) experiments/make_paper_figures.py

test: ## run the test suite
	$(PY) -m pytest

clean-pyc: ## remove bytecode caches
	find . -name '__pycache__' -type d -prune -exec rm -rf {} +
