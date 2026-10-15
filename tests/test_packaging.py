"""Packaging and layout checks. These run without the scientific stack."""

from __future__ import annotations

import compileall
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SKIP_DIRS = {".git", "_private", "brisbane_multivariate_assignment", "__pycache__",
             ".venv", "venv", "build", "dist", ".pytest_cache"}


def _iter_files(suffixes=None, name_glob=None):
    """Walk the repository, pruning directories we never publish."""
    for dirpath, dirnames, filenames in os.walk(ROOT):
        dirnames[:] = [
            d for d in dirnames if d not in SKIP_DIRS and not d.endswith(".egg-info")
        ]
        for filename in filenames:
            path = Path(dirpath) / filename
            if suffixes is not None and path.suffix not in suffixes:
                continue
            if name_glob is not None and not path.match(name_glob):
                continue
            yield path


def _pyproject_version() -> str:
    text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    try:
        import tomllib  # Python 3.11+
    except ModuleNotFoundError:
        match = re.search(r'^version\s*=\s*"([^"]+)"', text, re.MULTILINE)
        assert match, "no version field in pyproject.toml"
        return match.group(1)
    return tomllib.loads(text)["project"]["version"]


def test_package_imports_without_scientific_stack():
    sys.path.insert(0, str(ROOT / "src"))
    import qsys_audit

    assert qsys_audit.__version__
    assert callable(qsys_audit.default_config)


def test_version_matches_pyproject():
    sys.path.insert(0, str(ROOT / "src"))
    import qsys_audit

    assert qsys_audit.__version__ == _pyproject_version()


def test_packaged_config_matches_repository_config():
    root_config = (ROOT / "config.yaml").read_bytes()
    packaged = (ROOT / "src" / "qsys_audit" / "data" / "config.yaml").read_bytes()
    assert root_config == packaged, (
        "src/qsys_audit/data/config.yaml has drifted from config.yaml. "
        "Run: cp config.yaml src/qsys_audit/data/config.yaml"
    )


def test_console_script_targets_are_callable():
    sys.path.insert(0, str(ROOT / "src"))
    from qsys_audit import cli, experiments_runner

    assert callable(cli.main)
    assert callable(experiments_runner.main)


def test_cli_help_runs():
    env_src = str(ROOT / "src")
    proc = subprocess.run(
        [sys.executable, "-c",
         f"import sys; sys.path.insert(0, {env_src!r});"
         " from qsys_audit.cli import build_parser; build_parser().parse_args(['--help'])"],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr
    assert "--stress" in proc.stdout


@pytest.mark.parametrize(
    "stale", ["rebuttal_experiments", "rebuttal_common", "rebuttal_package"]
)
def test_no_stale_path_references(stale):
    this_file = Path(__file__).resolve()
    offenders = [
        str(p.relative_to(ROOT))
        for p in _iter_files({".py", ".sh", ".toml", ".cfg", ".md", ".txt", ".yaml", ".yml"})
        if p.resolve() != this_file
        and stale in p.read_text(encoding="utf-8", errors="replace")
    ]
    assert not offenders, f"{stale} still referenced in: {offenders}"


def test_confidential_material_is_not_outside_private():
    leaked = [
        str(p.relative_to(ROOT))
        for p in _iter_files(name_glob="Reviews_*")
    ]
    assert not leaked, f"review material outside _private/: {leaked}"


def test_gitignore_blocks_what_must_not_be_published():
    text = (ROOT / ".gitignore").read_text(encoding="utf-8")
    for entry in ("_private/", "Reviews_*.txt", "brisbane_multivariate_assignment/",
                  "__pycache__/", ".DS_Store"):
        assert entry in text, f"{entry} missing from .gitignore"


def test_every_module_byte_compiles():
    ok = True
    for directory in ("src", "experiments", "scripts", "tests"):
        ok &= bool(compileall.compile_dir(
            str(ROOT / directory), quiet=2, force=True, legacy=False
        ))
    assert compileall.compile_file(str(ROOT / "audit_cli.py"), quiet=2, force=True)
    assert ok


def test_verdict_classifier_matches_published_thresholds():
    sys.path.insert(0, str(ROOT / "src"))
    import yaml

    from qsys_audit.cli import classify_verdict

    cfg = yaml.safe_load((ROOT / "config.yaml").read_text(encoding="utf-8"))
    assert classify_verdict(0.0001, 0.0098, 0.0000, cfg) == "PASS"
    assert classify_verdict(0.2687, 0.1065, 0.2605, cfg) == "SOFT FAIL"
    assert classify_verdict(0.3377, 0.3907, 0.2618, cfg) == "HARD FAIL"
