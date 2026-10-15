"""Filesystem locations used by qsys_audit.

Resolution order for the configuration file:

1. the path passed explicitly (``--config`` on the command line),
2. ``$QSYS_AUDIT_CONFIG``,
3. ``./config.yaml`` in the current working directory,
4. ``config.yaml`` at the root of a source checkout,
5. the copy shipped inside the installed package.

Steps 3 to 5 make the package usable both from a cloned repository and from a
wheel installed into a virtual environment.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

PACKAGE_DIR = Path(__file__).resolve().parent
PACKAGED_CONFIG = PACKAGE_DIR / "data" / "config.yaml"


def repo_root() -> Optional[Path]:
    """Return the root of the source checkout, or None if running from a wheel."""
    for parent in PACKAGE_DIR.parents:
        if (parent / "pyproject.toml").is_file() and (parent / "config.yaml").is_file():
            return parent
    return None


def default_config() -> Path:
    """Locate config.yaml. Raises FileNotFoundError if nothing is found."""
    env = os.environ.get("QSYS_AUDIT_CONFIG")
    if env:
        candidate = Path(env).expanduser()
        if not candidate.is_file():
            raise FileNotFoundError(
                f"QSYS_AUDIT_CONFIG points at a missing file: {candidate}"
            )
        return candidate

    cwd_config = Path.cwd() / "config.yaml"
    if cwd_config.is_file():
        return cwd_config

    root = repo_root()
    if root is not None:
        return root / "config.yaml"

    if PACKAGED_CONFIG.is_file():
        return PACKAGED_CONFIG

    raise FileNotFoundError(
        "no config.yaml found. Pass --config explicitly or set QSYS_AUDIT_CONFIG."
    )


def experiments_dir() -> Optional[Path]:
    """Return the experiments/ directory of a source checkout, or None."""
    root = repo_root()
    if root is None:
        return None
    candidate = root / "experiments"
    return candidate if candidate.is_dir() else None
