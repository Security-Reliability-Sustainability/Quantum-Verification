#!/usr/bin/env python3
"""Backwards-compatible wrapper so that `python audit_cli.py ...` keeps working.

The implementation now lives in src/qsys_audit/cli.py and is also installed as
the `qsys-audit` console script.
"""

from __future__ import annotations

import os
import sys

_SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "src")
if os.path.isdir(_SRC) and _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from qsys_audit.cli import main  # noqa: E402

if __name__ == "__main__":
    main()
