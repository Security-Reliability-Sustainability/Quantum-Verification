"""Console-script shim that runs the full VeriXTalk experiment suite.

The experiment scripts ship with the source repository rather than inside the
wheel, because they read and write result CSVs and figures next to themselves.
This runner locates them in a source checkout and executes them in order.
"""

from __future__ import annotations

import runpy
import sys

from .paths import experiments_dir


def main() -> None:
    exp_dir = experiments_dir()
    if exp_dir is None:
        print(
            "The experiment scripts are not part of the installed wheel.\n"
            "Clone the repository and run:\n"
            "    python experiments/run_all.py",
            file=sys.stderr,
        )
        raise SystemExit(1)

    if str(exp_dir) not in sys.path:
        sys.path.insert(0, str(exp_dir))
    runpy.run_path(str(exp_dir / "run_all.py"), run_name="__main__")


if __name__ == "__main__":  # pragma: no cover
    main()
