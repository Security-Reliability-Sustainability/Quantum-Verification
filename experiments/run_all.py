"""Run all VeriXTalk experiments in order. Each script is runnable individually too:
    python experiments/exp1_rb.py    (etc.)

Equivalent console entry point after installation:
    qsys-audit-experiments
"""
import runpy
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

SCRIPTS = ["exp1_rb.py", "exp1b_irb_cnot.py", "exp2_shadow2q.py",
           "exp3_sensitivity.py", "exp4_D_vs_F.py"]

def main() -> None:
    for s in SCRIPTS:
        print("\n" + "=" * 70 + f"\nRUNNING {s}\n" + "=" * 70)
        runpy.run_path(str(HERE / s), run_name="__main__")
    print("\nAll experiments complete. See experiments/figs/ and experiments/*.csv, "
          "then docs/RESULTS.md.")


if __name__ == "__main__":
    main()
