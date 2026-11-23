"""Generate all figures."""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
MODULES = [
    "model_experiment_compare",
    "phir_gex",
    "destabilized_GFP_cost",
    "sensitivity_CtfCtx",
    "sensitivity_cost",
    "gphi_ratio",
    "fold_change_construct_condition",
    "destabilized_lacZ_cost",
]


def main():
    for module in MODULES:
        print(f"\n===== plotting/{module}.py =====", flush=True)
        subprocess.run([sys.executable, "-m", f"plotting.{module}"], cwd=ROOT, check=True)

    print("\nAll figures were written to", ROOT / "output")


if __name__ == "__main__":
    main()
