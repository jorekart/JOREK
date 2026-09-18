#!/usr/bin/env python3
"""Compare model-600 current and vorticity definition equations."""

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PROJECT_ROOT.parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from jorek_equations.model600_compare import compare_model600_simple_equations  # noqa: E402
from jorek_equations.source_compare import FortranComparisonError  # noqa: E402


def main():
    path = REPOSITORY_ROOT / "models/model600/mod_elt_matrix_fft.f90"
    try:
        compare_model600_simple_equations(path)
    except FortranComparisonError as error:
        print("Model-600 var_zj/var_w differences:")
        print(error)
        raise SystemExit(1)
    print("No differences found for model-600 var_zj and var_w equations")


if __name__ == "__main__":
    main()
