#!/usr/bin/env python3
"""Compare model-600 var_u psi/w AMAT columns."""

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PROJECT_ROOT.parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from jorek_equations.model600_compare import compare_model600_u_psi_w  # noqa: E402
from jorek_equations.source_compare import FortranComparisonError  # noqa: E402


def main():
    path = REPOSITORY_ROOT / "models/model600/mod_elt_matrix_fft.f90"
    try:
        compare_model600_u_psi_w(path)
    except FortranComparisonError as error:
        print("model600 var_u psi/w AMAT differences:")
        print(error)
        raise SystemExit(1)
    print("No differences found for model600 var_u psi/w AMAT blocks")


if __name__ == "__main__":
    main()
