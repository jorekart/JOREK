#!/usr/bin/env python3
"""Compare the single-temperature model-600 psi AMAT row only."""

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PROJECT_ROOT.parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from jorek_equations.model600_compare import (  # noqa: E402
    compare_model600_amat_psi_single_temperature,
)
from jorek_equations.source_compare import FortranComparisonError  # noqa: E402


def main():
    path = REPOSITORY_ROOT / "models/model600/mod_elt_matrix_fft.f90"
    try:
        compare_model600_amat_psi_single_temperature(path)
    except FortranComparisonError as error:
        print("Single-temperature amat(var_psi,:) differences:")
        print(error)
        raise SystemExit(1)
    print("No differences found for single-temperature amat(var_psi,:)")


if __name__ == "__main__":
    main()
