#!/usr/bin/env python3
"""Compare model-600 equation 1 (the psi equation) with Fortran."""

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PROJECT_ROOT.parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from jorek_equations.model600_compare import (  # noqa: E402
    compare_model600_amat_psi_psi,
    compare_model600_amat_psi_remaining,
    compare_model600_rhs_psi,
)
from jorek_equations.source_compare import FortranComparisonError  # noqa: E402


def main():
    path = REPOSITORY_ROOT / "models/model600/mod_elt_matrix_fft.f90"
    checks = (
        ("rhs(var_psi)", compare_model600_rhs_psi),
        ("amat(var_psi,var_psi)", compare_model600_amat_psi_psi),
        ("remaining amat(var_psi,:)", compare_model600_amat_psi_remaining),
    )
    failed = False
    for label, check in checks:
        try:
            check(path)
        except FortranComparisonError as error:
            print("{} differences:".format(label))
            print(error)
            failed = True
        else:
            print("No differences found for model-600 {}".format(label))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
