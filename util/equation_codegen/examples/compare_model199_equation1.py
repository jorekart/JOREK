#!/usr/bin/env python3
"""Compare generated model-199 equation 1 with the JOREK Fortran source."""

import argparse
from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PROJECT_ROOT.parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from jorek_equations.source_compare import (  # noqa: E402
    EQUATION_1_ASSIGNMENTS,
    FortranComparisonError,
    compare_assignment_maps,
    extract_fortran_assignments,
    generated_assignment_text,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "fortran_file",
        nargs="?",
        type=Path,
        default=REPOSITORY_ROOT / "models/model199/mod_elt_matrix_fft.f90",
        help="path to model199/mod_elt_matrix_fft.f90",
    )
    arguments = parser.parse_args()
    source = extract_fortran_assignments(
        arguments.fortran_file, EQUATION_1_ASSIGNMENTS
    )
    generated = generated_assignment_text()

    try:
        compare_assignment_maps(source, generated)
    except FortranComparisonError as error:
        print(str(error))
        raise SystemExit(1)
    print("PASS: all equation-1 RHS and AMAT expressions are symbolically equal")


if __name__ == "__main__":
    main()
