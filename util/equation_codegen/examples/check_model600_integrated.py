#!/usr/bin/env python3
"""Run the integrated model-600 checks for psi, u, zj, and w."""

from pathlib import Path
import argparse
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PROJECT_ROOT.parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from jorek_equations.model600_compare import (  # noqa: E402
    compare_model600_amat_psi_psi,
    compare_model600_amat_psi_remaining,
    compare_model600_rhs_psi,
    compare_model600_simple_equations,
    compare_model600_u_base_zj,
    compare_model600_u_psi_w,
    check_model600_u_neo_generation,
    check_model600_u_neutral_generation,
    check_model600_u_neutral_amat_generation,
    check_model600_u_neutral_source_generation,
    check_model600_u_delta_n_convection_generation,
    check_model600_u_impurity_pressure_generation,
    check_model600_u_full_amat_generation,
)
from jorek_equations.model600_markdown import export_model600_markdown  # noqa: E402
from jorek_equations.source_compare import FortranComparisonError  # noqa: E402


# Keep this disabled while the Markdown term reports are being inspected.
# Set to True to run the symbolic pass/fail checks as well.
RUN_CHECKS = False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--equation",
        choices=("psi", "u", "zj", "w"),
        action="append",
        help="export only this equation; repeat for multiple equations",
    )
    args = parser.parse_args()
    path = REPOSITORY_ROOT / "models/model600/mod_elt_matrix_fft.f90"
    reports = PROJECT_ROOT / "reports"
    failed = False
    if RUN_CHECKS:
        print("WARNING: NEO AMAT terms are currently excluded from comparison.")
        print("         Only NEO residual construction is checked.\n")
        checks = (
            ("psi RHS", compare_model600_rhs_psi),
            ("psi AMAT(var_psi,var_psi)", compare_model600_amat_psi_psi),
            ("psi remaining AMAT row", compare_model600_amat_psi_remaining),
            ("zj and w equations", compare_model600_simple_equations),
            ("u AMAT(var_u,var_zj)", compare_model600_u_base_zj),
            ("u AMAT(var_u,var_psi/var_w)", compare_model600_u_psi_w),
            ("u NEO residual generation", lambda path: check_model600_u_neo_generation()),
            ("u neutral residual generation", lambda path: check_model600_u_neutral_generation()),
            ("u neutral AMAT generation", lambda path: check_model600_u_neutral_amat_generation()),
            ("u particle/source contributions", lambda path: check_model600_u_neutral_source_generation()),
            ("u delta_n_convection switches", lambda path: check_model600_u_delta_n_convection_generation()),
            ("u impurity-pressure terms", lambda path: check_model600_u_impurity_pressure_generation()),
            ("u full non-NEO AMAT generation", lambda path: check_model600_u_full_amat_generation()),
        )
        for label, check in checks:
            try:
                check(path)
            except FortranComparisonError as error:
                print("\n{} differences:".format(label))
                print(error)
                failed = True
            else:
                print("OK: {}".format(label))
    else:
        print("Report-only mode: symbolic comparison checks are disabled (RUN_CHECKS=False).")
    source_report, generated_report = export_model600_markdown(
        path,
        reports / "model600_fortran_terms.md",
        reports / "model600_generated_terms.md",
        equations=args.equation,
    )
    print("\nWrote aligned term reports:")
    print("  {}".format(source_report))
    print("  {}".format(generated_report))
    print("Compare with: meld {} {}".format(source_report, generated_report))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
