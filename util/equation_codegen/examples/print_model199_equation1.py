#!/usr/bin/env python3
"""Print the symbolic model-199 equation-1 linearization."""

from pathlib import Path
import sys

import sympy as sp


# Make the example runnable directly from a source checkout.
PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from jorek_equations import coefficient, fortran  # noqa: E402
from jorek_equations.model199 import (  # noqa: E402
    FIELDS,
    induction_equation_1,
)


def print_terms(expression, indent="  "):
    """Print an expanded sum with one signed term per line."""

    terms = sp.Add.make_args(sp.expand(expression))
    for term in terms:
        negative = term.could_extract_minus_sign()
        magnitude = -term if negative else term
        sign = "-" if negative else "+"
        print("{}{} {}".format(indent, sign, sp.sstr(magnitude)))


def print_fortran_terms(expression, indent="  "):
    """Print an expanded sum in JOREK-style Fortran notation."""

    terms = sp.Add.make_args(sp.expand(expression))
    for term in terms:
        negative = term.could_extract_minus_sign()
        magnitude = -term if negative else term
        sign = "-" if negative else "+"
        print(
            "{}{} {}".format(
                indent,
                sign,
                fortran(
                    magnitude,
                    previous_names={"psi": "delta_g(mp,1,ms,mt)"},
                ),
            )
        )


def main() -> None:
    equation = induction_equation_1()
    timestep = coefficient("tstep")
    theta = coefficient("theta")
    zeta = coefficient("zeta")
    result = equation.linearize(
        fields=FIELDS,
        timestep=timestep,
        theta=theta,
        zeta=zeta,
    )

    print("Model 199, equation 1: induction")
    print("=" * 40)
    print("Weak A:")
    print_terms(equation.A)
    print("Weak B:")
    print_terms(equation.B)
    print("\nLinearized RHS:")
    print_terms(result.rhs)

    print("\nLinearized AMAT blocks:")
    for value in FIELDS:
        block = result.amat[value]
        if block != 0:
            print("  {}:".format(value.name))
            print_terms(block, indent="    ")

    print("\nFortran-style linearized terms:")
    print("RHS:")
    print_fortran_terms(result.rhs, indent="  ")
    for value in FIELDS:
        block = result.amat[value]
        if block != 0:
            print("{} AMAT:".format(value.name))
            print_fortran_terms(block, indent="  ")


if __name__ == "__main__":
    main()
