import unittest

import sympy as sp

from jorek_equations import (
    FieldRole,
    bracket,
    coefficient,
    dR,
    dZ,
    field,
    freeze,
    grad,
    test_function,
    variation,
)
from jorek_equations.model199 import FIELDS


class SymbolsAndOperatorsTest(unittest.TestCase):
    def test_model199_field_order_and_fortran_names(self):
        self.assertEqual(
            [value.name for value in FIELDS],
            ["psi", "u", "j", "omega", "rho", "T"],
        )
        self.assertEqual(FIELDS[2].fortran_current, "zj0")
        self.assertEqual(FIELDS[3].fortran_trial, "w")

    def test_roles_are_distinct(self):
        psi = field("psi_role_test")
        self.assertEqual(psi.current.role, FieldRole.CURRENT)
        self.assertEqual(psi.increment.role, FieldRole.DELTA)
        self.assertEqual(psi.previous_increment.role, FieldRole.PREVIOUS_DELTA)
        self.assertEqual(psi.trial.role, FieldRole.TRIAL)
        self.assertNotEqual(psi.current, psi.trial)

    def test_bare_field_is_abstract_until_linearization(self):
        psi = field("psi_abstract_test")
        u = field("u_abstract_test")
        expression = psi * u

        self.assertEqual(variation(expression, psi), psi.increment * u.current)
        self.assertEqual(variation(expression, u), psi.current * u.increment)

    def test_derivative_variation_commutes(self):
        psi = field("psi_derivative_test")
        actual = variation(dR(psi.current), psi)
        self.assertEqual(actual, dR(psi.increment))

    def test_bracket_product_rule(self):
        psi = field("psi_bracket_test")
        u = field("u_bracket_test")
        expression = bracket(psi.current, u.current)

        self.assertEqual(
            variation(expression, psi, direction=FieldRole.TRIAL),
            bracket(psi.trial, u.current),
        )
        self.assertEqual(
            variation(expression, u, direction=FieldRole.TRIAL),
            bracket(psi.current, u.trial),
        )

    def test_gradient_dot_product_rule(self):
        u = field("u_gradient_test")
        expression = sum(component**2 for component in grad(u.current))
        expected = 2 * (
            dR(u.current) * dR(u.trial)
            + dZ(u.current) * dZ(u.trial)
        )
        actual = variation(expression, u, direction=FieldRole.TRIAL)
        self.assertEqual(sp.expand(actual - expected), 0)

    def test_frozen_expression_has_zero_variation(self):
        rho = field("rho_frozen_test")
        u = field("u_frozen_test")
        expression = freeze(rho.current) * dR(u.current)

        self.assertEqual(variation(expression, rho), 0)
        self.assertEqual(
            variation(expression, u),
            freeze(rho.current) * dR(u.increment),
        )

    def test_test_function_and_coefficients_do_not_vary(self):
        psi = field("psi_test_function_test")
        v = test_function("v_test_function_test")
        Rloc = coefficient("R_test_function_test")
        expression = v * psi.current / Rloc

        expected = v * psi.trial / Rloc
        actual = variation(expression, psi, direction=FieldRole.TRIAL)
        self.assertEqual(actual, expected)


if __name__ == "__main__":
    unittest.main()
