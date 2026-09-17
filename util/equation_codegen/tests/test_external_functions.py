import unittest

import sympy as sp

from jorek_equations import (
    FieldRole,
    MissingDerivativeError,
    external_function,
    field,
    variation,
)


class ExternalFunctionsTest(unittest.TestCase):
    def test_active_function_uses_supplied_derivative(self):
        T = field("T_external_active_test")
        j = field("j_external_active_test")
        eta = external_function(
            "eta_external_active_test",
            arguments=("T",),
            derivatives={"T": "deta_dT_external_active_test"},
        )
        expression = eta(T.current) * j.current

        expected_T = (
            sp.Function("deta_dT_external_active_test")(T.current)
            * T.trial
            * j.current
        )
        expected_j = eta(T.current) * j.trial

        self.assertEqual(
            variation(expression, T, direction=FieldRole.TRIAL), expected_T
        )
        self.assertEqual(
            variation(expression, j, direction=FieldRole.TRIAL), expected_j
        )

    def test_piecewise_active_uses_branch_derivative_interface(self):
        T = field("T_external_piecewise_test")
        eta = external_function(
            "eta_external_piecewise_test",
            arguments=("T",),
            derivatives={"T": "deta_dT_external_piecewise_test"},
            policy="piecewise_active",
        )

        expected = (
            sp.Function("deta_dT_external_piecewise_test")(T.current)
            * T.increment
        )
        self.assertEqual(variation(eta(T.current), T), expected)

    def test_frozen_external_function_does_not_propagate_dependency(self):
        psi = field("psi_external_frozen_test")
        Dperp = external_function(
            "Dperp_external_frozen_test",
            arguments=("psi_norm",),
            derivatives={},
            policy="frozen",
        )

        self.assertEqual(variation(Dperp(psi.current), psi), 0)

    def test_missing_active_derivative_is_an_error(self):
        psi = field("psi_external_missing_test")
        profile = external_function(
            "profile_external_missing_test",
            arguments=("psi",),
            derivatives={},
            policy="active",
        )

        with self.assertRaises(MissingDerivativeError):
            variation(profile(psi.current), psi)


if __name__ == "__main__":
    unittest.main()

