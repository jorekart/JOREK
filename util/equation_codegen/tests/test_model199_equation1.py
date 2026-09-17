import unittest

import sympy as sp

from jorek_equations import coefficient, dR, dZ, dphi, expand_brackets, variation
from jorek_equations.model199 import (
    F0,
    FIELDS,
    eta,
    induction_equation_1,
    j,
    omega,
    psi,
    rho,
    T,
    u,
)


class Model199Equation1Test(unittest.TestCase):
    def setUp(self):
        self.timestep = coefficient("tstep_model199_eq1_test")
        self.theta = coefficient("theta_model199_eq1_test")
        self.zeta = coefficient("zeta_model199_eq1_test")
        self.xjac = coefficient("xjac")
        self.R = sp.Symbol("R")

    def test_field_order_is_the_model_order(self):
        self.assertEqual([value.name for value in FIELDS], [
            "psi", "u", "j", "omega", "rho", "T"
        ])

    def test_rhs_matches_fortran_equation_1_terms(self):
        equation = induction_equation_1()
        result = equation.linearize(
            fields=(psi, u, j, T),
            timestep=self.timestep,
            theta=self.theta,
            zeta=self.zeta,
        )
        v = equation.test
        current_source = sp.Symbol("current_source")
        eta_num = sp.Symbol("eta_num")
        eps_cyl = sp.Symbol("eps_cyl")

        expected_B = (
            v * eta(T.current) * (j.current - current_source) / self.R * self.xjac
            + v * self.xjac * (
                dR(psi.current) * dZ(u.current)
                - dZ(psi.current) * dR(u.current)
            )
            - v * eps_cyl * F0 / self.R * dphi(u.current) * self.xjac
            + eta_num * self.xjac * (
                dR(v) * dR(j.current) + dZ(v) * dZ(j.current)
            )
        )
        expected = self.timestep * expected_B + self.zeta * (
            v * psi.previous_increment / self.R * self.xjac
        )
        self.assertEqual(sp.expand(expand_brackets(result.rhs) - expected), 0)

    def test_definition_keeps_fields_abstract(self):
        equation = induction_equation_1()
        self.assertIn(psi.symbol, equation.A.atoms(type(psi.symbol)))
        self.assertIn(T.symbol, equation.B.atoms(type(T.symbol)))
        self.assertNotIn(psi.current, equation.A.atoms(type(psi.current)))

    def test_amat_blocks_match_fortran_equation_1_terms(self):
        equation = induction_equation_1()
        result = equation.linearize(
            fields=(psi, u, j, T),
            timestep=self.timestep,
            theta=self.theta,
            zeta=self.zeta,
        )
        v = equation.test
        eta_num = sp.Symbol("eta_num")
        eps_cyl = sp.Symbol("eps_cyl")
        current_source = sp.Symbol("current_source")

        expected_psi = (
            (1 + self.zeta) * v * psi.trial / self.R * self.xjac
            - self.theta * self.timestep * v * self.xjac * (
                dR(psi.trial) * dZ(u.current)
                - dZ(psi.trial) * dR(u.current)
            )
        )
        expected_u = (
            -self.theta * self.timestep * v * self.xjac * (
                dR(psi.current) * dZ(u.trial)
                - dZ(psi.current) * dR(u.trial)
            )
            + self.theta * self.timestep * v * eps_cyl * F0 / self.R
            * dphi(u.trial) * self.xjac
        )
        expected_j = -self.theta * self.timestep * (
            eta_num * self.xjac * (
                dR(v) * dR(j.trial) + dZ(v) * dZ(j.trial)
            )
            + eta(T.current) * v * j.trial / self.R * self.xjac
        )
        expected_T = -self.theta * self.timestep * (
            sp.Function("deta_dT")(T.current) * T.trial * v
            * (j.current - current_source) / self.R * self.xjac
        )

        self.assertEqual(
            sp.expand(expand_brackets(result.amat[psi]) - expected_psi), 0
        )
        self.assertEqual(
            sp.expand(expand_brackets(result.amat[u]) - expected_u), 0
        )
        self.assertEqual(sp.expand(result.amat[j] - expected_j), 0)
        self.assertEqual(sp.expand(result.amat[T] - expected_T), 0)

    def test_unrelated_model_fields_have_zero_blocks(self):
        equation = induction_equation_1()
        result = equation.linearize(
            fields=FIELDS,
            timestep=self.timestep,
            theta=self.theta,
            zeta=self.zeta,
        )
        self.assertEqual(result.amat[omega], 0)
        self.assertEqual(result.amat[rho], 0)


if __name__ == "__main__":
    unittest.main()
