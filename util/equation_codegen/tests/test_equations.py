import unittest

import sympy as sp

from jorek_equations import (
    ConstraintEquation,
    EvolutionEquation,
    coefficient,
    dR,
    external_function,
    field,
    test_function,
)


class EquationConventionTest(unittest.TestCase):
    def test_evolution_rhs_and_amat_signs(self):
        psi = field("psi_evolution_test")
        j = field("j_evolution_test")
        T = field("T_evolution_test")
        v = test_function("v_evolution_test")
        Rloc = coefficient("R_evolution_test")
        timestep = coefficient("tstep_evolution_test")
        theta = coefficient("theta_evolution_test")
        zeta = coefficient("zeta_evolution_test")
        eta = external_function(
            "eta_evolution_test",
            arguments=("T",),
            derivatives={"T": "deta_dT_evolution_test"},
        )

        A = v * psi.current / Rloc
        B = v * eta(T.current) * j.current / Rloc
        result = EvolutionEquation("induction_test", v, A, B).linearize(
            fields=(psi, j, T),
            timestep=timestep,
            theta=theta,
            zeta=zeta,
        )

        expected_rhs = (
            timestep * B
            + zeta * v * psi.previous_increment / Rloc
        )
        expected_psi = (1 + zeta) * v * psi.trial / Rloc
        expected_j = -theta * timestep * v * eta(T.current) * j.trial / Rloc
        expected_T = (
            -theta
            * timestep
            * v
            * sp.Function("deta_dT_evolution_test")(T.current)
            * T.trial
            * j.current
            / Rloc
        )

        self.assertEqual(sp.expand(result.rhs - expected_rhs), 0)
        self.assertEqual(sp.expand(result.amat[psi] - expected_psi), 0)
        self.assertEqual(sp.expand(result.amat[j] - expected_j), 0)
        self.assertEqual(sp.expand(result.amat[T] - expected_T), 0)

    def test_constraint_uses_newton_residual_sign(self):
        psi = field("psi_constraint_test")
        j = field("j_constraint_test")
        v = test_function("v_constraint_test")
        Rloc = coefficient("R_constraint_test")
        C = (dR(v) * dR(psi.current) + v * j.current) / Rloc

        result = ConstraintEquation("current_definition_test", v, C).linearize(
            fields=(psi, j)
        )

        self.assertEqual(result.rhs, -C)
        self.assertEqual(result.kind, "static")
        self.assertEqual(result.amat[psi], dR(v) * dR(psi.trial) / Rloc)
        self.assertEqual(result.amat[j], v * j.trial / Rloc)


if __name__ == "__main__":
    unittest.main()
