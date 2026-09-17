import unittest

from jorek_equations import dR, dphi, fortran
from jorek_equations.model199 import eta, j, psi, u


class FortranPrinterTest(unittest.TestCase):
    def test_model199_field_and_derivative_names(self):
        self.assertEqual(fortran(dR(j.current)), "zj0_x")
        self.assertEqual(fortran(dphi(u.trial)), "u_p")
        self.assertEqual(fortran(eta(psi.current)), "eta_T")

    def test_previous_increment_can_use_assembly_array_name(self):
        self.assertEqual(
            fortran(
                psi.previous_increment,
                previous_names={"psi": "delta_g(mp,1,ms,mt)"},
            ),
            "delta_g(mp,1,ms,mt)",
        )


if __name__ == "__main__":
    unittest.main()

