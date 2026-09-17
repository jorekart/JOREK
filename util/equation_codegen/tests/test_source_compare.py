import unittest
from pathlib import Path

from jorek_equations import coefficient, split_toroidal_channels
from jorek_equations.model199 import FIELDS, induction_equation_1, u
from jorek_equations.source_compare import (
    EQUATION_1_ASSIGNMENTS,
    FortranComparisonError,
    compare_assignment_maps,
    compare_model199_equation1,
    extract_fortran_assignments,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PROJECT_ROOT.parents[1]
FORTRAN_FILE = REPOSITORY_ROOT / "models/model199/mod_elt_matrix_fft.f90"


class SourceComparisonTest(unittest.TestCase):
    def test_extracts_all_equation_1_assignments(self):
        assignments = extract_fortran_assignments(
            FORTRAN_FILE, EQUATION_1_ASSIGNMENTS
        )
        self.assertEqual(set(assignments), set(EQUATION_1_ASSIGNMENTS))
        self.assertIn("eta_T", assignments["rhs_ij_1"])
        self.assertIn("deta_dT", assignments["amat_16"])

    def test_real_model199_equation_1_matches(self):
        compare_model199_equation1(FORTRAN_FILE)

    def test_toroidal_trial_term_is_in_n_channel(self):
        result = induction_equation_1().linearize(
            fields=FIELDS,
            timestep=coefficient("tstep"),
            theta=coefficient("theta"),
            zeta=coefficient("zeta"),
        )
        channels = split_toroidal_channels(result.amat[u])
        self.assertNotEqual(channels.p, 0)
        self.assertNotEqual(channels.n, 0)
        self.assertEqual(channels.k, 0)
        self.assertEqual(channels.kn, 0)

    def test_mismatch_raises_and_reports_difference(self):
        with self.assertRaises(FortranComparisonError) as raised:
            compare_assignment_maps(
                {"amat_test": "x + y"},
                {"amat_test": "x - y"},
            )
        message = str(raised.exception)
        self.assertIn("amat_test", message)
        self.assertIn("source-only:", message)
        self.assertIn("generated-only:", message)
        self.assertIn("      + y", message)
        self.assertIn("      - y", message)


if __name__ == "__main__":
    unittest.main()
