"""Checks for the report comparison helpers."""

import sys
import unittest
from pathlib import Path

import sympy as sp

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))
sys.path.insert(0, str(PROJECT_ROOT / "examples"))

from diff_model600_reports import to_element_basis, to_physical  # noqa: E402
from jorek_equations.model600_markdown import (  # noqa: E402
    SourceAssignment,
    _source_monomial_slots,
)


class CoordinateSpellingTest(unittest.TestCase):
    """The element routine writes poloidal brackets both ways."""

    def _both_spellings(self):
        a_s, a_t, b_s, b_t = sp.symbols("a_s a_t b_s b_t")
        a_x, a_y, b_x, b_y = sp.symbols("a_x a_y b_x b_y")
        xjac = sp.Symbol("xjac")
        return a_s * b_t - a_t * b_s, xjac * (a_x * b_y - a_y * b_x)

    def test_element_basis_identifies_the_two_spellings(self):
        element, physical = self._both_spellings()
        self.assertEqual(to_element_basis(element - physical), 0)

    def test_display_basis_identifies_the_two_spellings(self):
        element, physical = self._both_spellings()
        self.assertEqual(to_physical(element - physical), 0)

    def test_a_genuine_difference_survives(self):
        element, physical = self._both_spellings()
        self.assertNotEqual(to_element_basis(element - 2 * physical), 0)

    def test_variable_names_are_not_mangled(self):
        # ``var_t`` is an assignment name, not a derivative, and ``Sion_T``
        # ends in an upper-case T.
        expression = sp.Symbol("var_t") * sp.Symbol("Sion_T")
        self.assertEqual(to_element_basis(expression), expression)
        self.assertEqual(to_physical(expression), expression)


class CoefficientAlignmentTest(unittest.TestCase):
    """A coefficient mismatch must stay on one report line."""

    def test_same_structure_different_coefficient_is_paired(self):
        assignment = SourceAssignment("amat(var_u,var_u)", "u", 1, "2*a*b")
        slots = _source_monomial_slots(assignment_list(assignment), ["a*b"])
        self.assertEqual(len(slots), 1)
        _, source, generated = slots[0]
        self.assertEqual(source, "2*a*b")
        self.assertEqual(generated, "a*b")


def assignment_list(assignment):
    return [assignment]


if __name__ == "__main__":
    unittest.main()


class MatchingOrderTest(unittest.TestCase):
    """Exact matches must be resolved before the relaxed ones."""

    def test_exact_match_is_not_stolen_by_a_relaxed_one(self):
        # ``2*a*b`` has no exact partner and would, in a single pass, consume
        # ``a*b`` through the coefficient fallback, leaving the source line
        # that matches ``a*b`` exactly with an empty cell.
        assignment = SourceAssignment(
            "amat(var_u,var_u)", "u", 1, "2*a*b + a*b",
        )
        slots = _source_monomial_slots([assignment], ["a*b", "3*a*b"])
        rendered = [(source, generated) for _, source, generated in slots]
        self.assertEqual(len(rendered), 2)
        self.assertIn(("a*b", "a*b"), rendered)
        self.assertIn(("2*a*b", "3*a*b"), rendered)
        self.assertTrue(all(source and generated for source, generated in rendered))
