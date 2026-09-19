import unittest

from jorek_equations.model600_markdown import (
    SourceAssignment,
    ROWS,
    _canonical_display,
    _render_report,
    _split_top_level,
)


class Model600MarkdownTest(unittest.TestCase):
    def test_top_level_split_keeps_parenthesized_sum_together(self):
        self.assertEqual(
            _split_top_level("a + b*(c-d) - e"),
            ["a", "+ b*(c-d)", "- e"],
        )

    def test_missing_generated_term_is_an_empty_aligned_line(self):
        assignment = SourceAssignment("rhs_ij(var_psi)", "psi", 1, "a")
        slots = [(assignment, "a", "")]
        source = _render_report(slots, "source").splitlines()
        generated = _render_report(slots, "generated").splitlines()
        self.assertEqual(len(source), len(generated))
        self.assertIn("a", source)
        self.assertEqual(generated[source.index("a")], "")

    def test_supported_equation_rows_are_explicit(self):
        self.assertEqual(ROWS, ("psi", "u", "zj", "w", "rho", "vpar", "rhoimp"))

    def test_factor_order_is_shared_between_source_and_generated_forms(self):
        self.assertEqual(
            _canonical_display("F0 / BigR * v * u_p * xjac * theta * tstep"),
            _canonical_display("F0*theta*tstep*xjac*u_p*v/BigR"),
        )


if __name__ == "__main__":
    unittest.main()
