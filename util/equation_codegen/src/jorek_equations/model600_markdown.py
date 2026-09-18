"""Meld-friendly Markdown exports of model-600 source and generated terms."""

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
import re

import sympy as sp

from .equations import LinearizedEquation
from .fortran import fortran
from .model600 import (
    FIELDS,
    T,
    Te,
    current_constraint_equation_zj,
    induction_equation_1,
    j,
    momentum_equation_2,
    omega,
    psi,
    rho,
    u,
    vorticity_constraint_equation_w,
)
from .model600_compare import _normalize_model600_text
from .operators import SpatialDerivative, expand_derivatives, phi
from .source_compare import parse_fortran_expression
from .symbols import FieldRole, FieldValue, TestFunction, coefficient


ROWS = ("psi", "u", "zj", "w")
FIELD_NAMES = {
    field: name
    for field, name in zip(FIELDS, ("psi", "u", "zj", "w", "rho", "T", "vpar", "Ti", "Te", "rhon", "rhoimp"))
}
ASSIGNMENT_RE = re.compile(
    r"(?P<lhs>(?:rhs_ij(?:_k)?|amat(?:_n|_k|_kn|_nn)?)\s*\(\s*"
    r"var_(?P<row>psi|u|zj|w)\b[^=]*?\))\s*=\s*(?P<rhs>.*)$",
    re.I,
)


@dataclass(frozen=True)
class SourceAssignment:
    lhs: str
    row: str
    line: int
    expression: str


def _clean_lhs(lhs):
    return re.sub(r"\s+", "", lhs).lower()


def _read_source_assignments(path):
    """Read every relevant assignment, retaining repetitions and file order."""

    lines = Path(path).read_text(encoding="utf-8").splitlines()
    result = []
    index = 0
    while index < len(lines):
        code = lines[index].split("!", 1)[0]
        match = ASSIGNMENT_RE.search(code)
        if match is None:
            index += 1
            continue
        lhs = _clean_lhs(match.group("lhs"))
        pieces = []
        current = match.group("rhs")
        while True:
            current = current.strip()
            continued = current.endswith("&")
            pieces.append(current[:-1].strip() if continued else current)
            if not continued:
                break
            index += 1
            if index >= len(lines):
                raise ValueError("Unterminated assignment at line {}".format(index + 1))
            # JOREK often places a blank/comment line between continued
            # terms.  It is not the end of the Fortran statement.
            current = lines[index].split("!", 1)[0].strip()
            while not current:
                index += 1
                if index >= len(lines):
                    raise ValueError("Unterminated assignment at line {}".format(index + 1))
                current = lines[index].split("!", 1)[0].strip()
            if current.startswith("&"):
                current = current[1:]
        expression = " ".join(piece for piece in pieces if piece)
        # Repeated extension assignments have the form ``a = a + extension``.
        expression = re.sub(
            r"^{}\s*\+\s*".format(re.escape(lhs)), "", _clean_assignment_refs(expression),
            flags=re.I,
        )
        result.append(
            SourceAssignment(lhs, match.group("row").lower(), index + 2 - len(pieces), expression)
        )
        index += 1
    return result


def _clean_assignment_refs(expression):
    """Canonicalize only assignment references, leaving term spelling intact."""

    return re.sub(
        r"(?:rhs_ij(?:_k)?|amat(?:_n|_k|_kn|_nn)?)\s*\([^)]*\)",
        lambda match: _clean_lhs(match.group(0)),
        expression,
        flags=re.I,
    )


def _split_top_level(expression):
    """Split an expression at ordered, outermost additive signs."""

    terms = []
    start = 0
    depth = 0
    for index, char in enumerate(expression):
        if char == "(":
            depth += 1
        elif char == ")":
            depth = max(0, depth - 1)
        elif char in "+-" and depth == 0 and index > start:
            previous = expression[index - 1]
            if previous not in "eEdD":
                terms.append(expression[start:index].strip())
                start = index
    tail = expression[start:].strip()
    if tail:
        terms.append(tail)
    return terms


def _normalized_text(text, *, single_temperature=False, two_temperature=False):
    text = re.sub(
        r"factor\(\s*var_[a-z0-9_]+\s*,\s*[0-9]+\s*\)", "1", text,
        flags=re.I,
    )
    text = re.sub(r"\bBigR_x\b", "1", text, flags=re.I)
    return _normalize_model600_text(
        text,
        single_temperature=single_temperature,
        two_temperature=two_temperature,
    )


def _parsed_terms(text):
    """Return expanded monomials, or the original text when unsupported."""

    variants = []
    for single, two in ((False, False), (True, False), (False, True)):
        try:
            parsed = parse_fortran_expression(
                _normalized_text(text, single_temperature=single, two_temperature=two)
            )
        except Exception:
            continue
        terms = tuple(sp.Add.make_args(sp.expand(parsed)))
        if terms not in variants:
            variants.append(terms)
    return variants


def _channel(term):
    test_order = 0
    trial_order = 0
    for derivative in term.atoms(SpatialDerivative):
        if derivative.coordinate != phi:
            continue
        if derivative.expression.has(TestFunction):
            test_order += 1
        if any(
            value.role is FieldRole.TRIAL
            for value in derivative.expression.atoms(FieldValue)
        ):
            trial_order += 1
    if test_order == 0 and trial_order == 0:
        return ""
    if test_order == 0 and trial_order == 1:
        return "_n"
    if test_order == 1 and trial_order == 0:
        return "_k"
    if test_order == 1 and trial_order == 1:
        return "_kn"
    if test_order == 0 and trial_order == 2:
        return "_nn"
    raise ValueError("Unsupported toroidal channel in {}".format(term))


def _add_linearized(
    pools, row, linearized: LinearizedEquation, *, previous_names=None,
    include_rhs=True, include_amat=True, include_fields=None,
):
    if include_rhs:
        for term in sp.Add.make_args(sp.expand(expand_derivatives(linearized.rhs))):
            if term == 0:
                continue
            channel = _channel(term)
            lhs = "rhs_ij{}(var_{})".format(channel, row)
            pools[_clean_lhs(lhs)].append(fortran(term, previous_names=previous_names))
    if include_amat:
        expressions = linearized.amat.items()
    else:
        expressions = ()
    if include_fields is not None:
        include_fields = set(include_fields)
        expressions = (
            (field, expression)
            for field, expression in expressions
            if field in include_fields
        )
    for field, expression in expressions:
        column = FIELD_NAMES[field]
        for term in sp.Add.make_args(sp.expand(expand_derivatives(expression))):
            if term == 0:
                continue
            channel = _channel(term)
            lhs = "amat{}(var_{},var_{})".format(channel, row, column)
            pools[_clean_lhs(lhs)].append(fortran(term, previous_names=previous_names))


def _generated_pools(rows=ROWS):
    """Generate implemented non-NEO terms only for the requested rows."""

    rows = set(rows)
    pools = defaultdict(list)
    timestep = coefficient("tstep")
    theta = coefficient("theta")
    zeta = coefficient("zeta")
    if "psi" in rows:
        for with_tite in (False,):
            equation = induction_equation_1(with_TiTe=with_tite)
            result = equation.linearize(
                fields=FIELDS, timestep=timestep, theta=theta, zeta=zeta
            )
            _add_linearized(
                pools, "psi", result,
                previous_names={"psi": "delta_g(mp,var_psi,ms,mt)"},
                include_fields=tuple(field for field in FIELDS if field not in (psi, rho, Te)),
            )
        # The physical RHS is a single expression; only the AMAT has separate
        # one- and two-temperature source branches.
        equation = induction_equation_1(with_TiTe=True)
        result = equation.linearize(
            fields=FIELDS, timestep=timestep, theta=theta, zeta=zeta
        )
        _add_linearized(
            pools, "psi", result,
            previous_names={"psi": "delta_g(mp,var_psi,ms,mt)"},
            include_rhs=False,
            include_fields=(psi, rho, Te),
        )
    if "zj" in rows:
        _add_linearized(
            pools, "zj", current_constraint_equation_zj().linearize(fields=FIELDS)
        )
    if "w" in rows:
        _add_linearized(
            pools, "w", vorticity_constraint_equation_w().linearize(fields=FIELDS)
        )
    if "u" in rows:
        for with_tite in (False, True):
            equation = momentum_equation_2(
                with_TiTe=with_tite,
                include_diamagnetic=True,
                include_conservative=True,
                include_tgnum=True,
                include_neutrals=True,
                include_impurities=True,
                # NEO AMAT is deliberately excluded until its symbolic
                # expansion is supported by the integrated checker.
                include_neo=False,
            )
            result = equation.linearize(
                fields=FIELDS, timestep=timestep, theta=theta, zeta=zeta
            )
            _add_linearized(pools, "u", result, previous_names={"u": "delta_u"})
    # The two thermal variants share many terms.  A source term needs only one
    # corresponding generated slot, so discard exact textual duplicates.
    return {
        lhs: list(dict.fromkeys(terms))
        for lhs, terms in pools.items()
    }


def _canonical_candidates(text):
    candidates = set()
    for terms in _parsed_terms(text):
        if len(terms) == 1:
            candidates.add(terms[0])
    return candidates


def _source_term_order(source_piece, term, common_symbols=()):
    """Estimate textual order for expanded terms inside one source term."""

    text = source_piece.lower()
    positions = []
    common_symbols = set(common_symbols)
    for symbol in term.free_symbols:
        if symbol in common_symbols:
            continue
        name = str(symbol).lower()
        # A symbol appearing once is usually the factor that distinguishes this
        # expanded monomial (zj0/current_source/Jb in the induction example).
        if len(re.findall(r"\b{}\b".format(re.escape(name)), text)) == 1:
            positions.append(text.find(name))
    return min(positions) if positions else len(text)


def _factor_order(factor):
    """Use a stable JOREK-oriented order for multiplicative factors."""

    base = factor
    while getattr(base, "is_Pow", False):
        base = base.base
    name = str(base)
    prefixes = (
        "F0", "BigR", "Z", "v", "u", "psi", "zj", "w", "rho",
        "T", "Ti", "Te", "xjac", "theta", "tstep", "zeta",
    )
    prefix_rank = next(
        (index for index, prefix in enumerate(prefixes) if name.startswith(prefix)),
        len(prefixes),
    )
    return prefix_rank, name, str(factor)


def _canonical_monomial(term):
    """Render one monomial with the same factor order in both reports."""

    coefficient = sp.S.One
    factors = []
    for factor in sp.Mul.make_args(term):
        if factor.is_Number:
            coefficient *= factor
        else:
            factors.append(factor)
    factors.sort(key=_factor_order)
    numerator = [factor for factor in factors if not (factor.is_Pow and factor.exp.is_negative)]
    denominator = [factor.base ** (-factor.exp) for factor in factors if factor.is_Pow and factor.exp.is_negative]
    pieces = []
    if coefficient not in (1, -1) or (not numerator and not denominator):
        pieces.append(sp.sstr(coefficient))
    pieces.extend(sp.sstr(factor) for factor in numerator)
    result = "*".join(pieces) if pieces else "1"
    for factor in denominator:
        rendered = sp.sstr(factor)
        if not factor.is_Atom:
            rendered = "({})".format(rendered)
        result += "/{}".format(rendered)
    if coefficient == -1 and (numerator or denominator):
        result = "-{}".format(result)
    return result


def _canonical_display(text, *, two_temperature=False):
    """Normalize multiplication ordering without changing term grouping."""

    variants = _parsed_terms(text)
    if not variants:
        return text
    if two_temperature:
        terms = variants[0] if "te0" in text.lower() else variants[-1]
    else:
        terms = variants[0]
    if len(terms) == 1 and terms[0] == 0:
        return ""
    return " + ".join(_canonical_monomial(term) for term in terms)


def _source_slots(assignments, generated):
    """Align generated monomials to source monomials in source order."""

    slots = []
    used = defaultdict(set)
    # Parsing generated text is expensive, especially for the expanded u
    # tangent.  Build the canonical lookup once instead of once per source
    # term.
    canonical = {
        lhs: [_canonical_candidates(text) for text in terms]
        for lhs, terms in generated.items()
    }
    for assignment in assignments:
        for source_piece in _split_top_level(assignment.expression):
            variants = _parsed_terms(source_piece)
            if not variants:
                slots.append((assignment, source_piece, ""))
                continue
            # Keep one visual slot for each outer Fortran term.  A source term
            # such as ``eta*(zj0-current_source-Jb)`` expands into multiple
            # SymPy monomials; all matching generated monomials are folded back
            # into this one source-ordered slot.
            # The generic psi/rho blocks in model-600 use the two-temperature
            # ``Pe0`` convention, while the explicit ``var_T`` block uses the
            # single-temperature branch.  Select the matching source alias
            # expansion before ordering its monomials.
            use_two_temperature = (
                assignment.lhs in {
                    "amat(var_psi,var_psi)",
                    "amat(var_psi,var_rho)",
                    "amat_n(var_psi,var_rho)",
                }
                or "var_te" in assignment.lhs
            )
            if use_two_temperature:
                # Explicit Te0 terms are already in the desired branch; only
                # the abstract Pe0 aliases need the final two-temperature
                # expansion.
                variant_index = (
                    0 if "te0" in source_piece.lower() else len(variants) - 1
                )
            else:
                variant_index = 0
            source_terms_raw = variants[variant_index]
            common_symbols = set.intersection(
                *(set(term.free_symbols) for term in source_terms_raw)
            ) if source_terms_raw else set()
            source_terms = sorted(
                source_terms_raw,
                key=lambda term: _source_term_order(
                    source_piece, term, common_symbols
                ),
            )
            matched_terms = []
            for source_term in source_terms:
                matched = None
                for index, generated_text in enumerate(generated.get(assignment.lhs, ())):
                    if index in used[assignment.lhs]:
                        continue
                    if source_term in canonical[assignment.lhs][index]:
                        used[assignment.lhs].add(index)
                        matched = generated_text
                        break
                if matched is not None:
                    matched_terms.append(matched)
            slots.append((assignment, source_piece, " + ".join(matched_terms)))
    # Generated-only terms follow the source slots in their assignment block.
    for lhs, terms in generated.items():
        row_match = re.search(r"\(var_(psi|u|zj|w)", lhs)
        if row_match is None:
            continue
        synthetic = SourceAssignment(lhs, row_match.group(1), 0, "")
        for index, term in enumerate(terms):
            if index not in used[lhs]:
                slots.append((synthetic, "", term))
    return slots


def _render_report(slots, side, rows=ROWS):
    lines = [
        "# Model 600 equation terms",
        "",
        "> Rows currently exported: {}.".format(
            ", ".join("`{}`".format(row) for row in rows)
        ),
        "> NEO AMAT terms are not generated yet; their generated slots are blank.",
        "",
    ]
    for row in rows:
        lines.extend(("## var_{}".format(row), ""))
        row_slots = [slot for slot in slots if slot[0].row == row]
        lhs_order = list(dict.fromkeys(slot[0].lhs for slot in row_slots))
        for lhs in lhs_order:
            lines.extend(("### `{}`".format(lhs), ""))
            for assignment, source, generated in row_slots:
                if assignment.lhs != lhs:
                    continue
                value = source if side == "source" else generated
                value = " ".join(value.split())
                use_two_temperature = (
                    assignment.lhs in {
                        "amat(var_psi,var_psi)",
                        "amat(var_psi,var_rho)",
                        "amat_n(var_psi,var_rho)",
                    }
                    or "var_te" in assignment.lhs
                )
                value = _canonical_display(
                    value,
                    two_temperature=use_two_temperature and side == "source",
                )
                # Keep exactly one physical line per aligned term.  An
                # unavailable term is intentionally an actually empty line,
                # which makes its position immediately visible in Meld.
                lines.append(value)
            lines.append("")
    return "\n".join(lines) + "\n"


def export_model600_markdown(source_path, source_output, generated_output, equations=None):
    """Write aligned source/generated reports and return their paths.

    ``equations`` may contain any of ``psi``, ``u``, ``zj`` and ``w``.  If it
    is omitted, all four rows are exported.
    """

    rows = tuple(equations or ROWS)
    invalid = set(rows) - set(ROWS)
    if invalid:
        raise ValueError("Unknown model-600 equation(s): {}".format(", ".join(sorted(invalid))))
    assignments = [
        assignment for assignment in _read_source_assignments(source_path)
        if assignment.row in rows
    ]
    generated = _generated_pools(rows)
    slots = _source_slots(assignments, generated)
    source_output = Path(source_output)
    generated_output = Path(generated_output)
    source_output.parent.mkdir(parents=True, exist_ok=True)
    generated_output.parent.mkdir(parents=True, exist_ok=True)
    source_output.write_text(_render_report(slots, "source", rows), encoding="utf-8")
    generated_output.write_text(_render_report(slots, "generated", rows), encoding="utf-8")
    return source_output, generated_output
