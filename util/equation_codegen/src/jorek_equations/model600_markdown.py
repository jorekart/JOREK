"""Meld-friendly Markdown exports of model-600 source and generated terms."""

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
import re

import sympy as sp

from .equations import LinearizedEquation
from .linearization import resolve_current, variation
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
    result = _normalize_model600_text(
        text,
        single_temperature=single_temperature,
        two_temperature=two_temperature,
    )
    return re.sub(r"\bBigR_x\b", "1", result, flags=re.I)


def _parsed_terms(text):
    """Return expanded monomials, or the original text when unsupported."""

    variants = []
    for single, two in ((False, False), (True, False), (False, True)):
        try:
            parsed = parse_fortran_expression(
                _normalized_text(text, single_temperature=single, two_temperature=two)
            )
            # ``normalize_fortran_text`` expands several element-routine
            # aliases (notably ``r0_x_hat``) after the textual normalization
            # above.  Those expansions can re-introduce ``BigR_x``.  The
            # report follows the comparison convention BigR_x=1, so apply it
            # to the parsed tree as well as to the input text.
            parsed = parsed.subs(sp.Symbol("BigR_x"), sp.Integer(1))
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


def _raw_evolution_rhs(equation, fields, timestep, zeta):
    """Build RHS without SymPy's global additive expansion."""

    history = sp.Add(
        *(variation(equation.A, value, direction=FieldRole.PREVIOUS_DELTA) for value in fields)
    )
    terms = [
        timestep * resolve_current(term)
        for term in sp.Add.make_args(equation.B)
    ]
    terms.extend(zeta * term for term in sp.Add.make_args(history))
    return terms


def _add_raw_rhs(pools, row, expression, *, previous_names=None):
    for term in expression:
        if term == 0:
            continue
        channel = _channel(term)
        lhs = "rhs_ij{}(var_{})".format(channel, row)
        pools[_clean_lhs(lhs)].append(fortran(term, previous_names=previous_names))


def _generated_pools(rows=ROWS, requested_lhs=None):
    """Generate implemented non-NEO terms only for the requested rows."""

    rows = set(rows)
    requested_lhs = set(requested_lhs or ())
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
        only_u_rhs = requested_lhs and requested_lhs <= {"rhs_ij(var_u)"}
        for with_tite in (False, True):
            equation = momentum_equation_2(
                with_TiTe=with_tite,
                include_diamagnetic=True,
                include_conservative=True,
                include_tgnum=True,
                include_neutrals=True,
                include_impurities=not only_u_rhs,
                # NEO AMAT is deliberately excluded until its symbolic
                # expansion is supported by the integrated checker.
                include_neo=False,
            )
            # The source routine contains both the one- and two-temperature
            # RHS branches.  Retain both raw expressions so terms involving
            # ``Ti0``/``Te0`` can be aligned as well as the single-T terms.
            _add_raw_rhs(
                pools,
                "u",
                _raw_evolution_rhs(equation, FIELDS, timestep, zeta),
                previous_names={"u": "delta_u", "rho": "delta_rho_g"},
            )
            if not only_u_rhs:
                result = equation.linearize(
                    fields=FIELDS, timestep=timestep, theta=theta, zeta=zeta
                )
                _add_linearized(
                    pools, "u", result,
                    previous_names={"u": "delta_u", "rho": "delta_rho_g"},
                    include_rhs=False,
                )
            if only_u_rhs and with_tite:
                break
    # For the momentum report retain repeated terms: the Fortran routine has
    # separate base/extension assignments whose identical monomials must be
    # available more than once for source-row alignment.  Other rows retain
    # the historical duplicate suppression.
    if rows == {"u"}:
        return {lhs: list(terms) for lhs, terms in pools.items()}
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


def _expanded_candidates(text):
    candidates = set()
    for terms in _parsed_terms(text):
        candidates.update(terms)
    return candidates


def _u_rhs_order_key(text):
    """Order factored generated momentum terms like the Fortran RHS."""

    lowered = text.lower()
    # Specific extension signatures must precede broad inertia/conservative
    # signatures because their expanded products contain the same factors.
    specific = (
        ("fact_conservative_u" in lowered and "vpar0" in lowered, 19),
        ("particle_source" in lowered or "source_pellet" in lowered, 18),
        ("fact_conservative_u" in lowered and "tgnum_u" in lowered, 10),
        ("dvisco_dt" in lowered, 14),
        ("visco_t" in lowered and "v_xx" in lowered and "tauic" in lowered, 15),
        ("tauic" in lowered and "u0_xy" in lowered, 13),
        ("tauic" in lowered and "r0_y" in lowered, 12),
        ("visco_t" in lowered and "tauic" in lowered, 15),
    )
    for condition, rank in specific:
        if condition:
            return rank
    groups = (
        ("r0_x", "r0_y", "u0_x**2"),       # inertia
        ("w0*r0",),                          # advection
        ("zj0_s", "zj0_t"),                 # magnetic bracket
        ("visco_fact_old",),
        ("visco_fact_new", "w0"),
        ("u0_xpp", "u0_ypp"),
        ("zj0_p",),
        ("v_s", "T0", "r0_t"),             # pressure
        ("visco_num_T",),
        ("tgnum_u", "r0*w0"),
        ("fact_conservative_u*tgnum_u",),
        ("tauIC", "w0_s"),                  # diamagnetic
        ("tauIC", "Pi0_y"),
        ("Pi0_xx", "Pi0_yy"),
        ("dvisco_dT", "W_dia"),
        ("W_dia", "visco_T"),
        ("delta_u",),
        ("Sion_T", "Srec_T"),
        ("particle_source", "source_"),
        ("fact_conservative_u",),
        ("aux_P_", "aux_divP"),
    )
    for index, needles in enumerate(groups):
        if all(needle.lower() in lowered for needle in needles):
            return index
    for index, needles in enumerate(groups):
        if any(needle.lower() in lowered for needle in needles):
            return index
    return len(groups)


def _expanded_ordered_monomials(text, *, two_temperature=False):
    variants = _parsed_terms(text)
    if not variants:
        return []
    terms = variants[-1] if two_temperature and "te0" not in text.lower() else variants[0]
    common_symbols = set.intersection(
        *(set(term.free_symbols) for term in terms)
    ) if terms else set()
    return sorted(
        terms,
        key=lambda term: _source_term_order(text, term, common_symbols),
    )


def _source_monomial_slots(assignments, generated_terms):
    """Expand both sides and align every generated monomial to source order."""

    source_monomials = []
    for assignment in assignments:
        for piece in _split_top_level(assignment.expression):
            for term in _expanded_ordered_monomials(piece):
                if term == 0:
                    continue
                source_monomials.append((assignment, term))

    generated_by_key = defaultdict(list)
    for text in generated_terms:
        for term in _expanded_ordered_monomials(text):
            # The residual history in the element routine uses the physical
            # background density ``r0``; only the velocity tangent uses the
            # corrected coefficient ``r0_corr``.  The DSL shares one A form,
            # so apply that source convention when matching RHS monomials.
            if assignments and assignments[0].lhs == "rhs_ij(var_u)":
                term = term.xreplace({sp.Symbol("r0_corr"): sp.Symbol("r0")})
            # JOREK stores W_dia_rho as a logarithmic density derivative
            # (the explicit trial-density factor is absorbed in the work
            # variable), whereas the DSL variation exposes that factor.
            if (
                assignments
                and assignments[0].lhs == "amat(var_u,var_rho)"
                and term.has(sp.Symbol("W_dia_rho"))
                and term.has(sp.Symbol("rho"))
            ):
                term = sp.cancel(term / sp.Symbol("rho"))
            if (
                assignments
                and assignments[0].lhs == "amat(var_u,var_ti)"
                and term.has(sp.Symbol("W_dia_Ti"))
                and term.has(sp.Symbol("Ti"))
            ):
                term = sp.cancel(term / sp.Symbol("Ti"))
            if (
                assignments
                and assignments[0].lhs == "amat(var_u,var_t)"
                and term.has(sp.Symbol("W_dia_T"))
                and term.has(sp.Symbol("T"))
            ):
                term = sp.cancel(term / sp.Symbol("T"))
            if (
                assignments
                and assignments[0].lhs == "amat(var_u,var_rhoimp)"
                and term.has(sp.Symbol("alpha_e_bis"))
                and term.has(sp.Symbol("Te0"))
            ):
                # alpha_e_bis is supplied as d(alpha_e*Te)/dTe, so the work
                # variable already contains the explicit background Te0.
                term = sp.cancel(term / sp.Symbol("Te0"))
            generated_by_key[_canonical_monomial(term)].append(term)
            # SymPy combines the two identical active/background TGNUM
            # variations into a single 1/2 monomial.  The Fortran source keeps
            # them as two 1/4 assignments, so expose the split representation
            # to the row matcher as well.
            if (
                term.has(sp.Symbol("tgnum_u"))
                and term.is_Mul
                and abs(term.as_coeff_Mul()[0]) >= sp.Rational(1, 2)
            ):
                half = term / 2
                generated_by_key[_canonical_monomial(half)].extend((half, half))

    slots = []
    for assignment, source_term in source_monomials:
        source_display = _canonical_monomial(source_term)
        match_term = source_term
        if assignment.lhs == "amat(var_u,var_t)":
            # In the one-temperature source branch Ti0 is the common T0 and
            # W_dia_Ti is the derivative of W_dia with respect to that common
            # temperature.
            match_term = source_term.xreplace({
                sp.Symbol("Ti0_x"): sp.Symbol("T0_x"),
                sp.Symbol("Ti0_y"): sp.Symbol("T0_y"),
                sp.Symbol("W_dia_Ti"): sp.Symbol("W_dia_T"),
            })
        key = _canonical_monomial(match_term)
        matches = generated_by_key[key]
        generated = _canonical_monomial(matches.pop(0)) if matches else ""
        if generated and assignment.lhs == "amat(var_u,var_t)":
            # Render the accepted one-temperature aliases with the source
            # spelling so Meld shows a truly aligned report.
            generated = source_display
        slots.append((assignment, source_display, generated))
    return slots


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
    positional_consumed = defaultdict(int)
    # For a model-600 momentum report, use the same monomial alignment for
    # every RHS/AMAT block.  Each assignment is matched independently because
    # the source contains separate columns and toroidal channels.  This keeps
    # the generated line in the same row as its source monomial and leaves a
    # genuinely empty line when that contribution is not implemented.
    if assignments and all(assignment.row == "u" for assignment in assignments):
        by_lhs = defaultdict(list)
        for assignment in assignments:
            by_lhs[assignment.lhs].append(assignment)
        slots = []
        for lhs, lhs_assignments in by_lhs.items():
            slots.extend(
                _source_monomial_slots(lhs_assignments, generated.get(lhs, ()))
            )
        return slots
    # Parsing generated text is expensive, especially for the expanded u
    # tangent.  Build the canonical lookup once instead of once per source
    # term.
    canonical = {
        lhs: [_canonical_candidates(text) for text in terms]
        for lhs, terms in generated.items()
    }
    expanded = {
        lhs: [_expanded_candidates(text) for text in terms]
        for lhs, terms in generated.items()
    }
    for assignment in assignments:
        if assignment.lhs == "rhs_ij(var_u)":
            # The weak momentum RHS is intentionally kept factored.  Align its
            # outer Fortran terms positionally; expanding them into monomials
            # creates hundreds of artificial generated-only slots.
            source_pieces = _split_top_level(assignment.expression)
            all_generated_terms = generated.get(assignment.lhs, ())
            start = positional_consumed[assignment.lhs]
            generated_terms = sorted(
                all_generated_terms[start:], key=_u_rhs_order_key
            )
            positional_consumed[assignment.lhs] += len(generated_terms)
            generated_index = 0
            for source in source_pieces:
                take = 2 if "delta_u_x" in source and "delta_u_y" in source else 1
                selected = generated_terms[generated_index:generated_index + take]
                generated_index += take
                generated_term = " + ".join(selected)
                slots.append((assignment, source, generated_term))
            used[assignment.lhs].update(
                range(start, start + len(generated_terms))
            )
            continue
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
        if lhs == "rhs_ij(var_u)":
            # Focused u-RHS reports are source-slot aligned; do not append
            # SymPy-expanded leftovers as thousands of artificial lines.
            continue
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


def export_model600_markdown(
    source_path, source_output, generated_output, equations=None, assignments=None
):
    """Write aligned source/generated reports and return their paths.

    ``equations`` may contain any of ``psi``, ``u``, ``zj`` and ``w``.  If it
    is omitted, all four rows are exported.
    """

    rows = tuple(equations or ROWS)
    invalid = set(rows) - set(ROWS)
    if invalid:
        raise ValueError("Unknown model-600 equation(s): {}".format(", ".join(sorted(invalid))))
    assignment_filter = assignments
    assignments = [
        assignment for assignment in _read_source_assignments(source_path)
        if assignment.row in rows
    ]
    requested = {_clean_lhs(name) for name in assignment_filter} if assignment_filter else None
    generated = _generated_pools(rows, requested_lhs=requested)
    if assignment_filter is not None:
        assignments = [assignment for assignment in assignments if assignment.lhs in requested]
        generated = {
            lhs: terms for lhs, terms in generated.items() if lhs in requested
        }
    slots = _source_slots(assignments, generated)
    source_output = Path(source_output)
    generated_output = Path(generated_output)
    source_output.parent.mkdir(parents=True, exist_ok=True)
    generated_output.parent.mkdir(parents=True, exist_ok=True)
    source_output.write_text(_render_report(slots, "source", rows), encoding="utf-8")
    generated_output.write_text(_render_report(slots, "generated", rows), encoding="utf-8")
    return source_output, generated_output
