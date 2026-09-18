"""Compare generated symbolic terms with assignments in JOREK Fortran."""

from pathlib import Path
import re
from collections import Counter
from typing import Dict, Iterable, Mapping

import sympy as sp
from sympy.parsing.sympy_parser import parse_expr

from .channels import split_toroidal_channels
from .fortran import fortran
from .model199 import (
    FIELDS,
    T,
    current_constraint_equation_3,
    density_equation_5,
    induction_equation_1,
    j,
    momentum_equation_2,
    psi,
    rho,
    u,
    vorticity_constraint_equation_4,
    temperature_equation_6,
    xjac,
)
from .operators import expand_brackets, lower_brackets_to_element
from .symbols import coefficient


EQUATION_1_ASSIGNMENTS = (
    "rhs_ij_1",
    "amat_11",
    "amat_12",
    "amat_12_n",
    "amat_13",
    "amat_16",
)

EQUATION_2_ASSIGNMENTS = (
    "rhs_ij_2", "amat_21", "amat_22", "amat_23", "amat_23_n",
    "amat_24", "amat_25", "amat_26",
)
EQUATION_3_ASSIGNMENTS = ("rhs_ij_3", "amat_31", "amat_33")
EQUATION_4_ASSIGNMENTS = ("rhs_ij_4", "amat_42", "amat_44")
EQUATION_5_ASSIGNMENTS = (
    "rhs_ij_5", "rhs_ij_5_k", "amat_51", "amat_51_k", "amat_52",
    "amat_55", "amat_55_k", "amat_55_n", "amat_55_kn",
)
EQUATION_6_ASSIGNMENTS = (
    "rhs_ij_6", "rhs_ij_6_k", "amat_61", "amat_61_k", "amat_62",
    "amat_63", "amat_66", "amat_66_k", "amat_66_n", "amat_66_kn",
)


class FortranComparisonError(AssertionError):
    """Generated and source Fortran expressions differ."""


def format_expression_terms(expression) -> str:
    """Format an expression as one additive term per line.

    The comparison itself remains symbolic; this is only a diagnostic
    representation.  Expanding the outer addition makes it much easier to
    spot a missing sign or coefficient in the generated weak-form terms.
    """

    expanded = sp.expand(expression)
    if expanded == 0:
        return "    0"
    terms = sp.Add.make_args(expanded)
    lines = []
    for term in terms:
        if term.could_extract_minus_sign():
            lines.append("    - {}".format(str(-term)))
        else:
            lines.append("    + {}".format(str(term)))
    return "\n".join(lines)


def format_assignment(name: str, expression, label: str = "") -> str:
    """Return a readable, multiline assignment diagnostic."""

    heading = "{}{}:".format(label + " " if label else "", name)
    return heading + "\n" + format_expression_terms(expression)


def format_term_difference(source, generated) -> str:
    """Show unmatched additive terms without subtracting the expressions."""

    source_terms = Counter(sp.Add.make_args(sp.expand(source)))
    generated_terms = Counter(sp.Add.make_args(sp.expand(generated)))
    source_only = list((source_terms - generated_terms).elements())
    generated_only = list((generated_terms - source_terms).elements())
    lines = []
    if source_only:
        lines.append("    source-only:")
        lines.extend("      {}".format(_term_text(term)) for term in source_only)
    if generated_only:
        lines.append("    generated-only:")
        lines.extend("      {}".format(_term_text(term)) for term in generated_only)
    return "\n".join(lines) if lines else "    (no unmatched terms)"


def _term_text(term) -> str:
    """Render one term while retaining its sign."""

    if term.could_extract_minus_sign():
        return "- {}".format(str(-term))
    return "+ {}".format(str(term))


def extract_fortran_assignments(path, names: Iterable[str]) -> Dict[str, str]:
    """Extract continued scalar assignments from a Fortran source file."""

    requested = tuple(names)
    patterns = {
        name: re.compile(r"^\s*{}\s*=\s*(.*)$".format(re.escape(name)), re.I)
        for name in requested
    }
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    found = {}
    index = 0
    while index < len(lines):
        uncommented = lines[index].split("!", 1)[0]
        matched_name = None
        matched_rhs = None
        for name, pattern in patterns.items():
            match = pattern.match(uncommented)
            if match:
                matched_name = name
                matched_rhs = match.group(1)
                break
        if matched_name is None:
            index += 1
            continue

        pieces = []
        current = matched_rhs
        while True:
            current = current.strip()
            continued = current.endswith("&")
            pieces.append(current[:-1] if continued else current)
            if not continued:
                break
            index += 1
            if index >= len(lines):
                raise ValueError("Unterminated assignment for {}".format(matched_name))
            current = lines[index].split("!", 1)[0].lstrip()
            if current.startswith("&"):
                current = current[1:]
        found[matched_name] = " ".join(piece.strip() for piece in pieces)
        index += 1

    missing = [name for name in requested if name not in found]
    if missing:
        raise ValueError(
            "Assignments not found in {}: {}".format(path, ", ".join(missing))
        )
    return found


_D_LITERAL = re.compile(r"(?i)(?<![A-Za-z_])(\d+(?:\.\d*)?|\.\d+)d([+-]?\d+)")
_IDENTIFIER = re.compile(r"\b[A-Za-z_][A-Za-z0-9_]*\b")


def normalize_fortran_text(expression: str) -> str:
    """Normalize the small Fortran-expression subset used by the comparator."""

    result = expression.replace("&", " ")
    result = re.sub(r"\bBigr\b", "BigR", result, flags=re.I)
    result = re.sub(
        r"current_source\s*\(\s*ms\s*,\s*mt\s*\)",
        "current_source",
        result,
        flags=re.I,
    )
    for function_name in ("particle_source", "heat_source"):
        result = re.sub(
            r"{}\s*\(\s*ms\s*,\s*mt\s*\)".format(function_name),
            function_name,
            result,
            flags=re.I,
        )
    result = re.sub(
        r"delta_g\s*\(\s*mp\s*,\s*1\s*,\s*ms\s*,\s*mt\s*\)",
        "delta_g_1",
        result,
        flags=re.I,
    )
    result = re.sub(
        r"delta_g\s*\(\s*mp\s*,\s*([56])\s*,\s*ms\s*,\s*mt\s*\)",
        r"delta_g_\1",
        result,
        flags=re.I,
    )
    aliases = {
        "D_par_local": "D_par",
        "vv2": "(BigR**2*(u0_x**2+u0_y**2))",
        "r0_hat": "(BigR**2*r0)",
        "r0_x_hat": "(2*BigR*BigR_x*r0+BigR**2*r0_x)",
        "r0_y_hat": "(BigR**2*r0_y)",
        "rho_hat": "(BigR**2*rho)",
        "rho_x_hat": "(2*BigR*BigR_x*rho+BigR**2*rho_x)",
        "rho_y_hat": "(BigR**2*rho_y)",
        "P0_s": "(r0_s*T0+r0*T0_s)",
        "P0_t": "(r0_t*T0+r0*T0_t)",
        "BB2": "((F0**2+ps0_x**2+ps0_y**2)/BigR**2)",
        "BB2_psi": "(2*(psi_x*ps0_x+psi_y*ps0_y)/BigR**2)",
        "Bgrad_rho_star": "((v_x*ps0_y-v_y*ps0_x)/BigR)",
        "Bgrad_rho_k_star": "(F0*v_p/BigR**2)",
        "Bgrad_rho": "((F0*r0_p/BigR+r0_x*ps0_y-r0_y*ps0_x)/BigR)",
        "Bgrad_rho_star_psi": "((v_x*psi_y-v_y*psi_x)/BigR)",
        "Bgrad_rho_psi": "((r0_x*psi_y-r0_y*psi_x)/BigR)",
        "Bgrad_rho_rho": "((rho_x*ps0_y-rho_y*ps0_x)/BigR)",
        "Bgrad_rho_rho_n": "(F0*rho_p/BigR**2)",
        "Bgrad_T_star": "((v_x*ps0_y-v_y*ps0_x)/BigR)",
        "Bgrad_T_k_star": "(F0*v_p/BigR**2)",
        "Bgrad_T": "((F0*T0_p/BigR+T0_x*ps0_y-T0_y*ps0_x)/BigR)",
        "Bgrad_T_star_psi": "((v_x*psi_y-v_y*psi_x)/BigR)",
        "Bgrad_T_psi": "((T0_x*psi_y-T0_y*psi_x)/BigR)",
        "Bgrad_T_T": "((T_x*ps0_y-T_y*ps0_x)/BigR)",
        "Bgrad_T_T_n": "(F0*T_p/BigR**2)",
    }
    for name, replacement in aliases.items():
        result = re.sub(r"\b{}\b".format(name), replacement, result, flags=re.I)
    result = _D_LITERAL.sub(lambda match: match.group(1) + "e" + match.group(2), result)
    return " ".join(result.split())


def parse_fortran_expression(expression: str):
    """Parse a normalized scalar Fortran expression into an exact SymPy tree."""

    normalized = normalize_fortran_text(expression)
    names = set(_IDENTIFIER.findall(normalized))
    local_dict = {name: sp.Symbol(name) for name in names}
    parsed = parse_expr(normalized, local_dict=local_dict, evaluate=True)
    return sp.nsimplify(parsed)


def generated_model199_equation1_assignments() -> Dict[str, sp.Expr]:
    """Generate equation-1 terms in the same local convention as the source."""

    equation = induction_equation_1()
    result = equation.linearize(
        fields=FIELDS,
        timestep=coefficient("tstep"),
        theta=coefficient("theta"),
        zeta=coefficient("zeta"),
    )

    rhs_channels = split_toroidal_channels(result.rhs)
    psi_channels = split_toroidal_channels(result.amat[psi])
    u_channels = split_toroidal_channels(result.amat[u])
    j_channels = split_toroidal_channels(result.amat[j])
    T_channels = split_toroidal_channels(result.amat[T])

    unsupported = {
        "rhs": (rhs_channels.n, rhs_channels.k, rhs_channels.kn),
        "psi": (psi_channels.n, psi_channels.k, psi_channels.kn),
        "u": (u_channels.k, u_channels.kn),
        "j": (j_channels.n, j_channels.k, j_channels.kn),
        "T": (T_channels.n, T_channels.k, T_channels.kn),
    }
    for label, values in unsupported.items():
        if any(value != 0 for value in values):
            raise ValueError(
                "Unexpected toroidal channel in equation-1 {} block".format(label)
            )

    generated = {
        "rhs_ij_1": rhs_channels.p,
        "amat_11": psi_channels.p,
        "amat_12": u_channels.p,
        "amat_12_n": u_channels.n,
        "amat_13": j_channels.p,
        "amat_16": T_channels.p,
    }
    return {
        name: lower_brackets_to_element(expression, xjac)
        for name, expression in generated.items()
    }


def generated_assignment_text() -> Dict[str, str]:
    """Render generated terms with the identifiers used in model 199."""

    return {
        name: fortran(
            expression,
            previous_names={"psi": "delta_g(mp,1,ms,mt)"},
        )
        for name, expression in generated_model199_equation1_assignments().items()
    }


def generated_partial_model199_assignments() -> Dict[str, str]:
    """Render generated assignments for the currently implemented equations."""

    timestep = coefficient("tstep")
    theta = coefficient("theta")
    zeta = coefficient("zeta")
    linearized = momentum_equation_2().linearize(
        fields=FIELDS, timestep=timestep, theta=theta, zeta=zeta
    )
    channels = split_toroidal_channels(linearized.rhs)
    result_map = {"rhs_ij_2": channels.p}
    for name, field, channel in (
        ("amat_21", psi, "p"), ("amat_22", u, "p"),
        ("amat_23", j, "p"), ("amat_23_n", j, "n"),
        ("amat_24", FIELDS[3], "p"), ("amat_25", FIELDS[4], "p"),
        ("amat_26", T, "p"),
    ):
        result_map[name] = getattr(split_toroidal_channels(linearized.amat[field]), channel)
    for equation, rhs_name, block_fields in (
        (current_constraint_equation_3(), "rhs_ij_3", (("amat_31", psi), ("amat_33", j))),
        (vorticity_constraint_equation_4(), "rhs_ij_4", (("amat_42", u), ("amat_44", FIELDS[3]))),
    ):
        linearized = equation.linearize(fields=FIELDS)
        result_map[rhs_name] = linearized.rhs
        for assignment, field in block_fields:
            result_map[assignment] = linearized.amat[field]
    return {
        name: fortran(value, previous_names={"u": "delta_u"})
        for name, value in result_map.items()
    }


def generated_model199_transport_assignments() -> Dict[str, str]:
    """Render generated RHS/AMAT assignments for equations 5 and 6."""

    timestep = coefficient("tstep")
    theta = coefficient("theta")
    zeta = coefficient("zeta")
    output = {}
    for equation, fields, names in (
        (density_equation_5(), FIELDS, EQUATION_5_ASSIGNMENTS),
        (temperature_equation_6(), FIELDS, EQUATION_6_ASSIGNMENTS),
    ):
        linearized = equation.linearize(
            fields=fields, timestep=timestep, theta=theta, zeta=zeta
        )
        rhs = split_toroidal_channels(linearized.rhs)
        output[names[0]] = rhs.p
        output[names[1]] = rhs.k
        row_blocks = (("51", psi), ("52", u), ("55", rho)) if equation.name.endswith("density") else (("61", psi), ("62", u), ("63", j), ("66", T))
        for block, field in row_blocks:
            blocks = split_toroidal_channels(linearized.amat[field])
            output["amat_{}".format(block)] = blocks.p
            for channel in ("k", "n", "kn"):
                if blocks.__getattribute__(channel) != 0:
                    output["amat_{}_{}".format(block, channel)] = getattr(blocks, channel)
    previous_names = {
        "rho": "delta_g(mp,5,ms,mt)",
        "T": "delta_g(mp,6,ms,mt)",
    }
    return {
        name: fortran(expand_brackets(value), previous_names=previous_names)
        for name, value in output.items()
    }


def compare_assignment_maps(
    source: Mapping[str, str],
    generated: Mapping[str, str],
    *,
    set_BigR_x_one: bool = True,
) -> None:
    """Compare assignment maps and raise with symbolic differences.

    ``set_BigR_x_one`` reflects the usual cylindrical-coordinate convention
    and removes element-geometry ``BigR_x`` factors from this comparison.
    """

    messages = []
    for name in source:
        if name not in generated:
            messages.append("{}: missing generated assignment".format(name))
            continue
        source_expr = parse_fortran_expression(source[name])
        generated_expr = parse_fortran_expression(generated[name])
        if set_BigR_x_one:
            bigr_x = sp.Symbol("BigR_x")
            source_expr = source_expr.subs(bigr_x, 1)
            generated_expr = generated_expr.subs(bigr_x, 1)
        difference = sp.simplify(generated_expr - source_expr)
        if difference != 0:
            messages.append(
                "{name}:\n"
                "{difference}".format(
                    name=name,
                    difference=format_term_difference(source_expr, generated_expr),
                )
            )
    extra = sorted(set(generated) - set(source))
    if extra:
        messages.append("unexpected generated assignments: {}".format(", ".join(extra)))
    if messages:
        raise FortranComparisonError(
            "Equation comparison failed:\n" + "\n".join(messages)
        )


def compare_model199_equation1(path, *, set_BigR_x_one: bool = True) -> None:
    """Compare generated equation 1 with model199/mod_elt_matrix_fft.f90."""

    source = extract_fortran_assignments(path, EQUATION_1_ASSIGNMENTS)
    compare_assignment_maps(
        source, generated_assignment_text(), set_BigR_x_one=set_BigR_x_one
    )
