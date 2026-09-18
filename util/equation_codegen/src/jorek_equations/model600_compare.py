"""Comparison support for model-600 equation 1 (the psi equation)."""

import re
from pathlib import Path

from .channels import split_toroidal_channels
from .fortran import fortran
from .model600 import (
    FIELDS,
    Te,
    current_constraint_equation_zj,
    induction_equation_1,
    j,
    omega,
    rho,
    T,
    u,
    psi,
    vpar,
    vorticity_constraint_equation_w,
)
from .source_compare import (
    FortranComparisonError,
    compare_assignment_maps,
    parse_fortran_expression,
)
from .symbols import coefficient


MODEL600_EQ1_ASSIGNMENTS = (
    "rhs_var_psi",
    "amat_var_psi_var_psi",
    "amat_var_psi_var_u",
    "amat_n_var_psi_var_u",
    "amat_var_psi_var_zj",
    "amat_var_psi_var_rho",
    "amat_n_var_psi_var_rho",
    "amat_var_psi_var_T",
    "amat_n_var_psi_var_T",
    "amat_var_psi_var_Te",
    "amat_n_var_psi_var_Te",
    "amat_var_psi_var_rhoimp",
)


def _assignment_key(kind, row, column=None):
    if column is None:
        return "{}_{}".format(kind, row)
    return "{}_{}_{}".format(kind, row, column)


def extract_model600_equation1(path):
    """Extract descriptive indexed assignments for ``var_psi``."""

    text = Path(path).read_text(encoding="utf-8").splitlines()
    patterns = {
        "rhs_var_psi": re.compile(r"^\s*rhs_ij\(\s*var_psi\s*\)\s*=\s*(.*)$", re.I),
        "amat_var_psi_var_psi": re.compile(r"^\s*amat\(\s*var_psi\s*,\s*var_psi\s*\)\s*=\s*(.*)$", re.I),
        "amat_var_psi_var_u": re.compile(r"^\s*amat\(\s*var_psi\s*,\s*var_u\s*\)\s*=\s*(.*)$", re.I),
        "amat_n_var_psi_var_u": re.compile(r"^\s*amat_n\(\s*var_psi\s*,\s*var_u\s*\)\s*=\s*(.*)$", re.I),
        "amat_var_psi_var_zj": re.compile(r"^\s*amat\(\s*var_psi\s*,\s*var_zj\s*\)\s*=\s*(.*)$", re.I),
        "amat_var_psi_var_rho": re.compile(r"^\s*amat\(\s*var_psi\s*,\s*var_rho\s*\)\s*=\s*(.*)$", re.I),
        "amat_n_var_psi_var_rho": re.compile(r"^\s*amat_n\(\s*var_psi\s*,\s*var_rho\s*\)\s*=\s*(.*)$", re.I),
        "amat_var_psi_var_T": re.compile(r"^\s*amat\(\s*var_psi\s*,\s*var_T\s*\)\s*=\s*(.*)$", re.I),
        "amat_n_var_psi_var_T": re.compile(r"^\s*amat_n\(\s*var_psi\s*,\s*var_T\s*\)\s*=\s*(.*)$", re.I),
        "amat_var_psi_var_Te": re.compile(r"^\s*amat\(\s*var_psi\s*,\s*var_Te\s*\)\s*=\s*(.*)$", re.I),
        "amat_n_var_psi_var_Te": re.compile(r"^\s*amat_n\(\s*var_psi\s*,\s*var_Te\s*\)\s*=\s*(.*)$", re.I),
        "amat_var_psi_var_rhoimp": re.compile(r"^\s*amat\(\s*var_psi\s*,\s*var_rhoimp\s*\)\s*=\s*(.*)$", re.I),
    }
    found = {}
    index = 0
    while index < len(text):
        line = text[index].split("!", 1)[0]
        key = next((name for name, pattern in patterns.items() if pattern.match(line)), None)
        if key is None:
            index += 1
            continue
        current = patterns[key].match(line).group(1)
        pieces = []
        while True:
            current = current.strip()
            continued = current.endswith("&")
            pieces.append(current[:-1] if continued else current)
            if not continued:
                break
            index += 1
            while index < len(text) and not text[index].split("!", 1)[0].strip():
                index += 1
            if index >= len(text):
                raise ValueError("Unterminated assignment for {}".format(key))
            current = text[index].split("!", 1)[0].lstrip()
            if current.startswith("&"):
                current = current[1:]
        found[key] = " ".join(piece.strip() for piece in pieces)
        index += 1
    missing = [name for name in MODEL600_EQ1_ASSIGNMENTS if name not in found]
    if missing:
        raise ValueError("Assignments not found: {}".format(", ".join(missing)))
    return found


MODEL600_SIMPLE_ASSIGNMENTS = (
    "rhs_var_zj",
    "amat_var_zj_var_zj",
    "amat_var_zj_var_psi",
    "rhs_var_w",
    "amat_var_w_var_w",
    "amat_var_w_var_u",
)


def extract_model600_simple_equations(path):
    """Extract the current and vorticity definition assignments."""

    text = Path(path).read_text(encoding="utf-8").splitlines()
    patterns = {
        "rhs_var_zj": re.compile(r"^\s*rhs_ij\(\s*var_zj\s*\)\s*=\s*(.*)$", re.I),
        "amat_var_zj_var_zj": re.compile(r"^\s*amat\(\s*var_zj\s*,\s*var_zj\s*\)\s*=\s*(.*)$", re.I),
        "amat_var_zj_var_psi": re.compile(r"^\s*amat\(\s*var_zj\s*,\s*var_psi\s*\)\s*=\s*(.*)$", re.I),
        "rhs_var_w": re.compile(r"^\s*rhs_ij\(\s*var_w\s*\)\s*=\s*(.*)$", re.I),
        "amat_var_w_var_w": re.compile(r"^\s*amat\(\s*var_w\s*,\s*var_w\s*\)\s*=\s*(.*)$", re.I),
        "amat_var_w_var_u": re.compile(r"^\s*amat\(\s*var_w\s*,\s*var_u\s*\)\s*=\s*(.*)$", re.I),
    }
    found = {}
    for index, raw in enumerate(text):
        line = raw.split("!", 1)[0]
        for key, pattern in patterns.items():
            match = pattern.match(line)
            if match:
                found[key] = match.group(1).strip()
                break
    missing = [name for name in MODEL600_SIMPLE_ASSIGNMENTS if name not in found]
    if missing:
        raise ValueError("Assignments not found: {}".format(", ".join(missing)))
    return found


def _generated_model600_simple_equations():
    result = {}
    zj = current_constraint_equation_zj().linearize(fields=FIELDS)
    w = vorticity_constraint_equation_w().linearize(fields=FIELDS)
    result["rhs_var_zj"] = zj.rhs
    result["amat_var_zj_var_zj"] = zj.amat[j]
    result["amat_var_zj_var_psi"] = zj.amat[psi]
    result["rhs_var_w"] = w.rhs
    result["amat_var_w_var_w"] = w.amat[omega]
    result["amat_var_w_var_u"] = w.amat[u]
    return {name: fortran(value) for name, value in result.items()}


def compare_model600_simple_equations(path):
    """Compare all RHS/AMAT blocks for ``var_zj`` and ``var_w``."""

    source = {
        name: _normalize_model600_text(value)
        for name, value in extract_model600_simple_equations(path).items()
    }
    generated = {
        name: _normalize_model600_text(value)
        for name, value in _generated_model600_simple_equations().items()
    }
    compare_assignment_maps(source, generated)


def generated_model600_equation1():
    return _generated_model600_equation1(with_TiTe=False)


def _generated_model600_equation1(*, with_TiTe=False):
    equation = induction_equation_1(with_TiTe=with_TiTe)
    linearized = equation.linearize(
        fields=FIELDS,
        timestep=coefficient("tstep"),
        theta=coefficient("theta"),
        zeta=coefficient("zeta"),
    )
    result = {}
    result["rhs_var_psi"] = split_toroidal_channels(linearized.rhs).p
    thermal_name, thermal_field = ("Te", Te) if with_TiTe else ("T", T)
    blocks = (
        ("amat_var_psi_var_psi", psi, "p"),
        ("amat_var_psi_var_u", u, "p"),
        ("amat_n_var_psi_var_u", u, "n"),
        ("amat_var_psi_var_zj", j, "p"),
        ("amat_var_psi_var_rho", rho, "p"),
        ("amat_n_var_psi_var_rho", rho, "n"),
        ("amat_var_psi_var_{}".format(thermal_name), thermal_field, "p"),
        ("amat_n_var_psi_var_{}".format(thermal_name), thermal_field, "n"),
        ("amat_var_psi_var_rhoimp", FIELDS[-1], "p"),
    )
    for name, field, channel in blocks:
        result[name] = getattr(split_toroidal_channels(linearized.amat[field]), channel)
    return {
        name: fortran(value, previous_names={"psi": "delta_g(mp,var_psi,ms,mt)"})
        for name, value in result.items()
    }


def compare_model600_equation1(path):
    source_all = {
        name: _normalize_model600_text(value)
        for name, value in extract_model600_equation1(path).items()
    }
    common = {
        name for name in MODEL600_EQ1_ASSIGNMENTS
        if "_var_T" not in name and "_var_Te" not in name
    }
    differences = []
    for with_TiTe, thermal in ((False, "T"), (True, "Te")):
        active = common | {
            name for name in MODEL600_EQ1_ASSIGNMENTS if "_var_{}".format(thermal) in name
        }
        source = {name: source_all[name] for name in active}
        generated = {
            name: _normalize_model600_text(value)
            for name, value in _generated_model600_equation1(with_TiTe=with_TiTe).items()
            if name in active
        }
        try:
            compare_assignment_maps(source, generated)
        except FortranComparisonError as error:
            differences.append(
                "with_TiTe={}:\n{}".format(with_TiTe, error)
            )
    if differences:
        raise FortranComparisonError("Model-600 equation-1 differences:\n" + "\n".join(differences))


def compare_model600_rhs_psi(path):
    """Compare only the single physical ``rhs(var_psi)`` expression."""

    source = extract_model600_equation1(path)
    source = {"rhs_var_psi": _normalize_model600_text(source["rhs_var_psi"])}
    generated = _generated_model600_equation1(with_TiTe=False)
    generated = {"rhs_var_psi": _normalize_model600_text(generated["rhs_var_psi"])}
    compare_assignment_maps(source, generated)


def compare_model600_amat_psi_psi(path):
    """Compare the common ``amat(var_psi,var_psi)`` block."""

    source = extract_model600_equation1(path)
    source = {
        "amat_var_psi_var_psi": _normalize_model600_text(
            source["amat_var_psi_var_psi"]
        )
    }
    generated = _generated_model600_equation1(with_TiTe=False)
    generated = {
        "amat_var_psi_var_psi": _normalize_model600_text(
            generated["amat_var_psi_var_psi"]
        )
    }
    compare_assignment_maps(source, generated)


def compare_model600_amat_psi_remaining(path):
    """Compare all remaining equation-1 AMAT blocks for both thermal branches."""

    source_all = extract_model600_equation1(path)
    differences = []
    for with_TiTe, thermal in ((False, "T"), (True, "Te")):
        names = (
            "amat_var_psi_var_u", "amat_n_var_psi_var_u",
            "amat_var_psi_var_zj", "amat_var_psi_var_rho",
            "amat_n_var_psi_var_rho",
            "amat_var_psi_var_{}".format(thermal),
            "amat_n_var_psi_var_{}".format(thermal),
            "amat_var_psi_var_rhoimp",
        )
        source = {
            name: _normalize_model600_text(
                source_all[name],
                single_temperature=not with_TiTe,
                two_temperature=with_TiTe,
            ) for name in names
        }
        generated_all = _generated_model600_equation1(with_TiTe=with_TiTe)
        generated = {
            name: _normalize_model600_text(
                generated_all[name],
                single_temperature=not with_TiTe,
                two_temperature=with_TiTe,
            ) for name in names
        }
        try:
            compare_assignment_maps(source, generated)
        except Exception as error:
            differences.append("with_TiTe={}:\n{}".format(with_TiTe, error))
    if differences:
        raise FortranComparisonError("Model-600 remaining AMAT differences:\n" + "\n".join(differences))


def compare_model600_amat_psi_single_temperature(path):
    """Compare only the single-temperature ``amat(var_psi,:)`` row.

    This intentionally excludes the two-temperature branch so the diamagnetic
    density/temperature tangent can be inspected in isolation.
    """

    source_all = extract_model600_equation1(path)
    names = (
        "amat_var_psi_var_psi",
        "amat_var_psi_var_u",
        "amat_n_var_psi_var_u",
        "amat_var_psi_var_zj",
        "amat_var_psi_var_rho",
        "amat_n_var_psi_var_rho",
        "amat_var_psi_var_T",
        "amat_n_var_psi_var_T",
        "amat_var_psi_var_rhoimp",
    )
    source = {
        name: _normalize_model600_text(source_all[name], single_temperature=True)
        for name in names
    }
    generated_all = _generated_model600_equation1(with_TiTe=False)
    generated = {
        name: _normalize_model600_text(generated_all[name], single_temperature=True)
        for name in names
    }
    compare_assignment_maps(source, generated)


def _normalize_model600_text(
    expression, *, single_temperature=False, two_temperature=False
):
    # Convert the physical-coordinate cross product used by the source to its
    # equivalent element-coordinate bracket before expanding aliases.
    expression = re.sub(
        r"ps0_x\s*\*\s*Pe0_y\s*-\s*ps0_y\s*\*\s*Pe0_x",
        "(ps0_s*Pe0_t-ps0_t*Pe0_s)/xjac",
        expression,
        flags=re.I,
    )
    expression = re.sub(
        r"factor\(\s*var_psi\s*,\s*([1-6])\s*\)",
        "1",
        expression,
        flags=re.I,
    )
    expression = re.sub(
        r"factor\(\s*var_(?:zj|w)\s*,\s*([1-6])\s*\)",
        "1",
        expression,
        flags=re.I,
    )
    expression = re.sub(
        r"delta_g\(\s*mp\s*,\s*var_psi\s*,\s*ms\s*,\s*mt\s*\)",
        "delta_psi",
        expression,
        flags=re.I,
    )
    expression = re.sub(
        r"current_source\(\s*ms\s*,\s*mt\s*\)",
        "current_source",
        expression,
        flags=re.I,
    )
    # The symbolic DSL prints the supplied density correction as a function
    # call, while the element routine stores its value and derivative in the
    # work variables ``r0_corr`` and ``dr0_corr_dn``.
    expression = re.sub(
        r"corr_neg_dens\(\s*r0\s*\)", "r0_corr", expression, flags=re.I
    )
    expression = re.sub(r"\bcorr_neg_dens\b", "r0_corr", expression, flags=re.I)
    expression = re.sub(
        r"dr0_corr_dn\(\s*r0\s*\)", "dr0_corr_dn", expression, flags=re.I
    )
    # The element routine's displayed tangent is written for the default
    # correction branch where d(r0_corr)/d(rho)=1; compare that convention
    # with the explicit DSL derivative call.
    expression = re.sub(r"\bdr0_corr_dn\b", "1", expression)
    pressure_temperature = "Te0" if two_temperature else "T0"
    aliases = {
        "BB2": "((F0**2+ps0_x**2+ps0_y**2)/BigR**2)",
        "Pe0": "(r0*{})".format(pressure_temperature),
        "Pe0_s": "(r0_s*{}+r0*{}_s)".format(pressure_temperature, pressure_temperature),
        "Pe0_t": "(r0_t*{}+r0*{}_t)".format(pressure_temperature, pressure_temperature),
        "Pe0_p": "(r0_p*{}+r0*{}_p)".format(pressure_temperature, pressure_temperature),
        "Pe0_x": "(r0_x*{}+r0*{}_x)".format(pressure_temperature, pressure_temperature),
        "Pe0_y": "(r0_y*{}+r0*{}_y)".format(pressure_temperature, pressure_temperature),
    }
    for name, replacement in aliases.items():
        expression = re.sub(r"\b{}\b".format(name), replacement, expression)
    if single_temperature:
        expression = re.sub(r"\bTe0([_a-z]*)\b", r"T0\1", expression)
    return expression
