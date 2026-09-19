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
    momentum_equation_2,
    j,
    omega,
    rho,
    rhon,
    rhoimp,
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

MODEL600_U_BASE_ASSIGNMENTS = (
    "amat_var_u_var_zj",
    "amat_n_var_u_var_zj",
)

MODEL600_U_PSI_W_ASSIGNMENTS = (
    "amat_var_u_var_psi",
    "amat_var_u_var_w",
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


def _extract_model600_u_base(path):
    text = Path(path).read_text(encoding="utf-8").splitlines()
    patterns = {
        "amat_var_u_var_zj": re.compile(
            r"^\s*amat\(\s*var_u\s*,\s*var_zj\s*\)\s*=\s*(.*)$", re.I
        ),
        "amat_n_var_u_var_zj": re.compile(
            r"^\s*amat_n\(\s*var_u\s*,\s*var_zj\s*\)\s*=\s*(.*)$", re.I
        ),
    }
    found = {}
    index = 0
    while index < len(text):
        line = text[index].split("!", 1)[0]
        match_item = next(
            (
                (name, match)
                for name, pattern in patterns.items()
                for match in (pattern.match(line),)
                if match is not None
            ),
            None,
        )
        if match_item is None:
            index += 1
            continue
        name, match = match_item
        if name in found:
            index += 1
            continue
        pieces = []
        current = match.group(1)
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
                break
            current = text[index].split("!", 1)[0].strip()
            if current.startswith("&"):
                current = current[1:]
        found[name] = " ".join(pieces)
        index += 1
    return found


def compare_model600_u_base_zj(path):
    """Compare the unconditional ``var_u`` current-coupling blocks."""

    source = {
        name: _normalize_model600_text(value)
        for name, value in _extract_model600_u_base(path).items()
    }
    equation = momentum_equation_2()
    linearized = equation.linearize(
        fields=FIELDS,
        timestep=coefficient("tstep"),
        theta=coefficient("theta"),
        zeta=coefficient("zeta"),
    )
    generated = {
        "amat_var_u_var_zj": fortran(linearized.amat[j]),
        "amat_n_var_u_var_zj": fortran(linearized.amat[j]),
    }
    # The p and n channels are separated after linearization.
    generated["amat_var_u_var_zj"] = fortran(
        split_toroidal_channels(linearized.amat[j]).p
    )
    generated["amat_n_var_u_var_zj"] = fortran(
        split_toroidal_channels(linearized.amat[j]).n
    )
    generated = {
        name: _normalize_model600_text(value) for name, value in generated.items()
    }
    compare_assignment_maps(source, generated)


def _extract_model600_u_assignments(path):
    text = Path(path).read_text(encoding="utf-8").splitlines()
    patterns = {
        "amat_var_u_var_psi": re.compile(
            r"^\s*amat\(\s*var_u\s*,\s*var_psi\s*\)\s*=\s*(.*)$", re.I
        ),
        "amat_var_u_var_w": re.compile(
            r"^\s*amat\(\s*var_u\s*,\s*var_w\s*\)\s*=\s*(.*)$", re.I
        ),
    }
    found = {}
    index = 0
    while index < len(text):
        line = text[index].split("!", 1)[0]
        match_item = next(
            (
                (name, match)
                for name, pattern in patterns.items()
                for match in (pattern.match(line),)
                if match is not None
            ),
            None,
        )
        if match_item is None:
            index += 1
            continue
        name, match = match_item
        if name in found:
            index += 1
            continue
        pieces = []
        current = match.group(1)
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
                break
            current = text[index].split("!", 1)[0].strip()
            if current.startswith("&"):
                current = current[1:]
        found[name] = " ".join(pieces)
        index += 1
    return found


def compare_model600_u_psi_w(path):
    """Compare the base-generated psi and vorticity AMAT columns."""

    source = {
        name: _normalize_model600_text(value)
        for name, value in _extract_model600_u_assignments(path).items()
    }
    equation = momentum_equation_2(
        include_diamagnetic=True,
        include_conservative=True,
        include_tgnum=True,
    )
    linearized = equation.linearize(
        fields=(psi, omega),
        timestep=coefficient("tstep"),
        theta=coefficient("theta"),
        zeta=coefficient("zeta"),
    )
    generated = {
        "amat_var_u_var_psi": fortran(
            split_toroidal_channels(linearized.amat[psi]).p
        ),
        "amat_var_u_var_w": fortran(
            split_toroidal_channels(linearized.amat[omega]).p
        ),
    }
    generated = {
        name: _normalize_model600_text(value) for name, value in generated.items()
    }
    compare_assignment_maps(source, generated)


def check_model600_u_neo_generation():
    """Ensure the optional NEO momentum residual can be constructed."""

    momentum_equation_2(include_neo=True)


def check_model600_u_neutral_generation():
    """Ensure the optional neutral momentum residual can be constructed."""

    momentum_equation_2(include_neutrals=True)


def check_model600_u_neutral_amat_generation():
    """Build neutral-only AMAT columns for both thermal configurations."""

    for with_TiTe, thermal in ((False, T), (True, Te)):
        equation = momentum_equation_2(
            with_TiTe=with_TiTe,
            include_neutrals=True,
            neutral_only=True,
        )
        linearized = equation.linearize(
            fields=(rho, thermal, rhon, rhoimp),
            timestep=coefficient("tstep"),
            theta=coefficient("theta"),
            zeta=coefficient("zeta"),
        )
        for field in (rho, thermal, rhon, rhoimp):
            fortran(linearized.amat[field])


def check_model600_u_neutral_source_generation():
    """Verify all external particle/source contributions reach the residual."""

    equation = momentum_equation_2(include_neutrals=True, neutral_only=True)
    rendered = fortran(equation.B)
    required = (
        "particle_source",
        "source_pellet",
        "source_bg_drift",
        "source_imp_drift",
    )
    missing = [name for name in required if name not in rendered]
    if missing:
        raise AssertionError(
            "Neutral source terms missing from var_u residual: {}".format(
                ", ".join(missing)
            )
        )


def check_model600_u_delta_n_convection_generation():
    """Verify neutral residual/tangents retain the convection switch."""

    for with_TiTe, thermal in ((False, T), (True, Te)):
        equation = momentum_equation_2(
            with_TiTe=with_TiTe,
            include_neutrals=True,
            neutral_only=True,
        )
        residual = fortran(equation.B)
        if "delta_n_convection" not in residual:
            raise AssertionError("delta_n_convection missing from neutral residual")
        linearized = equation.linearize(
            fields=(rho, thermal, rhon, rhoimp),
            timestep=coefficient("tstep"),
            theta=coefficient("theta"),
            zeta=coefficient("zeta"),
        )
        if not any(
            "delta_n_convection" in fortran(linearized.amat[field])
            for field in (rho, thermal, rhon, rhoimp)
        ):
            raise AssertionError("delta_n_convection missing from neutral AMAT")


def check_model600_u_impurity_pressure_generation():
    """Ensure impurity pressure reaches both residual temperature branches."""

    for with_TiTe in (False, True):
        equation = momentum_equation_2(
            with_TiTe=with_TiTe,
            include_impurities=True,
        )
        rendered = fortran(equation.B)
        if "rhoimp" not in rendered:
            raise AssertionError("rhoimp missing from impurity-pressure residual")


def check_model600_u_full_amat_generation():
    """Generate every non-NEO var_u AMAT column and FFT channel."""

    for with_TiTe in (False, True):
        equation = momentum_equation_2(
            with_TiTe=with_TiTe,
            include_diamagnetic=True,
            include_conservative=True,
            include_tgnum=True,
            include_neutrals=True,
            include_impurities=True,
        )
        linearized = equation.linearize(
            fields=FIELDS,
            timestep=coefficient("tstep"),
            theta=coefficient("theta"),
            zeta=coefficient("zeta"),
        )
        for field in FIELDS:
            # The current channel splitter intentionally rejects second
            # toroidal derivatives; those belong to JOREK's ``amat_nn``
            # channel.  Still render the complete column so nn terms are
            # validated for printer support.
            fortran(linearized.amat[field])


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
    # The element routine uses short background-density names, whereas the
    # symbolic model uses the descriptive field names.  They denote the same
    # frozen quantities in the single-temperature momentum equation.
    expression = re.sub(r"\brn0(?=\b|_)", "rhon0", expression, flags=re.I)
    expression = re.sub(r"\brimp0(?=\b|_)", "rhoimp0", expression, flags=re.I)
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
    # ``delta_g(mp,var_X,ms,mt)`` is the previous Newton increment of field X.
    # Give every field the same short spelling the DSL prints.
    expression = re.sub(
        r"delta_g\s*\(\s*mp\s*,\s*var_([a-z0-9_]+)\s*,\s*ms\s*,\s*mt\s*\)",
        r"delta_\1",
        expression,
        flags=re.I,
    )
    expression = re.sub(
        r"current_source\(\s*ms\s*,\s*mt\s*\)",
        "current_source",
        expression,
        flags=re.I,
    )
    # NEO profile work arrays are scalar values at the element/toroidal
    # point, despite being called as indexed Fortran functions.
    for profile in ("amu_neo_prof", "aki_neo_prof"):
        expression = re.sub(
            r"\b{}\s*\(\s*ms\s*,\s*mt\s*\)".format(profile),
            profile,
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
    if single_temperature:
        # The one-temperature model evolves the total temperature.  Do not
        # identify either Ti0 or Te0 with T0 individually; only their sum
        # (and the corresponding directional/spatial derivatives) is T0.
        for suffix in ("", "_s", "_t", "_p", "_x", "_y", "_xx", "_yy", "_xy"):
            expression = re.sub(
                r"\bTi0{}\s*\+\s*Te0{}\b".format(suffix, suffix),
                "T0{}".format(suffix), expression, flags=re.I,
            )
            expression = re.sub(
                r"\bTe0{}\s*\+\s*Ti0{}\b".format(suffix, suffix),
                "T0{}".format(suffix), expression, flags=re.I,
            )
        expression = re.sub(
            r"\balpha_imp_T\b", "(alpha_imp*T0)", expression, flags=re.I,
        )
    # ``alpha_e_T`` is the DSL print name of the closure value alpha_e*Te0.
    # It is spelled the same way in both temperature models; the single
    # temperature substitution Te0 = T0/2 is applied at the end.
    expression = re.sub(
        r"\balpha_e_T\b", "(alpha_e*Te0)", expression, flags=re.I,
    )
    pressure_temperature = "Te0"
    ion_temperature = "Ti0"
    # ``construct_pressure`` always includes impurity pressure when the
    # impurity extension is active.  Reports compare that full definition,
    # rather than treating P0/Pi0/Pe0 as opaque element work variables.
    # ``construct_pressure`` is written once for both temperature models: it
    # always builds the species pressures from Ti0/Te0 and the per-species
    # impurity coefficients.  The one-temperature branch reaches it with
    # Ti0=Te0=T0/2, which the substitutions at the end of this function
    # apply.  Expanding the aliases in the two-species basis therefore
    # reproduces the element routine exactly in both branches.
    ion_alpha = "alpha_i"
    ion_density = "(r0+rhoimp0*alpha_i)"
    electron_density = "(r0+rhoimp0*alpha_e)"
    electron_density_tangent = "(r0+rhoimp0*alpha_e_bis)"

    pi0 = "({}*{})".format(ion_density, ion_temperature)
    pe0 = "({}*{})".format(electron_density, pressure_temperature)
    pi_derivatives = {
        "s": "((r0_s+rhoimp0_s*{a})*{t}+{d}*{t}_s)",
        "t": "((r0_t+rhoimp0_t*{a})*{t}+{d}*{t}_t)",
        "p": "((r0_p+rhoimp0_p*{a})*{t}+{d}*{t}_p)",
        "x": "((r0_x+rhoimp0_x*{a})*{t}+{d}*{t}_x)",
        "y": "((r0_y+rhoimp0_y*{a})*{t}+{d}*{t}_y)",
        "xx": "((r0_xx+rhoimp0_xx*{a})*{t}+2*(r0_x+rhoimp0_x*{a})*{t}_x+{d}*{t}_xx)",
        "yy": "((r0_yy+rhoimp0_yy*{a})*{t}+2*(r0_y+rhoimp0_y*{a})*{t}_y+{d}*{t}_yy)",
        "xy": "((r0_xy+rhoimp0_xy*{a})*{t}+(r0_x+rhoimp0_x*{a})*{t}_y+(r0_y+rhoimp0_y*{a})*{t}_x+{d}*{t}_xy)",
    }
    pi_derivatives = {
        key: value.format(a=ion_alpha, d=ion_density, t=ion_temperature)
        for key, value in pi_derivatives.items()
    }
    pe_derivatives = {
        "s": "((r0_s+rhoimp0_s*alpha_e)*{t}+{dt}*{t}_s)",
        "t": "((r0_t+rhoimp0_t*alpha_e)*{t}+{dt}*{t}_t)",
        "p": "((r0_p+rhoimp0_p*alpha_e)*{t}+{dt}*{t}_p)",
        "x": "((r0_x+rhoimp0_x*alpha_e)*{t}+{dt}*{t}_x)",
        "y": "((r0_y+rhoimp0_y*alpha_e)*{t}+{dt}*{t}_y)",
    }
    pe_derivatives = {
        key: value.format(d=electron_density, dt=electron_density_tangent,
                          t=pressure_temperature)
        for key, value in pe_derivatives.items()
    }
    p0 = "({}+{})".format(pi0, pe0)
    p0_derivatives = {
        key: "({}+{})".format(pi_derivatives[key], pe_derivatives[key])
        for key in ("s", "t", "p")
    }
    # The implicit heating floor of the ion and electron energy equations is
    # written inline with ``min`` and ``exp``.  Name the two pieces so that
    # the expression parser sees ordinary work values, and resolve the
    # derivative of the exponential through its own value.
    expression = re.sub(
        r"\bmin\s*\(\s*Ti0\s*,\s*Tie_min_neg\s*\)", "Ti0_floor",
        expression, flags=re.I,
    )
    expression = re.sub(
        r"\bexp\s*\(\s*\(\s*Ti0_floor\s*-\s*Tie_min_neg\s*\)\s*/\s*"
        r"\(\s*0\.5(?:d0)?\s*\*\s*Tie_min_neg\s*\)\s*\)",
        "Ti_floor_exp", expression, flags=re.I,
    )
    expression = re.sub(
        r"\bdTi_floor_exp\b", "(Ti_floor_exp/(0.5*Tie_min_neg))",
        expression, flags=re.I,
    )
    expression = re.sub(r"\bdTi_floor\b", "1", expression, flags=re.I)
    # Impurity negative-density correction, mirroring ``corr_neg_dens``.
    expression = re.sub(
        r"corr_neg_dens_imp\(\s*rhoimp0?\s*\)", "rhoimp0_corr",
        expression, flags=re.I,
    )
    expression = re.sub(
        r"\bcorr_neg_dens_imp\b", "rhoimp0_corr", expression, flags=re.I,
    )
    expression = re.sub(r"\bdrimp0_corr_dn\b", "1", expression, flags=re.I)
    expression = re.sub(
        r"\brimp0_corr\b", "rhoimp0_corr", expression, flags=re.I,
    )
    # ``sqrt`` is the only Fortran intrinsic appearing in the exported volume
    # terms.  Rewrite it as a power so the expression parser needs no function
    # table.
    expression = re.sub(
        r"\bsqrt\s*\(([^()]*)\)", r"((\1)**(1/2))", expression, flags=re.I,
    )
    # Fortran is case insensitive; the density and parallel-velocity blocks
    # spell the parallel-velocity trial function ``Vpar`` while the DSL prints
    # ``vpar``.  ``Vpar0`` is covered by its own alias below.
    expression = re.sub(
        r"\bVpar(_[a-z]+)?\b", r"vpar\1", expression, flags=re.I,
    )
    aliases = {
        "BB2": "((F0**2+ps0_x**2+ps0_y**2)/BigR**2)",
        "psi_grad2": "(ps0_x**2+ps0_y**2)",
        # Parallel-gradient work values of the density equation.  The density
        # ones are already covered by ``normalize_fortran_text``; the impurity
        # ones must be expanded after ``rimp0`` has been renamed.
        "Bgrad_rhoimp": "((F0*rhoimp0_p/BigR+rhoimp0_x*ps0_y-rhoimp0_y*ps0_x)/BigR)",
        "Bgrad_rhoimp_psi": "((rhoimp0_x*psi_y-rhoimp0_y*psi_x)/BigR)",
        "Bgrad_rhoimp_rhoimp": "((rhoimp_x*ps0_y-rhoimp_y*ps0_x)/BigR)",
        "Bgrad_rhoimp_rhoimp_n": "(F0*rhoimp_p/BigR**2)",
        # Parallel-velocity work values.
        # Ion-temperature parallel-gradient work values.
        "Bgrad_Ti": "((F0*Ti0_p/BigR+Ti0_x*ps0_y-Ti0_y*ps0_x)/BigR)",
        "Bgrad_Ti_psi": "((Ti0_x*psi_y-Ti0_y*psi_x)/BigR)",
        "Bgrad_Ti_Ti": "((Ti_x*ps0_y-Ti_y*ps0_x)/BigR)",
        "Bgrad_Ti_Ti_n": "(F0*Ti_p/BigR**2)",
        "Bgrad_vpar": "((F0*vpar0_p/BigR+vpar0_x*ps0_y-vpar0_y*ps0_x)/BigR)",
        "Bgrad_vpar_psi": "((vpar0_x*psi_y-vpar0_y*psi_x)/BigR)",
        "Bgrad_vpar_vpar": "((vpar_x*ps0_y-vpar_y*ps0_x)/BigR)",
        # Prescribed rotation profile: a flux function whose flux derivative
        # is stored as ``dV_dpsi_source``.
        "Vt0_x": "(dV_dpsi_source*ps0_x)",
        "Vt0_y": "(dV_dpsi_source*ps0_y)",
        "Vt_x_psi": "(dV_dpsi_source*psi_x)",
        "Vt_y_psi": "(dV_dpsi_source*psi_y)",
        # Fortran is case insensitive: amat_n(var_vpar,var_Ti) spells the
        # background density ``R0``.
        "R0": "r0",
        "Btheta2": "((ps0_x**2+ps0_y**2)/BigR**2)",
        "Btheta2_psi": "(2*(psi_x*ps0_x+psi_y*ps0_y)/BigR**2)",
        "Vpar0": "vpar0",

        "Pe0": pe0,
        "Pe0_s": pe_derivatives["s"],
        "Pe0_t": pe_derivatives["t"],
        "Pe0_p": pe_derivatives["p"],
        "Pe0_x": pe_derivatives["x"],
        "Pe0_y": pe_derivatives["y"],
        "Pi0": pi0,
        "Pi0_s": pi_derivatives["s"],
        "Pi0_t": pi_derivatives["t"],
        "Pi0_y": pi_derivatives["y"],
        # Cartesian derivatives used by the diamagnetic viscosity block.
        # Keep these as explicit product rules so source and DSL expressions
        # are expanded to the same monomials by the report matcher.
        "Pi0_x": pi_derivatives["x"],
        "Pi0_xx": pi_derivatives["xx"],
        "Pi0_yy": pi_derivatives["yy"],
        "Pi0_xy": pi_derivatives["xy"],
        "P0": p0,
        "P0_s": p0_derivatives["s"],
        "P0_t": p0_derivatives["t"],
        "P0_p": p0_derivatives["p"],
    }
    for name, replacement in aliases.items():
        expression = re.sub(
            r"\b{}\b".format(name), replacement, expression, flags=re.I
        )
    if single_temperature:
        # The element routine sets Ti0 = Te0 = T0/2 in the one-temperature
        # branch (the evolved T0 is the total temperature), so a species
        # temperature reaching this point is half of the evolved one.
        expression = re.sub(r"\bTe0([_a-z]*)\b", r"(T0\1/2)", expression)
        expression = re.sub(r"\bTi0([_a-z]*)\b", r"(T0\1/2)", expression)
        # The one-temperature impurity closure stores the means of the
        # per-species coefficients:  alpha_imp = (alpha_i+alpha_e)/2,
        # alpha_imp_bis = (alpha_i+alpha_e_bis)/2 (alpha_i does not depend on
        # temperature, so it is its own "bis"), and alpha_imp_tri =
        # alpha_e_tri/4, the extra quarter coming from dTe0/dT0 = 1/2 applied
        # twice.  Express everything in the two-species basis so that
        # assignments shared between the branches and assignments written for
        # one branch only use the same symbols.
        expression = re.sub(
            r"\balpha_imp_tri\b", "(alpha_e_tri/4)", expression, flags=re.I,
        )
        expression = re.sub(
            r"\balpha_imp_bis\b", "((alpha_i+alpha_e_bis)/2)",
            expression, flags=re.I,
        )
        expression = re.sub(
            r"\balpha_imp\b", "((alpha_i+alpha_e)/2)", expression, flags=re.I,
        )
        # W_dia_Ti is assembled from Pi0_*_Ti, which the element routine
        # builds from the *trial* temperature basis function without the
        # dTi0/dT0 = 1/2 chain factor.  It is therefore twice the derivative
        # of W_dia with respect to the evolved one-temperature field.
        expression = re.sub(
            r"\bW_dia_Ti\b", "(2*W_dia_T)", expression, flags=re.I,
        )
    return expression
