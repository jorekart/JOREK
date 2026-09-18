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
    ion_temperature = "Ti0" if two_temperature else "T0"
    aliases = {
        "BB2": "((F0**2+ps0_x**2+ps0_y**2)/BigR**2)",
        "Pe0": "(r0*{})".format(pressure_temperature),
        "Pe0_s": "(r0_s*{}+r0*{}_s)".format(pressure_temperature, pressure_temperature),
        "Pe0_t": "(r0_t*{}+r0*{}_t)".format(pressure_temperature, pressure_temperature),
        "Pe0_p": "(r0_p*{}+r0*{}_p)".format(pressure_temperature, pressure_temperature),
        "Pe0_x": "(r0_x*{}+r0*{}_x)".format(pressure_temperature, pressure_temperature),
        "Pe0_y": "(r0_y*{}+r0*{}_y)".format(pressure_temperature, pressure_temperature),
        "Pi0": "(r0*{})".format(ion_temperature),
        "Pi0_s": "(r0_s*{}+r0*{}_s)".format(ion_temperature, ion_temperature),
        "Pi0_t": "(r0_t*{}+r0*{}_t)".format(ion_temperature, ion_temperature),
        "Pi0_y": "(r0_y*{}+r0*{}_y)".format(ion_temperature, ion_temperature),
    }
    for name, replacement in aliases.items():
        expression = re.sub(r"\b{}\b".format(name), replacement, expression)
    if single_temperature:
        expression = re.sub(r"\bTe0([_a-z]*)\b", r"T0\1", expression)
    return expression
