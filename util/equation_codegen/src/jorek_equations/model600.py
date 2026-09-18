"""Model-600 weak equations, extending the model-199 DSL definitions."""

import sympy as sp

from .equations import ConstraintEquation, EvolutionEquation
from .external import external_function
from .model199 import FIELDS as MODEL199_FIELDS
from .model199 import j, omega, psi, rho, T, u, xjac
from .operators import R, dR, dZ, dphi, ds, dt, element_bracket, dot, grad
from .symbols import coefficient, field, freeze, test_function


# Model-600 extension fields.  The model199 declarations are reused so their
# current/trial and Fortran naming conventions remain identical.
vpar = field("vpar", fortran_current="vpar0", fortran_trial="vpar")
Ti = field("Ti", fortran_current="Ti0", fortran_trial="Ti")
Te = field("Te", fortran_current="Te0", fortran_trial="Te")
rhon = field("rhon", fortran_current="rn0", fortran_trial="rhon")
rhoimp = field("rhoimp", fortran_current="rimp0", fortran_trial="rhoimp")

FIELDS = MODEL199_FIELDS + (vpar, Ti, Te, rhon, rhoimp)


Jb = coefficient("Jb")
aux_jre_ind = coefficient("aux_jre_ind")
tauIC = coefficient("tauIC")
r0_corr = coefficient("r0_corr")
F0 = coefficient("F0")
factor_psi = tuple(coefficient("factor_psi_{}".format(index)) for index in range(1, 7))

eta = external_function(
    "eta600",
    arguments=("T", "rho", "rhoimp"),
    derivatives={
        "T": "deta_dT",
        "rho": "deta_dr0",
        "rhoimp": "deta_drimp0",
    },
    policy="piecewise_active",
    fortran_name="eta_T",
)
eta_num_T = external_function(
    "eta_num600",
    arguments=("T",),
    derivatives={"T": "deta_num_dT"},
    policy="piecewise_active",
    fortran_name="eta_num_T",
)

visco = external_function(
    "visco600",
    arguments=("T",),
    derivatives={"T": "dvisco_dT"},
    policy="piecewise_active",
    fortran_name="visco_T",
)
visco_num = external_function(
    "visco_num600",
    arguments=("T",),
    derivatives={"T": "dvisco_num_dT"},
    policy="piecewise_active",
    fortran_name="visco_num_T",
)
visco_fact_old = coefficient("visco_fact_old")
visco_fact_new = coefficient("visco_fact_new")
dvisco_dT = coefficient("dvisco_dT")
fact_conservative_u = coefficient("fact_conservative_u")
delta_n_convection = coefficient("delta_n_convection")
tgnum_u = coefficient("tgnum_u")
amu_neo_prof = coefficient("amu_neo_prof")
aki_neo_prof = coefficient("aki_neo_prof")
Btheta2 = coefficient("Btheta2")
epsil = coefficient("epsil")
alpha_e = coefficient("alpha_e")
alpha_i = coefficient("alpha_i")
alpha_imp = coefficient("alpha_imp")

dvisco_state = external_function(
    "dvisco_state600", arguments=("T",),
    derivatives={"T": "d2visco_dT2"},
    policy="piecewise_active", fortran_name="dvisco_dT",
)
Sion_T = coefficient("Sion_T")
Srec_T = coefficient("Srec_T")
particle_source = coefficient("particle_source")
source_pellet = coefficient("source_pellet")
source_bg_drift = coefficient("source_bg_drift")
source_imp_drift = coefficient("source_imp_drift")
aux_P_par_re = coefficient("aux_P_par_re")
aux_P_perp_re = coefficient("aux_P_perp_re")
aux_divPIR_perp = coefficient("aux_divPIR_perp")
aux_divPIZ_perp = coefficient("aux_divPIZ_perp")

# Neutral rates and electron fraction are state-dependent work variables in
# model 600.  Declaring their interfaces lets the Newton generator produce the
# ``dS*_dT`` and ``dalpha_e_dT`` AMAT terms instead of treating the rates as
# frozen coefficients.
Sion_rate = external_function(
    "Sion_rate600", arguments=("T",), derivatives={"T": "dSion_dT"},
    policy="piecewise_active", fortran_name="Sion_T",
)
Srec_rate = external_function(
    "Srec_rate600", arguments=("T",), derivatives={"T": "dSrec_dT"},
    policy="piecewise_active", fortran_name="Srec_T",
)
alpha_e_state = external_function(
    "alpha_e600", arguments=("T",), derivatives={"T": "dalpha_e_dT"},
    policy="piecewise_active", fortran_name="alpha_e",
)

# Charge-state coefficients used in the impurity pressure have their own
# first/second temperature derivatives in the element routine (``bis`` and
# ``tri``).  Register the derivative calls as external functions too, so a
# second directional derivative remains well defined during AMAT generation.
alpha_e_bis_state = external_function(
    "alpha_e_bis", arguments=("T",), derivatives={"T": "alpha_e_tri"},
    policy="piecewise_active", fortran_name="alpha_e_bis",
)
alpha_e_pressure = external_function(
    "alpha_e_pressure600", arguments=("T",), derivatives={"T": "alpha_e_bis"},
    policy="piecewise_active", fortran_name="alpha_e",
)
alpha_imp_bis_state = external_function(
    "alpha_imp_bis", arguments=("T",), derivatives={"T": "alpha_imp_tri"},
    policy="piecewise_active", fortran_name="alpha_imp_bis",
)
alpha_imp_pressure = external_function(
    "alpha_imp_pressure600", arguments=("T",),
    derivatives={"T": "alpha_imp_bis"},
    policy="piecewise_active", fortran_name="alpha_imp",
)

# Density correction used by the model-600 resistive/pressure terms.  It is
# active in the Newton tangent; the Fortran routine supplies both its value
# and ``dr0_corr_dn``.  Keeping it as an external function (rather than a
# frozen coefficient) lets the DSL generate the corresponding density terms.
corr_neg_dens = external_function(
    "corr_neg_dens",
    arguments=("rho",),
    derivatives={"rho": "dr0_corr_dn"},
    policy="piecewise_active",
    fortran_name="corr_neg_dens",
)

# ``W_dia`` is assembled in the element routine and its state derivatives are
# supplied as work variables.  Declare that interface explicitly instead of
# expanding the long pressure-gradient definition inside every tangent.  This
# keeps the generated weak form identical to the Fortran decomposition while
# retaining the required derivatives for AMAT generation.
W_dia_single = external_function(
    "W_dia_single600",
    arguments=("rho", "T"),
    derivatives={"rho": "W_dia_rho", "T": "W_dia_T"},
    policy="piecewise_active",
    fortran_name="W_dia",
)
W_dia_two = external_function(
    "W_dia_two600",
    arguments=("rho", "Ti", "Te"),
    derivatives={
        "rho": "W_dia_rho",
        "Ti": "W_dia_Ti",
        "Te": "W_dia_Te",
    },
    policy="piecewise_active",
    fortran_name="W_dia",
)


def induction_equation_1(
    *,
    include_pressure_coupling: bool = True,
    include_runaway_coupling: bool = True,
    with_TiTe: bool = False,
):
    """Return model-600 ``rhs(var_psi)`` in integrated weak form.

    The two keyword switches correspond to the optional pressure and kinetic
    runaway-electron terms in ``mod_elt_matrix_fft.f90``.  The ``factor``
    switches are represented by independent coefficients so term-by-term
    generation can preserve the source routine's decomposition.
    """

    v = test_function("v")
    thermal = Te if with_TiTe else T
    # In model600 ``r0_corr`` is the corrected density, not an independent
    # frozen coefficient.  The equation below uses the correction only in
    # the pressure coupling; resistivity still depends on the physical rho.
    corrected_rho = corr_neg_dens(rho)
    psi_term = (
        v * eta(thermal, rho, rhoimp)
        * (j - coefficient("current_source") - Jb) / R * xjac
    )
    advection = v * element_bracket(psi, u) - v * F0 / R * dphi(u) * xjac
    resistive = eta_num_T(thermal) * dot(grad(v), grad(j)) * xjac
    base_B = psi_term + advection + resistive
    B = base_B
    single_temperature_tangent_B = base_B
    if include_pressure_coupling:
        Pe = rho * thermal
        diamagnetic_core = (
            -v * tauIC / (corrected_rho * _parallel_norm(psi))
            * F0**2 / R**2 * element_bracket(psi, Pe)
            + v * tauIC / (corrected_rho * _parallel_norm(psi))
            * F0**3 / R**3 * dphi(Pe) * xjac
        )
        # The physical RHS uses 2*tauIC.  In the legacy one-temperature
        # Jacobian, however, the T tangent is coded with tauIC; retain that
        # convention explicitly without changing the RHS or other tangents.
        B += 2 * diamagnetic_core
        single_temperature_tangent_B += diamagnetic_core
    if include_runaway_coupling:
        runaway = -v * eta(thermal, rho, rhoimp) * aux_jre_ind / R * xjac
        B += runaway
        single_temperature_tangent_B += runaway
    A = v * psi / R * xjac
    overrides = {T: single_temperature_tangent_B} if not with_TiTe else {}
    return EvolutionEquation(
        "model600_induction", v, A, B,
        amat_variation_overrides=overrides,
    )


def current_constraint_equation_zj():
    """Model-600 current-definition equation for ``var_zj``."""

    v = test_function("v")
    C = (dot(grad(v), grad(psi)) + v * j) / R * xjac
    return ConstraintEquation("model600_current_constraint", v, C, kind="static")


def vorticity_constraint_equation_w():
    """Model-600 vorticity-definition equation for ``var_w``."""

    v = test_function("v")
    C = (dot(grad(v), grad(u)) + v * omega) * R * xjac
    return ConstraintEquation("model600_vorticity_constraint", v, C, kind="static")


def _laplacian(value):
    return dR(dR(value)) + dZ(dZ(value)) + dR(value) / R


def _rho_hat(value=rho):
    return R**2 * value


def _velocity_norm(value=u):
    return R**2 * (dR(value)**2 + dZ(value)**2)


def _poloidal_cross(left, right):
    return dR(left) * dZ(right) - dZ(left) * dR(right)


def _diamagnetic_pressure(*, with_TiTe=False, include_impurities=False):
    """Ion pressure used by the diamagnetic momentum terms."""

    ion_temperature = Ti if with_TiTe else T
    pressure = rho * ion_temperature
    if include_impurities:
        impurity_alpha = alpha_i if with_TiTe else alpha_imp
        pressure += rhoimp * impurity_alpha * ion_temperature
    return pressure


def _diamagnetic_viscosity(*, with_TiTe=False, include_impurities=False):
    """Return the externally supplied model-600 ``W_dia`` coefficient."""

    if with_TiTe:
        return W_dia_two(rho, Ti, Te)
    return W_dia_single(rho, T)


def momentum_equation_2(
    *,
    with_TiTe=False,
    include_extensions=False,
    include_diamagnetic=False,
    include_conservative=False,
    include_tgnum=False,
    include_neo=False,
    neo_only=False,
    include_neutrals=False,
    include_impurities=False,
    neutral_only=False,
    include_auxiliary=True,
):
    """Return the unconditional model-600 perpendicular-momentum equation.

    ``include_extensions`` is deliberately separate: NEO, diamagnetic
    viscosity, conservative-form, neutral, impurity, and parallel-velocity
    terms are added in subsequent model slices after their source assignments
    have independent comparisons.
    """

    v = test_function("v")
    thermal = Te if with_TiTe else T
    alpha_e_value = alpha_e_state(thermal)
    sion_rate = Sion_rate(thermal)
    srec_rate = Srec_rate(thermal)
    # In the two-temperature branch p0 is the sum of ion and electron
    # pressures; the single-temperature branch uses the unified T field.
    if with_TiTe:
        pressure = rho * (Ti + Te)
        pressure_tangent = rho * (Ti + Te)
        if include_impurities:
            pressure += rhoimp * (
                alpha_i * Ti + alpha_e_pressure(Te) * Te
            )
            pressure_tangent += rhoimp * (
                alpha_i * Ti + alpha_e_pressure(Te)
            )
    else:
        pressure = rho * T
        pressure_tangent = pressure
        if include_impurities:
            pressure += rhoimp * alpha_imp_pressure(T) * T
            pressure_tangent += rhoimp * alpha_imp_pressure(T)
    rho_hat = _rho_hat()
    grad_v = grad(v)
    grad_omega = grad(omega)
    lap_v = _laplacian(v)
    lap_omega = _laplacian(omega)

    # Base perpendicular-momentum terms, grouped as they appear in the
    # weak-form derivation.
    inertia = (
        -sp.Rational(1, 2) * _velocity_norm()
        * (grad_v[0] * dZ(rho_hat) - grad_v[1] * dR(rho_hat))
        * xjac
    )
    advection = -rho_hat * R**2 * omega * element_bracket(v, u)
    magnetic = v * element_bracket(psi, j) - v * F0 / R * dphi(j) * xjac
    pressure_term = R**2 * element_bracket(v, pressure)
    pressure_tangent_term = R**2 * element_bracket(v, pressure_tangent)
    viscosity = (
        -visco(thermal) * R**3 * visco_fact_old
        * dot(grad_v, grad_omega) * xjac
        -2 * visco(thermal) * R**2 * visco_fact_new * omega * dR(v) * xjac
        -visco(thermal) * R * visco_fact_new
        * (dR(v) * dR(dphi(dphi(u))) + dZ(v) * dZ(dphi(dphi(u))))
        * xjac
        -visco_num(thermal) * lap_v * lap_omega * xjac
    )
    B = inertia + advection + magnetic + pressure_term + viscosity
    diamagnetic = sp.S.Zero
    diamagnetic_viscosity = sp.S.Zero
    if neo_only:
        B = sp.S.Zero
    if neutral_only:
        B = sp.S.Zero
    if include_extensions:
        include_diamagnetic = True
        include_conservative = True
        include_tgnum = True
        include_neo = True
    if include_diamagnetic:
        pi = _diamagnetic_pressure(
            with_TiTe=with_TiTe, include_impurities=include_impurities
        )
        W_dia = _diamagnetic_viscosity(
            with_TiTe=with_TiTe, include_impurities=include_impurities
        )
        diamagnetic = (
            -v * tauIC * 2 * R**4 * element_bracket(pi, omega)
            -tauIC * 2 * R**3 * dZ(pi) * dot(grad_v, grad(u)) * xjac
            -v * tauIC * 2 * R**4 * (
                dR(dZ(u)) * (dR(dR(pi)) - dZ(dZ(pi)))
                - dR(dZ(pi)) * (dR(dR(u)) - dZ(dZ(u)))
            ) * xjac
        )
        diamagnetic_viscosity = (
            dvisco_state(thermal) * R * W_dia
            * dot(grad(Ti if with_TiTe else T), grad_v) * xjac
            + visco(thermal) * R * W_dia * lap_v * xjac
        )
        B += diamagnetic + diamagnetic_viscosity
    if include_conservative:
        grad_u = grad(u)
        velocity_contraction = dot(grad_v, grad_u)
        rho_hat_x = dR(rho_hat)
        rho_hat_y = dZ(rho_hat)
        conservative = fact_conservative_u * (
            -R**2 * (rho_hat_x * dZ(u) - rho_hat_y * dR(u))
            * velocity_contraction * xjac
            +R * F0 * (rho * dphi(vpar) + vpar * dphi(rho))
            * velocity_contraction * xjac
            +R**2 * rho * (dR(vpar) * dZ(psi) - dZ(vpar) * dR(psi))
            * velocity_contraction * xjac
            +R**2 * vpar * (dR(rho) * dZ(psi) - dZ(rho) * dR(psi))
            * velocity_contraction * xjac
        )
        B += conservative
    if include_neutrals:
        source_contraction = dot(grad_v, grad(u)) * xjac
        neutral_sources = (
            (rho + alpha_e_value * rhoimp) * rhon * sion_rate
            - (rho + alpha_e_value * rhoimp) * (rho - rhoimp) * srec_rate
        )
    if include_auxiliary:
        # Runaway/auxiliary pressure contributions are supplied as scalar
        # work variables by JOREK and enter the weak form directly.
        B += (
            -R * v * (aux_P_par_re + aux_P_perp_re) * xjac
            +R**2 * (
                -aux_divPIR_perp * dZ(v)
                +aux_divPIZ_perp * dR(v)
            ) * xjac
        )
        B += (
            (1 - delta_n_convection) * R**3
            * neutral_sources * source_contraction
            + (1 - fact_conservative_u) * R**3
            * (particle_source + source_pellet + source_bg_drift + source_imp_drift)
            * source_contraction
        )
    if include_tgnum:
        tstep = coefficient("tstep")
        # The source splits each TGNUM product into one active and one
        # background bracket.  Preserve that split explicitly: differentiating
        # both brackets at once creates spurious products of two trial
        # velocities and loses cross terms such as
        # ``v_y*u0_x*u_y*w0_x``.
        velocity_cross = _poloidal_cross(v, u)
        velocity_cross_background = _poloidal_cross(v, freeze(u))
        background_cross = _poloidal_cross(omega, u)
        background_cross_background = _poloidal_cross(omega, freeze(u))
        rho_background_cross = _poloidal_cross(rho_hat, u)
        rho_background_cross_background = _poloidal_cross(rho_hat, freeze(u))
        tgnum = (
            -tgnum_u * sp.Rational(1, 4) * rho_hat * R**3
            * background_cross * velocity_cross_background * xjac * tstep
            -tgnum_u * sp.Rational(1, 4) * rho_hat * R**3
            * background_cross_background * velocity_cross * xjac * tstep
            -tgnum_u * sp.Rational(1, 4) * omega * R**3
            * rho_background_cross * velocity_cross_background * xjac * tstep
            * fact_conservative_u
            -tgnum_u * sp.Rational(1, 4) * omega * R**3
            * rho_background_cross_background * velocity_cross * xjac * tstep
            * fact_conservative_u
        )
        B += tgnum
    if include_neo:
        grad_psi = grad(psi)
        grad_u = grad(u)
        grad_pi = grad(
            _diamagnetic_pressure(
                with_TiTe=with_TiTe, include_impurities=include_impurities
            )
        )
        grad_ti = grad(Ti)
        btheta2 = (grad_psi[0]**2 + grad_psi[1]**2) / R**2
        neo_force = (
            rho * dot(grad_psi, grad_u)
            +tauIC * 2 * dot(grad_psi, grad_pi)
            +aki_neo_prof * tauIC * 2 * rho * dot(grad_psi, grad_ti)
            -rho * vpar * btheta2
        )
        B += (
            amu_neo_prof * _parallel_norm(psi)
            / (btheta2 + epsil)**2
            * dot(grad_psi, grad_v) * neo_force * R * xjac
        )
    # The element routine uses the corrected density ``r0_corr`` in the
    # velocity mass matrix.  Keep the correction evaluated at the current
    # state, while freezing its density argument for the tangent; this gives
    # the expected corrected-density coefficient in ``amat(var_u,var_u)``
    # without introducing an artificial density derivative there.
    A = -R * _rho_hat(corr_neg_dens(freeze(rho))) * dot(grad(v), grad(u)) * xjac
    # Conservative momentum form contributes an additional mass-like term.
    # It is part of A (and therefore carries the (1+zeta) factor in AMAT),
    # rather than B; omitting it loses the ``rho`` and ``rho_corr`` terms in
    # the u/u and u/rho blocks.
    if include_conservative:
        A += -fact_conservative_u * R**3 * rho * dot(grad(v), grad(u)) * xjac
    # The legacy one-temperature element tangent uses the undoubled
    # diamagnetic contribution for the temperature column, while the residual
    # (and the other columns) retain the physical factor two.
    overrides = {}
    if include_impurities:
        tangent_B = B - pressure_term + pressure_tangent_term
        if with_TiTe:
            overrides[Te] = tangent_B
        else:
            overrides[T] = tangent_B
    if not with_TiTe and include_diamagnetic:
        overrides[T] = overrides.get(T, B) - sp.Rational(1, 2) * (
            diamagnetic
        )
    return EvolutionEquation(
        "model600_momentum",
        v,
        A,
        B,
        kind="evolution",
        amat_variation_overrides=overrides,
    )


vpar_psi_equation_1 = induction_equation_1


def _parallel_norm(flux):
    return (F0**2 + grad(flux)[0]**2 + grad(flux)[1]**2) / R**2
