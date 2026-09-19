"""Model-600 weak equations, extending the model-199 DSL definitions."""

import sympy as sp

from .equations import ConstraintEquation, EvolutionEquation
from .external import external_function
from .model199 import FIELDS as MODEL199_FIELDS
from .model199 import GAMMA, gamma, j, omega, psi, rho, T, u, xjac
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
    # The element routine computes dvisco_num_dT but deliberately does not
    # place that derivative in the momentum AMAT temperature columns.
    derivatives={},
    policy="frozen",
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

# ``alpha_i`` is supplied as an element work value in the two-temperature
# branch.  It has no spatial derivative in the assembled source terms; an
# argument-free frozen external represents that convention more faithfully
# than a plain coefficient under derivative expansion.
alpha_i_state = external_function(
    "alpha_i600", arguments=(), derivatives={}, policy="frozen",
    fortran_name="alpha_i",
)

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
alpha_e_temperature = external_function(
    "alpha_e_temperature600", arguments=("T",), derivatives={"T": "alpha_e_bis"},
    policy="piecewise_active", fortran_name="alpha_e_T",
)
alpha_imp_bis_state = external_function(
    "alpha_imp_bis", arguments=("T",), derivatives={"T": "alpha_imp_tri"},
    policy="piecewise_active", fortran_name="alpha_imp_bis",
)
alpha_imp_temperature = external_function(
    "alpha_imp_temperature600", arguments=("T",),
    derivatives={"T": "alpha_imp_bis"},
    policy="piecewise_active", fortran_name="alpha_imp_T",
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
    arguments=("rho", "Ti"),
    derivatives={
        "rho": "W_dia_rho",
        "Ti": "W_dia_Ti",
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
    if include_pressure_coupling:
        # ``construct_pressure`` builds Pe0 = (r0 + rimp0*alpha_e)*Te0 in both
        # temperature models.  The one-temperature model evolves the total
        # temperature and supplies Te0 = T0/2, so the electron temperature is
        # half of the evolved field there.  alpha_e_temperature carries the
        # closure's d(alpha_e*Te)/dTe = alpha_e_bis convention.
        electron_temperature = Te if with_TiTe else T / 2
        Pe = (
            rho * electron_temperature
            + rhoimp * alpha_e_temperature(electron_temperature)
        )
        diamagnetic_core = (
            -v * tauIC / (corrected_rho * _parallel_norm(psi))
            * F0**2 / R**2 * element_bracket(psi, Pe)
            + v * tauIC / (corrected_rho * _parallel_norm(psi))
            * F0**3 / R**3 * dphi(Pe) * xjac
        )
        B += 2 * diamagnetic_core
    if include_runaway_coupling:
        runaway = -v * eta(thermal, rho, rhoimp) * aux_jre_ind / R * xjac
        B += runaway
    A = v * psi / R * xjac
    # The one-temperature electron temperature Te0 = T0/2 is carried by the
    # residual itself, so the temperature column needs no tangent override.
    return EvolutionEquation("model600_induction", v, A, B)


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
    """Ion pressure used by the diamagnetic momentum terms.

    ``construct_pressure`` builds Pi0 from the ion temperature and the ion
    impurity coefficient in both temperature models.  The one-temperature
    model evolves the total temperature and supplies Ti0 = T0/2, so the ion
    temperature is half of the evolved field there.  ``alpha_i`` carries no
    temperature dependence in either model.
    """

    ion_temperature = Ti if with_TiTe else T / 2
    pressure = rho * ion_temperature
    if include_impurities:
        pressure += rhoimp * alpha_i_state() * ion_temperature
    return pressure


def _diamagnetic_viscosity(*, with_TiTe=False, include_impurities=False):
    """Return the externally supplied model-600 ``W_dia`` coefficient."""

    # W_dia is built from the ion pressure only, so it carries no explicit
    # electron-temperature dependence.
    if with_TiTe:
        return W_dia_two(rho, Ti)
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
                alpha_i_state() * Ti + alpha_e_temperature(Te)
            )
            pressure_tangent += rhoimp * (
                alpha_i_state() * Ti + alpha_e_bis_state(Te)
            )
    else:
        pressure = rho * T
        pressure_tangent = pressure
        if include_impurities:
            pressure += rhoimp * alpha_imp_temperature(T)
            pressure_tangent += rhoimp * alpha_imp_bis_state(T)
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
        def diamagnetic_from(pressure):
            return (
                -v * tauIC * 2 * R**4 * element_bracket(pressure, omega)
                -tauIC * 2 * R**3 * dZ(pressure) * dot(grad_v, grad(u)) * xjac
                -v * tauIC * 2 * R**4 * (
                    dR(dZ(u)) * (dR(dR(pressure)) - dZ(dZ(pressure)))
                    - dR(dZ(pressure)) * (dR(dR(u)) - dZ(dZ(u)))
                ) * xjac
            )

        diamagnetic = diamagnetic_from(pi)
        diamagnetic_viscosity = (
            dvisco_state(thermal) * R * W_dia
            * dot(grad(Ti if with_TiTe else T / 2), grad_v) * xjac
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
        # The Taylor-Galerkin residual is a single product of two poloidal
        # brackets per contribution.  The Fortran Jacobian writes the two
        # summands of the product rule separately, which is exactly what the
        # directional derivative of this compact form returns; writing the
        # residual itself in split form would count it twice and would double
        # the omega and density tangents.
        velocity_cross = _poloidal_cross(v, u)
        tgnum = (
            -tgnum_u * sp.Rational(1, 4) * rho_hat * R**3
            * _poloidal_cross(omega, u) * velocity_cross * xjac * tstep
            -tgnum_u * sp.Rational(1, 4) * omega * R**3
            * _poloidal_cross(rho_hat, u) * velocity_cross * xjac * tstep
            * fact_conservative_u
        )
        B += tgnum
    if include_neo:
        grad_psi = grad(psi)
        grad_u = grad(u)
        # JOREK uses the common-temperature pressure in the one-temperature
        # branch, while the two-temperature branch uses the ion pressure and
        # an additional ion-temperature contribution.  Keeping this choice
        # explicit is important: the NEO tangent is formed from the same
        # branch-specific residual, not from a universally doubled term.
        if with_TiTe:
            # The NEO work variable ``Pi0`` in JOREK is the main-ion pressure
            # ``r0*Ti0``.  Impurity pressure enters the separate impurity
            # extensions, not this NEO force, even when those extensions are
            # enabled for the surrounding momentum equation.
            grad_pi = grad(rho * Ti)
            grad_ti = grad(Ti)
            neo_temperature_factor = sp.Integer(2)
        else:
            grad_pi = grad(rho * T)
            grad_ti = grad(T)
            # The residual uses the same factor-two diamagnetic pressure
            # convention as the source NEO RHS.  The one-temperature AMAT
            # override below adjusts only the temperature tangent to the
            # legacy undoubled form used by the element routine.
            neo_temperature_factor = sp.Integer(2)
        btheta2 = (grad_psi[0]**2 + grad_psi[1]**2) / R**2
        neo_force = (
            rho * dot(grad_psi, grad_u)
            +tauIC * neo_temperature_factor * dot(grad_psi, grad_pi)
            +aki_neo_prof * tauIC * neo_temperature_factor
            * rho * dot(grad_psi, grad_ti)
            -rho * vpar * btheta2
        )
        neo_contribution = (
            # ``BB2`` is a work variable in the element routine and is frozen
            # in the Newton tangent.  The poloidal ``Btheta2`` factor is
            # differentiated separately below, but the numerator must not
            # contribute a ``BB2_psi`` variation.
            amu_neo_prof * _parallel_norm(freeze(psi))
            / (btheta2 + epsil)**2
            * dot(grad_psi, grad_v) * neo_force * R * xjac
        )
        B += neo_contribution
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
        # JOREK's conservative history contribution is evaluated with the
        # background velocity.  It contributes to the rho-history tangent,
        # but not to amat(var_u,var_u).
        A += (
            -fact_conservative_u * R**3 * rho
            * dot(grad(v), grad(freeze(u))) * xjac
        )
    # The one-temperature ion temperature Ti0 = T0/2 is now carried by the
    # residual itself, so the temperature column needs no tangent override.
    overrides = {}
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


# ---------------------------------------------------------------------------
# Density equation (``var_rho``)
# ---------------------------------------------------------------------------

# Transport coefficients supplied as element work values.  They are profile
# quantities evaluated at the background state and are frozen in the Newton
# tangent, exactly as the element routine treats them.
D_par_local = coefficient("D_par_local")
D_par_local_imp = coefficient("D_par_local_imp")
D_par_sc_num = coefficient("D_par_sc_num")
D_par_imp_sc_num = coefficient("D_par_imp_sc_num")
D_perp_num_psin = coefficient("D_perp_num_psin")
D_prof = coefficient("D_prof")
D_prof_imp = coefficient("D_prof_imp")
tau_sc = coefficient("tau_sc")
tgnum_rho = coefficient("tgnum_rho")
V_prof_pinch = coefficient("V_prof_pinch")
aux_rho0 = coefficient("aux_rho0")


def _b_dot_grad(value, flux=psi):
    """Return the element routine's ``Bgrad_*`` work value for ``value``.

    ``Bgrad_rho = (F0/BigR*r0_p + r0_x*ps0_y - r0_y*ps0_x)/BigR``.  The element
    routine stores the toroidal and poloidal halves of the *test* function
    version separately (``Bgrad_rho_k_star`` and ``Bgrad_rho_star``) because
    they belong to different FFT channels; the channel split is recovered
    automatically from the toroidal derivative order.
    """

    return (
        F0 / R * dphi(value)
        + dR(value) * dZ(flux) - dZ(value) * dR(flux)
    ) / R


def _b_dot_grad_element(value, flux=psi):
    """``_b_dot_grad`` written with the element-coordinate poloidal bracket.

    The two forms are identical, because ``a_s b_t - a_t b_s`` equals
    ``xjac*(a_x b_y - a_y b_x)``.  Which one to use is purely a matter of
    matching the spelling of the element routine, which writes the parallel
    pressure gradient, the parallel kinetic-energy flux and the ``tgnum_vpar``
    terms in element coordinates.
    """

    return (
        F0 / R * dphi(value) - element_bracket(flux, value) / xjac
    ) / R


def _perpendicular_diffusion(test, value):
    """Weak perpendicular diffusion operator, including its toroidal part."""

    return dot(grad(test), grad(value)) + dphi(test) * dphi(value) / R**2


def density_equation_rho(
    *,
    with_TiTe=False,
    include_diamagnetic=True,
    include_parallel_velocity=True,
    include_neutrals=True,
    include_impurities=True,
    include_tgnum=True,
    include_pinch=True,
    include_auxiliary=True,
):
    """Return the model-600 density equation for ``var_rho``."""

    v = test_function("v")
    thermal = Te if with_TiTe else T
    alpha_e_value = alpha_e_state(thermal)
    sion_rate = Sion_rate(thermal)
    srec_rate = Srec_rate(thermal)
    bb2 = _parallel_norm(psi)
    # ``(D_par_local + D_par_sc_num*tau_sc) - D_prof`` is the parallel
    # diffusivity in excess of the perpendicular one; the element routine
    # writes the same grouping for the impurity species.
    d_par_excess = D_par_local + D_par_sc_num * tau_sc - D_prof
    d_par_excess_imp = D_par_local_imp + D_par_imp_sc_num * tau_sc - D_prof_imp

    B = (
        v * R * (
            particle_source + source_pellet + source_bg_drift + source_imp_drift
        ) * xjac
        + v * R**2 * element_bracket(rho, u)
        + v * 2 * R * rho * dZ(u) * xjac
        - d_par_excess * R / bb2 * _b_dot_grad(v)
        * (_b_dot_grad(rho) - _b_dot_grad(rhoimp)) * xjac
        - d_par_excess_imp * R / bb2 * _b_dot_grad(v)
        * _b_dot_grad(rhoimp) * xjac
        - D_prof * R * _perpendicular_diffusion(v, rho - rhoimp) * xjac
        - D_prof_imp * R * _perpendicular_diffusion(v, rhoimp) * xjac
        - D_perp_num_psin * _laplacian(v) * _laplacian(rho) * R * xjac
    )
    if include_parallel_velocity:
        B += (
            -v * F0 / R * vpar * dphi(rho) * xjac
            - v * vpar * element_bracket(rho, psi)
            - v * F0 / R * rho * dphi(vpar) * xjac
            - v * rho * element_bracket(vpar, psi)
        )
    if include_diamagnetic:
        pi = _diamagnetic_pressure(
            with_TiTe=with_TiTe, include_impurities=include_impurities
        )
        B += v * 2 * tauIC * 2 * dZ(pi) * R * xjac
    if include_neutrals:
        electron_density = rho + alpha_e_value * rhoimp
        B += (
            v * electron_density * rhon * R * sion_rate * xjac
            - v * electron_density * (rho - rhoimp) * R * srec_rate * xjac
        )
    if include_tgnum:
        timestep = coefficient("tstep")
        B += (
            -tgnum_rho * sp.Rational(1, 4) * R**3
            * _poloidal_cross(rho, u) * _poloidal_cross(v, u) * xjac * timestep
            - tgnum_rho * sp.Rational(1, 4) * R * vpar**2
            * _b_dot_grad(rho) * _b_dot_grad(v) * xjac * timestep
        )
    if include_auxiliary:
        B += v * R * aux_rho0 * xjac
    if include_pinch:
        # Weak form of -div(rho * V_pinch), with
        # V_pinch = -V_prof_pinch * grad(psi)/|grad(psi)|.  The pinch
        # direction is built from the evolved flux, so it is varied like any
        # other state quantity; the element routine freezes it (see
        # JOREK_FINDINGS.md).
        psi_gradient = grad(psi)
        B += (
            -V_prof_pinch
            / sp.sqrt(psi_gradient[0]**2 + psi_gradient[1]**2)
            * dot(grad(v), psi_gradient) * rho * R * xjac
        )
    A = v * rho * R * xjac
    return EvolutionEquation("model600_density", v, A, B, kind="evolution")


# ---------------------------------------------------------------------------
# Parallel velocity equation (``var_vpar``)
# ---------------------------------------------------------------------------

visco_par = coefficient("visco_par")
visco_par_sc_num = coefficient("visco_par_sc_num")
visco_par_num = coefficient("visco_par_num")
visco_par_par = coefficient("visco_par_par")
tgnum_vpar = coefficient("tgnum_vpar")
aux_mom_par0 = coefficient("aux_mom_par0")

# Prescribed rotation profile subtracted from the parallel-velocity gradient
# by the perpendicular parallel viscosity.  It is a flux function whose flux
# derivative the element routine stores as the frozen profile work value
# ``dV_dpsi_source``:  ``Vt0_x = dV_dpsi_source*ps0_x`` and
# ``Vt_x_psi = dV_dpsi_source*psi_x``.  Only the flux is varied, so
# ``grad(Vt) = dV_dpsi_source*grad(psi)`` is the faithful form.
dV_dpsi_source = coefficient("dV_dpsi_source")


def parallel_velocity_equation_vpar(
    *,
    with_TiTe=False,
    element_brackets=True,
    include_conservative=True,
    include_tgnum=True,
    include_neutrals=True,
    include_impurities=True,
    include_sources=True,
    include_pinch=True,
    include_auxiliary=True,
):
    """Return the model-600 parallel-velocity equation for ``var_vpar``.

    The perpendicular parallel viscosity follows the default
    ``normalized_velocity_profile = .true.`` branch of the element routine.

    ``element_brackets`` selects how the poloidal bracket of the parallel
    kinetic-energy flux is spelled.  The two forms are mathematically
    identical, but the element routine writes this particular group in element
    coordinates in ``rhs_ij(var_vpar)`` and in the ``rho`` and ``vpar``
    columns and in physical coordinates in ``amat(var_vpar,var_psi)``, so the
    exporter builds both and aligns each assignment against the spelling it
    happens to use.  The pressure gradient and the Taylor-Galerkin terms are
    written in element coordinates throughout and need no such switch.
    """

    v = test_function("v")
    thermal = Te if with_TiTe else T
    alpha_e_value = alpha_e_state(thermal)
    sion_rate = Sion_rate(thermal)
    srec_rate = Srec_rate(thermal)
    bb2 = _parallel_norm(psi)
    if with_TiTe:
        pressure = rho * (Ti + Te)
        if include_impurities:
            pressure += rhoimp * (
                alpha_i_state() * Ti + alpha_e_temperature(Te)
            )
    else:
        pressure = rho * T
        if include_impurities:
            pressure += rhoimp * alpha_imp_temperature(T)
    rho_hat = _rho_hat()
    visco_par_eff = visco_par + visco_par_sc_num * tau_sc
    kinetic_gradient = (
        _b_dot_grad_element if element_brackets else _b_dot_grad
    )
    rotation_shear = tuple(
        parallel - dV_dpsi_source * flux
        for parallel, flux in zip(grad(vpar), grad(psi))
    )

    B = (
        # Parallel pressure gradient.
        -v * R * _b_dot_grad_element(pressure) * xjac
        # Parallel advection of the kinetic energy, 0.5*v_par**2*B**2.
        + sp.Rational(1, 2) * vpar**2 * bb2 * R * xjac
        * (rho * kinetic_gradient(v) + v * kinetic_gradient(rho))
        # Numerical and physical parallel viscosities.
        - visco_par_num * _laplacian(v) * _laplacian(vpar) * R * xjac
        - visco_par_par * F0**2 / (R * bb2)
        * _b_dot_grad(vpar) * _b_dot_grad(v) * xjac
        - visco_par_eff * dot(grad(v), rotation_shear) * R * xjac
    )
    if include_sources:
        B += (
            -v * (
                particle_source + source_pellet + source_bg_drift
                + source_imp_drift
            ) * vpar * bb2 * R * xjac * (1 - fact_conservative_u)
        )
    if include_conservative:
        # -(d_t rho + div(rho v)) v_par B**2 R; the d_t rho part belongs to A.
        B += fact_conservative_u * v * vpar * bb2 * xjac * (
            _poloidal_cross(rho_hat, u) - R * _b_dot_grad(rho * vpar)
        )
    if include_neutrals:
        electron_density = rho + alpha_e_value * rhoimp
        B += (1 - delta_n_convection) * (
            -v * electron_density * rhon * sion_rate * vpar * bb2 * R * xjac
            + v * electron_density * (rho - rhoimp) * srec_rate
            * vpar * bb2 * R * xjac
        )
    if include_tgnum:
        timestep = coefficient("tstep")
        B += (
            -tgnum_vpar * sp.Rational(1, 4) * rho * vpar**2 * bb2 * R
            * _b_dot_grad_element(vpar) * _b_dot_grad_element(v) * xjac * timestep
            - tgnum_vpar * sp.Rational(1, 4) * v * vpar**2 * bb2 * R
            * (1 - fact_conservative_u)
            * _b_dot_grad_element(vpar) * _b_dot_grad_element(rho)
            * xjac * timestep
            - tgnum_vpar * sp.Rational(1, 4) * vpar**3 * bb2 * R
            * fact_conservative_u
            * _b_dot_grad_element(rho) * _b_dot_grad_element(v) * xjac * timestep
        )
    if include_auxiliary:
        B += (
            -v * aux_rho0 * vpar * bb2 * R * (1 - fact_conservative_u) * xjac
            + v * R * aux_mom_par0 * xjac
        )
    if include_pinch:
        psi_gradient = grad(psi)
        B += (
            V_prof_pinch
            / sp.sqrt(psi_gradient[0]**2 + psi_gradient[1]**2)
            * dot(psi_gradient, grad(vpar)) * rho * v * R * xjac
        )
    # The parallel momentum density is rho*v_par*B**2*R.  The element routine
    # freezes the density correction in the tangent, as in the perpendicular
    # momentum equation.
    A = v * corr_neg_dens(freeze(rho)) * vpar * bb2 * R * xjac
    if include_conservative:
        A += fact_conservative_u * v * rho * freeze(vpar) * bb2 * R * xjac
    return EvolutionEquation("model600_parallel_velocity", v, A, B)


# ---------------------------------------------------------------------------
# Impurity density equation (``var_rhoimp``)
# ---------------------------------------------------------------------------

Dn_perp_num = coefficient("Dn_perp_num")
tgnum_rhoimp = coefficient("tgnum_rhoimp")


def impurity_density_equation_rhoimp(
    *,
    with_TiTe=False,
    include_parallel_velocity=True,
    include_tgnum=True,
    include_auxiliary=True,
):
    """Return the model-600 impurity-density equation for ``var_rhoimp``.

    The impurity species is advected and diffused like the main one, but it
    has no atomic sources and no diamagnetic contribution; the element routine
    marks the latter with an explicit placeholder comment.
    """

    v = test_function("v")
    bb2 = _parallel_norm(psi)
    d_par_excess_imp = D_par_local_imp + D_par_imp_sc_num * tau_sc - D_prof_imp

    B = (
        -d_par_excess_imp * R / bb2 * _b_dot_grad(v)
        * _b_dot_grad(rhoimp) * xjac
        - D_prof_imp * R * _perpendicular_diffusion(v, rhoimp) * xjac
        + v * R**2 * element_bracket(rhoimp, u)
        + v * 2 * R * rhoimp * dZ(u) * xjac
        - Dn_perp_num * _laplacian(v) * _laplacian(rhoimp) * R * xjac
    )
    if include_parallel_velocity:
        B += (
            -v * F0 / R * vpar * dphi(rhoimp) * xjac
            - v * vpar * element_bracket(rhoimp, psi)
            - v * F0 / R * rhoimp * dphi(vpar) * xjac
            - v * rhoimp * element_bracket(vpar, psi)
        )
    if include_tgnum:
        timestep = coefficient("tstep")
        B += (
            -tgnum_rhoimp * sp.Rational(1, 4) * R**3
            * _poloidal_cross(rhoimp, u) * _poloidal_cross(v, u)
            * xjac * timestep
            - tgnum_rhoimp * sp.Rational(1, 4) * R * vpar**2
            * _b_dot_grad(rhoimp) * _b_dot_grad(v) * xjac * timestep
        )
    if include_auxiliary:
        B += R * v * source_imp_drift * xjac
    A = v * rhoimp * R * xjac
    return EvolutionEquation("model600_impurity_density", v, A, B)


# ---------------------------------------------------------------------------
# Ion energy equation (``var_Ti``)
# ---------------------------------------------------------------------------

heat_source_i = coefficient("heat_source_i")
ZK_i_perp_num_psin = coefficient("ZK_i_perp_num_psin")
tgnum_Ti = coefficient("tgnum_Ti")
visco_par_heating = coefficient("visco_par_heating")
implicit_heat_source = coefficient("implicit_heat_source")
Tie_min_neg = coefficient("Tie_min_neg")
aux_E0_Ti = coefficient("aux_E0_Ti")

# Parallel ion heat conductivity and the perpendicular profile.  The element
# routine supplies their temperature and density derivatives as work values.
ZKi_par = external_function(
    "ZKi_par600", arguments=("Ti",), derivatives={"Ti": "dZKi_par_dT"},
    policy="piecewise_active", fortran_name="ZKi_par_T",
)
ZKi_perp = external_function(
    "ZKi_prof600", arguments=("rho",), derivatives={"rho": "dZKi_prof_drho"},
    policy="piecewise_active", fortran_name="ZKi_prof",
)
# The viscous-heating coefficient is evaluated at the electron temperature.
visco_heating = external_function(
    "visco_heating600", arguments=("Te",),
    derivatives={"Te": "dvisco_dT_heating"},
    policy="piecewise_active", fortran_name="visco_T_heating",
)
# Ion-electron energy exchange.
Ti_e_exchange = external_function(
    "dTi_e600", arguments=("Ti", "Te", "rho", "rhoimp"),
    derivatives={
        "Ti": "ddTi_e_dTi", "Te": "ddTi_e_dTe",
        "rho": "ddTi_e_drho", "rhoimp": "ddTi_e_drhoimp",
    },
    policy="piecewise_active", fortran_name="dTi_e",
)
# Implicit heating floor for small temperatures.  The element routine writes
# ``min(Ti0,Tie_min_neg)`` and the exponential of its normalized distance to
# the floor inline; both are registered here so the tangent stays well
# defined.  ``_normalize_model600_text`` maps the two spellings onto these
# names and resolves their derivatives.
Ti_floor = external_function(
    "Ti_floor600", arguments=("Ti",), derivatives={"Ti": "dTi_floor"},
    policy="piecewise_active", fortran_name="Ti0_floor",
)
Ti_floor_exp = external_function(
    "Ti_floor_exp600", arguments=("Ti",),
    derivatives={"Ti": "dTi_floor_exp"},
    policy="piecewise_active", fortran_name="Ti_floor_exp",
)
# Negative-density correction of the impurity species.
corr_neg_dens_imp = external_function(
    "corr_neg_dens_imp", arguments=("rhoimp",),
    derivatives={"rhoimp": "drimp0_corr_dn"},
    policy="piecewise_active", fortran_name="corr_neg_dens_imp",
)


def ion_energy_equation_Ti(
    *,
    element_brackets=True,
    include_impurities=True,
    include_parallel_velocity=True,
    include_tgnum=True,
    include_friction=True,
    include_heating=True,
    include_auxiliary=True,
):
    """Return the model-600 ion energy equation for ``var_Ti``.

    ``element_brackets`` selects the spelling of the poloidal advection
    bracket, as in the parallel-velocity equation: the element routine writes
    it in element coordinates in the residual and in the ``rho``, ``Ti`` and
    ``rhoimp`` columns, and in physical coordinates in ``amat(var_Ti,var_u)``.
    """

    v = test_function("v")
    alpha_i = alpha_i_state()
    ion_density = rho + alpha_i * rhoimp if include_impurities else rho
    corrected_ion_density = (
        corr_neg_dens(rho) + alpha_i * corr_neg_dens_imp(rhoimp)
        if include_impurities else corr_neg_dens(rho)
    )
    electron_density = (
        rho + alpha_e_state(Te) * rhoimp if include_impurities else rho
    )
    ion_pressure = ion_density * Ti
    bb2 = _parallel_norm(psi)
    conduction_excess = ZKi_par(Ti) - ZKi_perp(rho)
    advection_bracket = element_bracket if element_brackets else (
        lambda left, right: xjac * _poloidal_cross(left, right)
    )

    B = (
        v * R * heat_source_i * xjac
        # Advection and compression of the ion pressure.
        + v * R**2 * advection_bracket(ion_pressure, u)
        + 2 * GAMMA * v * R * ion_pressure * dZ(u) * xjac
        # Perpendicular and parallel heat conduction.
        - conduction_excess * R / bb2
        * _b_dot_grad(v) * _b_dot_grad(Ti) * xjac
        - ZKi_perp(rho) * R * _perpendicular_diffusion(v, Ti) * xjac
        - ZK_i_perp_num_psin * _laplacian(v) * _laplacian(Ti) * R * xjac
    )
    if include_parallel_velocity:
        B += (
            -v * F0 / R * vpar * dphi(ion_pressure) * xjac
            - v * vpar * element_bracket(ion_pressure, psi)
            - GAMMA * v * ion_pressure * (
                element_bracket(vpar, psi) + F0 / R * dphi(vpar) * xjac
            )
        )
    if include_tgnum:
        timestep = coefficient("tstep")
        B += (
            -tgnum_Ti * sp.Rational(1, 4) * R**3
            * _poloidal_cross(ion_pressure, u) * _poloidal_cross(v, u)
            * xjac * timestep
            - tgnum_Ti * sp.Rational(1, 4) * R * vpar**2
            * _b_dot_grad(ion_pressure) * _b_dot_grad(v) * xjac * timestep
        )
    if include_friction:
        # Kinetic energy released by the particle sources.
        released = (
            electron_density * rhon * Sion_rate(Te)
            + particle_source + source_pellet
            + source_bg_drift + source_imp_drift
        )
        B += (
            v * R * (GAMMA - 1) / 2
            * (vpar**2 * bb2 + _velocity_norm(u)) * released * xjac
        )
    if include_heating:
        grad_vpar = grad(vpar)
        B += (GAMMA - 1) * R * visco_par_heating * xjac * (
            v * dot(grad_vpar, grad_vpar) + vpar * dot(grad(v), grad_vpar)
        )
        heating = visco_heating(Te)
        B += (
            -(GAMMA - 1) * v * heating * R**3 * visco_fact_old
            * dot(grad(u), grad(omega)) * xjac
            - (GAMMA - 1) * v * heating * 2 * R**2 * visco_fact_new
            * omega * dR(u) * xjac
            - (GAMMA - 1) * v * heating * R * visco_fact_new
            * (
                dR(u) * dR(dphi(dphi(u))) + dZ(u) * dZ(dphi(dphi(u)))
            ) * xjac
        )
    B += (
        # Ion-electron energy exchange.
        v * R * Ti_e_exchange(Ti, Te, rho, rhoimp) * xjac
        # Implicit heating floor for small temperatures.
        + implicit_heat_source * (gamma - 1) * v * R * xjac * (
            sp.Rational(1, 2) * Tie_min_neg * (1 + Ti_floor_exp(Ti))
            - Ti_floor(Ti)
        )
        # Recombination sink.
        - v * Ti * R * corr_neg_dens(rho)**2 * Srec_rate(Te) * xjac
    )
    if include_auxiliary:
        B += (
            v * R * aux_E0_Ti * xjac
            + (gamma - 1) * sp.Rational(1, 2) * v * aux_rho0
            * vpar**2 * bb2 * R * xjac
            - (gamma - 1) * v * aux_mom_par0 * vpar * R * xjac
        )
    A = v * corrected_ion_density * Ti * R * xjac
    return EvolutionEquation("model600_ion_energy", v, A, B)
