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


# -------------------------------------------------------------------------
# Fields
# -------------------------------------------------------------------------

vpar = field("vpar", fortran_current="vpar0", fortran_trial="vpar")

Ti = field("Ti", fortran_current="Ti0", fortran_trial="Ti")

Te = field("Te", fortran_current="Te0", fortran_trial="Te")

rhon = field("rhon", fortran_current="rn0", fortran_trial="rhon")

rhoimp = field("rhoimp", fortran_current="rimp0", fortran_trial="rhoimp")

FIELDS = MODEL199_FIELDS + (vpar, Ti, Te, rhon, rhoimp)

# -------------------------------------------------------------------------
# Scalar work values
# -------------------------------------------------------------------------
# Values the element routine computes once per quadrature point and
# freezes in the Newton tangent.
Jb = coefficient("Jb")
aux_jre_ind = coefficient("aux_jre_ind")
tauIC = coefficient("tauIC")
tstep = coefficient("tstep")
r0_corr = coefficient("r0_corr")
F0 = coefficient("F0")

factor_psi = tuple(coefficient("factor_psi_{}".format(index)) for index in range(1, 7))

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
visco_par = coefficient("visco_par")
visco_par_sc_num = coefficient("visco_par_sc_num")
visco_par_num = coefficient("visco_par_num")
visco_par_par = coefficient("visco_par_par")
tgnum_vpar = coefficient("tgnum_vpar")
aux_mom_par0 = coefficient("aux_mom_par0")
dV_dpsi_source = coefficient("dV_dpsi_source")
Dn_perp_num = coefficient("Dn_perp_num")
tgnum_rhoimp = coefficient("tgnum_rhoimp")
heat_source_i = coefficient("heat_source_i")
ZK_i_perp_num_psin = coefficient("ZK_i_perp_num_psin")
tgnum_Ti = coefficient("tgnum_Ti")
visco_par_heating = coefficient("visco_par_heating")
implicit_heat_source = coefficient("implicit_heat_source")
Tie_min_neg = coefficient("Tie_min_neg")
aux_E0_Ti = coefficient("aux_E0_Ti")
heat_source_e = coefficient("heat_source_e")
ZK_e_perp_num_psin = coefficient("ZK_e_perp_num_psin")
tgnum_Te = coefficient("tgnum_Te")
ksi_ion_norm = coefficient("ksi_ion_norm")
aux_jre = coefficient("aux_jre")
aux_E0_Te = coefficient("aux_E0_Te")
power_dens_teleport_ju = coefficient("power_dens_teleport_ju")
heat_source_total = coefficient("heat_source")
ZK_perp_num_psin = coefficient("ZK_perp_num_psin")
tgnum_T = coefficient("tgnum_T")
T_min_neg = coefficient("T_min_neg")
aux_E0 = coefficient("aux_E0")


# -------------------------------------------------------------------------
# State-dependent work values
# -------------------------------------------------------------------------
# Quantities the element routine supplies together with the derivatives
# the Jacobian needs; declaring the interface lets the generator
# differentiate them by name.
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

alpha_i_state = external_function(
    "alpha_i600", arguments=(), derivatives={}, policy="frozen",
    fortran_name="alpha_i",
)

dvisco_state = external_function(
    "dvisco_state600", arguments=("T",),
    derivatives={"T": "d2visco_dT2"},
    policy="piecewise_active", fortran_name="dvisco_dT",
)

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

corr_neg_dens = external_function(
    "corr_neg_dens",
    arguments=("rho",),
    derivatives={"rho": "dr0_corr_dn"},
    policy="piecewise_active",
    fortran_name="corr_neg_dens",
)

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

ZKi_par = external_function(
    "ZKi_par600", arguments=("Ti",), derivatives={"Ti": "dZKi_par_dT"},
    policy="piecewise_active", fortran_name="ZKi_par_T",
)

ZKi_perp = external_function(
    "ZKi_prof600", arguments=("rho",), derivatives={"rho": "dZKi_prof_drho"},
    policy="piecewise_active", fortran_name="ZKi_prof",
)

visco_heating = external_function(
    "visco_heating600", arguments=("Te",),
    derivatives={"Te": "dvisco_dT_heating"},
    policy="piecewise_active", fortran_name="visco_T_heating",
)

Ti_e_exchange = external_function(
    "dTi_e600", arguments=("Ti", "Te", "rho", "rhoimp"),
    derivatives={
        "Ti": "ddTi_e_dTi", "Te": "ddTi_e_dTe",
        "rho": "ddTi_e_drho", "rhoimp": "ddTi_e_drhoimp",
    },
    policy="piecewise_active", fortran_name="dTi_e",
)

Ti_floor = external_function(
    "Ti_floor600", arguments=("Ti",), derivatives={"Ti": "dTi_floor"},
    policy="piecewise_active", fortran_name="Ti0_floor",
)

Ti_floor_exp = external_function(
    "Ti_floor_exp600", arguments=("Ti",),
    derivatives={"Ti": "dTi_floor_exp"},
    policy="piecewise_active", fortran_name="Ti_floor_exp",
)

corr_neg_dens_imp = external_function(
    "corr_neg_dens_imp", arguments=("rhoimp",),
    derivatives={"rhoimp": "drimp0_corr_dn"},
    policy="piecewise_active", fortran_name="corr_neg_dens_imp",
)

E_ion_bg_state = external_function(
    "E_ion_bg600", arguments=(), derivatives={}, policy="frozen",
    fortran_name="E_ion_bg",
)

ZKe_par = external_function(
    "ZKe_par600", arguments=("Te",), derivatives={"Te": "dZKe_par_dT"},
    policy="piecewise_active", fortran_name="ZKe_par_T",
)

ZKe_perp = external_function(
    "ZKe_prof600", arguments=("rho",), derivatives={"rho": "dZKe_prof_drho"},
    policy="piecewise_active", fortran_name="ZKe_prof",
)

eta_ohm_e = external_function(
    "eta_ohm600", arguments=("Te", "rho", "rhoimp"),
    derivatives={
        "Te": "deta_dT_ohm", "rho": "deta_dr0_ohm",
        "rhoimp": "deta_drimp0_ohm",
    },
    policy="piecewise_active", fortran_name="eta_T_ohm",
)

dE_ion_dT_state = external_function(
    "dE_ion_dT", arguments=("Te",), derivatives={}, policy="frozen",
)

E_ion_state = external_function(
    "E_ion600", arguments=("Te",), derivatives={"Te": "dE_ion_dT"},
    policy="piecewise_active", fortran_name="E_ion",
)

LradDrays = external_function(
    "LradDrays600", arguments=("Te",), derivatives={"Te": "dLradDrays_dT"},
    policy="piecewise_active", fortran_name="LradDrays_T",
)

LradDcont = external_function(
    "LradDcont600", arguments=("Te",),
    derivatives={"Te": "dLradDcont_dT_corr"},
    policy="piecewise_active", fortran_name="LradDcont_corr",
)

frad_bg_state = external_function(
    "frad_bg600", arguments=("Te",), derivatives={"Te": "dfrad_bg_dT"},
    policy="piecewise_active", fortran_name="frad_bg",
)

Lrad_state = external_function(
    "Lrad600", arguments=("Te",), derivatives={"Te": "dLrad_dT"},
    policy="piecewise_active", fortran_name="Lrad",
)

Te_i_exchange = external_function(
    "dTe_i600", arguments=("Ti", "Te", "rho", "rhoimp"),
    derivatives={
        "Ti": "ddTe_i_dTi", "Te": "ddTe_i_dTe",
        "rho": "ddTe_i_drho", "rhoimp": "ddTe_i_drhoimp",
    },
    policy="piecewise_active", fortran_name="dTe_i",
)

Te_floor = external_function(
    "Te_floor600", arguments=("Te",), derivatives={"Te": "dTe_floor"},
    policy="piecewise_active", fortran_name="Te0_floor",
)

Te_floor_exp = external_function(
    "Te_floor_exp600", arguments=("Te",),
    derivatives={"Te": "dTe_floor_exp"},
    policy="piecewise_active", fortran_name="Te_floor_exp",
)

corr_neg_dens_n = external_function(
    "corr_neg_dens_n", arguments=("rhon",),
    derivatives={"rhon": "drn0_corr_dn"},
    policy="piecewise_active", fortran_name="corr_neg_dens_n",
)

ZK_par = external_function(
    "ZK_par600", arguments=("T",), derivatives={"T": "dZK_par_dT"},
    policy="piecewise_active", fortran_name="ZK_par_T",
)

ZK_perp = external_function(
    "ZK_prof600", arguments=("rho",), derivatives={"rho": "dZK_prof_drho"},
    policy="piecewise_active", fortran_name="ZK_prof",
)

T_floor = external_function(
    "T_floor600", arguments=("T",), derivatives={"T": "dT_floor"},
    policy="piecewise_active", fortran_name="T0_floor",
)

T_floor_exp = external_function(
    "T_floor_exp600", arguments=("T",), derivatives={"T": "dT_floor_exp"},
    policy="piecewise_active", fortran_name="T_floor_exp",
)


# -------------------------------------------------------------------------
# Operators and shared term groups
# -------------------------------------------------------------------------
def _laplacian(value):
    return dR(dR(value)) + dZ(dZ(value)) + dR(value) / R

def _rho_hat(value=rho):
    return R**2 * value

def _velocity_norm(value=u):
    return R**2 * (dR(value)**2 + dZ(value)**2)

def _poloidal_cross(left, right):
    return dR(left) * dZ(right) - dZ(left) * dR(right)

def _B2(flux):
    """Total magnetic field squared (B^2)"""
    return (F0**2 + grad(flux)[0]**2 + grad(flux)[1]**2) / R**2

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

def _species_advection(v, density, *, bracket=element_bracket):
    """Advection and compression of a particle density."""

    return (
        v * R**2 * bracket(density, u)
        + v * 2 * R * density * dZ(u) * xjac
    )

def _parallel_convection(v, density):
    """Convection of a particle density along the field by ``vpar``."""

    return (
        -v * F0 / R * vpar * dphi(density) * xjac
        - v * vpar * element_bracket(density, psi)
        - v * F0 / R * density * dphi(vpar) * xjac
        - v * density * element_bracket(vpar, psi)
    )

def _species_tgnum(v, density, factor):
    """Taylor-Galerkin stabilization of a particle-density equation."""

    return (
        -factor * sp.Rational(1, 4) * R**3
        * _poloidal_cross(density, u) * _poloidal_cross(v, u) * xjac * tstep
        - factor * sp.Rational(1, 4) * R * vpar**2
        * _b_dot_grad(density) * _b_dot_grad(v) * xjac * tstep
    )

def _parallel_diffusion(v, density, excess):
    """Anisotropic diffusion of a particle density along the field."""

    return (
        -excess * R / _B2(psi)
        * _b_dot_grad(v) * _b_dot_grad(density) * xjac
    )

def _diamagnetic_pressure(*, with_TiTe=False):
    """Ion pressure ``Pi = (rho + alpha_i*rhoimp)*Ti`` of ``construct_pressure``.

    The one-temperature model evolves the total temperature and supplies
    Ti0 = T0/2, so the ion temperature is half of the evolved field there.
    ``alpha_i`` carries no temperature dependence in either model.
    """

    ion_temperature = Ti if with_TiTe else T / 2
    return (rho + rhoimp * alpha_i_state()) * ion_temperature

def _diamagnetic_viscosity(*, with_TiTe=False):
    """The element routine's ``W_dia``.

    It is built from the ion pressure alone, so it carries no explicit
    electron-temperature dependence.
    """

    return W_dia_two(rho, Ti) if with_TiTe else W_dia_single(rho, T)


def _physical_bracket(left, right):
    """The poloidal bracket written in (R,Z) rather than element coordinates."""

    return xjac * _poloidal_cross(left, right)


def _parallel_divergence(value, flux=psi):
    """``xjac*R`` times the parallel gradient, as the energy equations spell it."""

    return element_bracket(value, flux) + F0 / R * dphi(value) * xjac


def _pressure_transport(v, pressure, *, bracket):
    """Advection, compression and parallel convection of a pressure."""

    return (
        v * R**2 * bracket(pressure, u)
        + 2 * GAMMA * v * R * pressure * dZ(u) * xjac
        - v * F0 / R * vpar * dphi(pressure) * xjac
        - v * vpar * element_bracket(pressure, psi)
        - GAMMA * v * pressure * _parallel_divergence(vpar)
    )


def _heat_conduction(v, temperature, parallel, perpendicular, numerical):
    """Parallel, perpendicular and numerical conduction of one temperature."""

    return (
        -(parallel - perpendicular) * R / _B2(psi)
        * _b_dot_grad(v) * _b_dot_grad(temperature) * xjac
        - perpendicular * R * _perpendicular_diffusion(v, temperature) * xjac
        - numerical * _laplacian(v) * _laplacian(temperature) * R * xjac
    )


def _energy_tgnum(v, pressure, factor):
    """Taylor-Galerkin stabilization of an energy equation."""

    return (
        -factor * sp.Rational(1, 4) * R**3
        * _poloidal_cross(pressure, u) * _poloidal_cross(v, u) * xjac * tstep
        - factor * sp.Rational(1, 4) * R * vpar**2
        * _b_dot_grad(pressure) * _b_dot_grad(v) * xjac * tstep
    )


def _heating_floor(v, exponential, floor, minimum):
    """Implicit heating that keeps a temperature away from its floor."""

    return implicit_heat_source * (gamma - 1) * v * R * xjac * (
        sp.Rational(1, 2) * minimum * (1 + exponential) - floor
    )


def _released_kinetic_energy(thermal):
    """Particle sources whose kinetic energy is released into the ions."""

    return (
        (rho + alpha_e_state(thermal) * rhoimp) * rhon * Sion_rate(thermal)
        + particle_source + source_pellet + source_bg_drift + source_imp_drift
    )


def _friction_heating(v, released):
    """Kinetic energy handed to the ions by the particle sources."""

    return (
        v * R * (GAMMA - 1) / 2
        * (vpar**2 * _B2(psi) + _velocity_norm(u)) * released * xjac
    )


def _viscous_heating(v, heating):
    """Parallel and perpendicular viscous heating."""

    grad_vpar = grad(vpar)
    return (
        (GAMMA - 1) * R * visco_par_heating * xjac * (
            v * dot(grad_vpar, grad_vpar) + vpar * dot(grad(v), grad_vpar)
        )
        - (GAMMA - 1) * v * heating * R**3 * visco_fact_old
        * dot(grad(u), grad(omega)) * xjac
        - (GAMMA - 1) * v * heating * 2 * R**2 * visco_fact_new
        * omega * dR(u) * xjac
        - (GAMMA - 1) * v * heating * R * visco_fact_new
        * (dR(u) * dR(dphi(dphi(u))) + dZ(u) * dZ(dphi(dphi(u)))) * xjac
    )


def _kinetic_coupling(v):
    """Energy and momentum handed over by the kinetic neutral/impurity model."""

    return (
        (gamma - 1) * sp.Rational(1, 2) * v * aux_rho0
        * vpar**2 * _B2(psi) * R * xjac
        - (gamma - 1) * v * aux_mom_par0 * vpar * R * xjac
    )


def _ohmic_heating(v, resistivity):
    """Ohmic dissipation of the toroidal current."""

    return v * (GAMMA - 1) * resistivity * ((j - aux_jre) / R)**2 * R * xjac


def _radiation_sinks(v, electron_density, temperature):
    """Line, continuum, background and impurity radiation."""

    return (
        -v * R * electron_density * corr_neg_dens_n(rhon)
        * LradDrays(temperature) * xjac
        - v * R * electron_density
        * (corr_neg_dens(rho) - corr_neg_dens_imp(rhoimp))
        * LradDcont(temperature) * xjac
        - v * R * electron_density * frad_bg_state(temperature) * xjac
        - v * R * electron_density * corr_neg_dens_imp(rhoimp)
        * Lrad_state(temperature) * xjac
    )


def _ionization_energy(temperature):
    """Potential energy stored in the impurity and background ionization."""

    return E_ion_state(temperature) * rhoimp + E_ion_bg_state() * (rho - rhoimp)


def _ionization_energy_transport(v, energy, temperature):
    """Transport and diffusive flux of the ionization potential energy.

    The element routine keeps this group in element coordinates even where the
    pressure advection right above it uses physical ones.
    """

    d_par_excess = D_par_local + D_par_sc_num * tau_sc - D_prof
    d_par_excess_imp = D_par_local_imp + D_par_imp_sc_num * tau_sc - D_prof_imp
    bb2 = _B2(psi)
    return (GAMMA - 1) * (
        v * R**2 * element_bracket(energy, u)
        + 2 * v * R * energy * dZ(u) * xjac
        - v * F0 / R * vpar * dphi(energy) * xjac
        - v * vpar * element_bracket(energy, psi)
        - v * energy * _parallel_divergence(vpar)
        - E_ion_state(temperature) * d_par_excess_imp * R / bb2
        * _b_dot_grad(v) * _b_dot_grad(rhoimp) * xjac
        - E_ion_state(temperature) * D_prof_imp * R
        * _perpendicular_diffusion(v, rhoimp) * xjac
        - E_ion_bg_state() * d_par_excess * R / bb2
        * _b_dot_grad(v) * _b_dot_grad(rho - rhoimp) * xjac
        - E_ion_bg_state() * D_prof * R
        * _perpendicular_diffusion(v, rho - rhoimp) * xjac
    )


# -------------------------------------------------------------------------
# Equations
# -------------------------------------------------------------------------
# In the order of the element routine.
def induction_equation_1(
    *,
    include_diamag: bool = True,
    include_runaway_coupling: bool = True,
    with_TiTe: bool = False,
):
    """Weak form of the induction equation for var_psi."""

    v = test_function("v")
    T_or_Te = Te if with_TiTe else T

    # Time derivative
    A = v * psi / R * xjac    

    # RHS
    B = (
        # -B.grad u
        + v * element_bracket(psi, u) - v * F0 / R * dphi(u) * xjac

        # eta*j
        + v * eta(T_or_Te, rho, rhoimp)*(j - coefficient("current_source") - Jb) / R * xjac

        # hyper-resistivity
        + eta_num_T(T_or_Te) * dot(grad(v), grad(j)) * xjac
    )

    if include_runaway_coupling:
        # -eta*j_RE
        B += -v * eta(T_or_Te, rho, rhoimp) * aux_jre_ind / R * xjac


    if include_diamag:
        # In model600 ``r0_corr`` is the corrected density > 0, to avoid
        # divergence in diamag term denominator.
        rho_corr = corr_neg_dens(rho)
        # ``construct_pressure`` builds Pe0 = (r0 + rimp0*alpha_e)*Te0 in both
        # temperature models.  The one-temperature model evolves the total
        # temperature and supplies Te0 = T0/2, so the electron temperature is
        # half of the evolved field there.  alpha_e_temperature carries the
        # closure's d(alpha_e*Te)/dTe = alpha_e_bis convention.
        Te_gen = Te if with_TiTe else T / 2

        Pe = rho * Te_gen + rhoimp * alpha_e_temperature(Te_gen)

        # B.grad Pe diamagnetic term
        B += 2*(
            - v * tauIC / (rho_corr * _B2(psi)) * F0**2 / R**2 * element_bracket(psi, Pe)
            + v * tauIC / (rho_corr * _B2(psi)) * F0**3 / R**3 * dphi(Pe) * xjac
        )

    return EvolutionEquation("model600_induction", v, A, B)

def momentum_equation_2(
    *,
    with_TiTe=False,
    include_neo=False,
):
    """Model-600 perpendicular momentum equation (``var_u``).

    ``include_neo`` adds the neoclassical friction; the reports omit it by
    default because the element routine guards it with its own ``NEO``
    switch.
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
        pressure += rhoimp * (alpha_i_state() * Ti + alpha_e_temperature(Te))
        pressure_tangent += rhoimp * (
            alpha_i_state() * Ti + alpha_e_bis_state(Te)
        )
    else:
        pressure = rho * T
        pressure_tangent = pressure
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
    pi = _diamagnetic_pressure(with_TiTe=with_TiTe)
    W_dia = _diamagnetic_viscosity(with_TiTe=with_TiTe)
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
    source_contraction = dot(grad_v, grad(u)) * xjac
    neutral_sources = (
        (rho + alpha_e_value * rhoimp) * rhon * sion_rate
        - (rho + alpha_e_value * rhoimp) * (rho - rhoimp) * srec_rate
    )
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
            amu_neo_prof * _B2(freeze(psi))
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

def density_equation_rho(*, with_TiTe=False):
    """Model-600 density equation (``var_rho``)."""

    v = test_function("v")
    thermal = Te if with_TiTe else T
    electron_density = rho + alpha_e_state(thermal) * rhoimp
    # The parallel diffusivity in excess of the perpendicular one; the element
    # routine writes the same grouping for the impurity species.
    d_par_excess = D_par_local + D_par_sc_num * tau_sc - D_prof
    d_par_excess_imp = D_par_local_imp + D_par_imp_sc_num * tau_sc - D_prof_imp
    psi_gradient = grad(psi)

    B = (
        v * R * (
            particle_source + source_pellet + source_bg_drift + source_imp_drift
        ) * xjac
        + _species_advection(v, rho)
        + _parallel_convection(v, rho)
        + _species_tgnum(v, rho, tgnum_rho)
        # The main species diffuses with the excess of the bulk diffusivity and
        # the impurity one with its own.
        + _parallel_diffusion(v, rho - rhoimp, d_par_excess)
        + _parallel_diffusion(v, rhoimp, d_par_excess_imp)
        - D_prof * R * _perpendicular_diffusion(v, rho - rhoimp) * xjac
        - D_prof_imp * R * _perpendicular_diffusion(v, rhoimp) * xjac
        - D_perp_num_psin * _laplacian(v) * _laplacian(rho) * R * xjac
        # Diamagnetic drift, atomic sources and the kinetic coupling.
        + v * 2 * tauIC * 2 * dZ(_diamagnetic_pressure(with_TiTe=with_TiTe))
        * R * xjac
        + v * electron_density * rhon * R * Sion_rate(thermal) * xjac
        - v * electron_density * (rho - rhoimp) * R * Srec_rate(thermal) * xjac
        + v * R * aux_rho0 * xjac
        # Weak form of -div(rho*V_pinch) with
        # V_pinch = -V_prof_pinch*grad(psi)/|grad(psi)|.  The pinch direction
        # is built from the evolved flux, so it is varied like any other state
        # quantity; the element routine freezes it (see JOREK_FINDINGS.md).
        - V_prof_pinch / sp.sqrt(psi_gradient[0]**2 + psi_gradient[1]**2)
        * dot(grad(v), psi_gradient) * rho * R * xjac
    )
    A = v * rho * R * xjac
    return EvolutionEquation("model600_density", v, A, B)

def parallel_velocity_equation_vpar(*, with_TiTe=False, element_brackets=True):
    """Model-600 parallel velocity equation (``var_vpar``).

    The perpendicular parallel viscosity follows the default
    ``normalized_velocity_profile = .true.`` branch of the element routine.

    ``element_brackets`` selects how the poloidal bracket of the parallel
    kinetic-energy flux is spelled.  The element routine writes this one group
    in element coordinates in the residual and in the ``rho`` and ``vpar``
    columns, and in physical coordinates in ``amat(var_vpar,var_psi)``.  The
    pressure gradient and the Taylor-Galerkin terms use element coordinates
    throughout and need no such switch.
    """

    v = test_function("v")
    thermal = Te if with_TiTe else T
    bb2 = _B2(psi)
    pressure = rho * (Ti + Te) if with_TiTe else rho * T
    if with_TiTe:
        pressure += rhoimp * (alpha_i_state() * Ti + alpha_e_temperature(Te))
    else:
        pressure += rhoimp * alpha_imp_temperature(T)
    electron_density = rho + alpha_e_state(thermal) * rhoimp
    kinetic_gradient = (
        _b_dot_grad_element if element_brackets else _b_dot_grad
    )
    # Only the flux is varied in the prescribed rotation profile, so
    # ``grad(Vt) = dV_dpsi_source*grad(psi)`` is the faithful form.
    rotation_shear = tuple(
        parallel - dV_dpsi_source * flux
        for parallel, flux in zip(grad(vpar), grad(psi))
    )
    visco_par_eff = visco_par + visco_par_sc_num * tau_sc
    rho_hat = _rho_hat()
    psi_gradient = grad(psi)

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
        # Momentum carried by the particle sources; the conservative form
        # moves part of it into A.
        - v * (
            particle_source + source_pellet + source_bg_drift
            + source_imp_drift
        ) * vpar * bb2 * R * xjac * (1 - fact_conservative_u)
        # -(d_t rho + div(rho v)) v_par B**2 R; the d_t rho part belongs to A.
        + fact_conservative_u * v * vpar * bb2 * xjac * (
            _poloidal_cross(rho_hat, u) - R * _b_dot_grad(rho * vpar)
        )
        + (1 - delta_n_convection) * (
            -v * electron_density * rhon * Sion_rate(thermal)
            * vpar * bb2 * R * xjac
            + v * electron_density * (rho - rhoimp) * Srec_rate(thermal)
            * vpar * bb2 * R * xjac
        )
        # Taylor-Galerkin stabilization.
        - tgnum_vpar * sp.Rational(1, 4) * rho * vpar**2 * bb2 * R
        * _b_dot_grad_element(vpar) * _b_dot_grad_element(v) * xjac * tstep
        - tgnum_vpar * sp.Rational(1, 4) * v * vpar**2 * bb2 * R
        * (1 - fact_conservative_u)
        * _b_dot_grad_element(vpar) * _b_dot_grad_element(rho) * xjac * tstep
        - tgnum_vpar * sp.Rational(1, 4) * vpar**3 * bb2 * R
        * fact_conservative_u
        * _b_dot_grad_element(rho) * _b_dot_grad_element(v) * xjac * tstep
        # Kinetic coupling and the inward pinch.
        - v * aux_rho0 * vpar * bb2 * R * (1 - fact_conservative_u) * xjac
        + v * R * aux_mom_par0 * xjac
        + V_prof_pinch / sp.sqrt(psi_gradient[0]**2 + psi_gradient[1]**2)
        * dot(psi_gradient, grad(vpar)) * rho * v * R * xjac
    )
    # The parallel momentum density is rho*v_par*B**2*R.  The element routine
    # freezes the density correction in the tangent, as in the perpendicular
    # momentum equation.
    A = (
        v * corr_neg_dens(freeze(rho)) * vpar * bb2 * R * xjac
        + fact_conservative_u * v * rho * freeze(vpar) * bb2 * R * xjac
    )
    return EvolutionEquation("model600_parallel_velocity", v, A, B)

def impurity_density_equation_rhoimp(*, with_TiTe=False):
    """Model-600 impurity-density equation (``var_rhoimp``).

    The impurity species is advected and diffused like the main one, but it
    has no atomic sources and no diamagnetic contribution; the element routine
    marks the latter with an explicit placeholder comment.
    """

    v = test_function("v")
    d_par_excess_imp = D_par_local_imp + D_par_imp_sc_num * tau_sc - D_prof_imp

    B = (
        _species_advection(v, rhoimp)
        + _parallel_convection(v, rhoimp)
        + _species_tgnum(v, rhoimp, tgnum_rhoimp)
        + _parallel_diffusion(v, rhoimp, d_par_excess_imp)
        - D_prof_imp * R * _perpendicular_diffusion(v, rhoimp) * xjac
        - Dn_perp_num * _laplacian(v) * _laplacian(rhoimp) * R * xjac
        + R * v * source_imp_drift * xjac
    )
    A = v * rhoimp * R * xjac
    return EvolutionEquation("model600_impurity_density", v, A, B)

def ion_energy_equation_Ti(*, element_brackets=True):
    """Model-600 ion energy equation (``var_Ti``).

    ``element_brackets`` selects the spelling of the poloidal advection
    bracket: the element routine writes it in element coordinates in the
    residual and in the ``rho``, ``Ti`` and ``rhoimp`` columns, and in
    physical coordinates in ``amat(var_Ti,var_u)``.
    """

    v = test_function("v")
    bracket = element_bracket if element_brackets else _physical_bracket
    ion_density = rho + alpha_i_state() * rhoimp
    ion_pressure = ion_density * Ti

    B = (
        v * R * heat_source_i * xjac
        + _pressure_transport(v, ion_pressure, bracket=bracket)
        + _heat_conduction(v, Ti, ZKi_par(Ti), ZKi_perp(rho), ZK_i_perp_num_psin)
        + _energy_tgnum(v, ion_pressure, tgnum_Ti)
        + _friction_heating(v, _released_kinetic_energy(Te))
        + _viscous_heating(v, visco_heating(Te))
        + _heating_floor(v, Ti_floor_exp(Ti), Ti_floor(Ti), Tie_min_neg)
        + _kinetic_coupling(v)
        # Ion-electron energy exchange, recombination sink, kinetic coupling.
        + v * R * Ti_e_exchange(Ti, Te, rho, rhoimp) * xjac
        - v * Ti * R * corr_neg_dens(rho)**2 * Srec_rate(Te) * xjac
        + v * R * aux_E0_Ti * xjac
    )
    A = v * (
        corr_neg_dens(rho) + alpha_i_state() * corr_neg_dens_imp(rhoimp)
    ) * Ti * R * xjac
    return EvolutionEquation("model600_ion_energy", v, A, B)

def electron_energy_equation_Te(*, element_brackets=True):
    """Model-600 electron energy equation (``var_Te``).

    The electron pressure follows ``construct_pressure``:
    ``Pe = rho*Te + rhoimp*alpha_e(Te)*Te``, so a gradient of it carries
    ``alpha_e`` on the density part and ``alpha_e_bis`` on the temperature
    part, exactly as the element routine spells it out term by term.

    ``element_brackets`` applies to the electron-pressure advection only; the
    ionization-energy advection is always written in element coordinates,
    including in ``amat(var_Te,var_u)`` where the pressure bracket right above
    it uses physical ones.
    """

    v = test_function("v")
    bracket = element_bracket if element_brackets else _physical_bracket
    electron_pressure = rho * Te + rhoimp * alpha_e_temperature(Te)
    electron_density = corr_neg_dens(rho) + alpha_e_state(Te) * corr_neg_dens_imp(rhoimp)
    ionization_energy = _ionization_energy(Te)

    B = (
        v * R * heat_source_e * xjac
        + _pressure_transport(v, electron_pressure, bracket=bracket)
        + _heat_conduction(v, Te, ZKe_par(Te), ZKe_perp(rho), ZK_e_perp_num_psin)
        + _energy_tgnum(v, electron_pressure, tgnum_Te)
        + _ohmic_heating(v, eta_ohm_e(Te, rho, rhoimp))
        + _radiation_sinks(v, electron_density, Te)
        + _heating_floor(v, Te_floor_exp(Te), Te_floor(Te), Tie_min_neg)
        + _ionization_energy_transport(v, ionization_energy, Te)
        # Ionization sink, ion-electron exchange, plasmoid-drift teleportation.
        - v * R * ksi_ion_norm * (rho + alpha_e_state(Te) * rhoimp)
        * rhon * Sion_rate(Te) * xjac
        + v * R * Te_i_exchange(Ti, Te, rho, rhoimp) * xjac
        + v * R * power_dens_teleport_ju * xjac
        + v * R * aux_E0_Te * xjac
    )
    A = v * R * xjac * (
        corr_neg_dens(rho) * Te + corr_neg_dens_imp(rhoimp) * alpha_e_temperature(Te)
        + (GAMMA - 1) * ionization_energy
    )
    return EvolutionEquation("model600_electron_energy", v, A, B)

def total_energy_equation_T(*, element_brackets=True):
    """Model-600 single-temperature energy equation (``var_T``).

    It carries the ion and the electron contributions at once: advection and
    conduction of the total pressure, the friction and viscous-heating terms
    of the ion equation, and the Ohmic, radiation and ionization-energy terms
    of the electron equation.  The total pressure follows the same convention
    as the two-temperature ones, ``P = rho*T + rhoimp*alpha_imp(T)*T``.
    """

    v = test_function("v")
    bracket = element_bracket if element_brackets else _physical_bracket
    pressure = rho * T + rhoimp * alpha_imp_temperature(T)
    electron_density = corr_neg_dens(rho) + alpha_e_state(T) * corr_neg_dens_imp(rhoimp)
    ionization_energy = _ionization_energy(T)

    B = (
        v * R * heat_source_total * xjac
        + _pressure_transport(v, pressure, bracket=bracket)
        + _heat_conduction(v, T, ZK_par(T), ZK_perp(rho), ZK_perp_num_psin)
        + _energy_tgnum(v, pressure, tgnum_T)
        + _ohmic_heating(v, eta_ohm_e(T, rho, rhoimp))
        + _friction_heating(v, _released_kinetic_energy(T))
        + _viscous_heating(v, visco_heating(T))
        + _radiation_sinks(v, electron_density, T)
        + _heating_floor(v, T_floor_exp(T), T_floor(T), T_min_neg)
        + _ionization_energy_transport(v, ionization_energy, T)
        + _kinetic_coupling(v)
        # Ionization and recombination sinks, teleportation, kinetic coupling.
        - v * R * ksi_ion_norm * (rho + alpha_e_state(T) * rhoimp)
        * rhon * Sion_rate(T) * xjac
        - (GAMMA - 1) * v * sp.Rational(1, 2) * T * R
        * corr_neg_dens(rho)**2 * Srec_rate(T) * xjac
        + v * R * power_dens_teleport_ju * xjac
        + v * R * aux_E0 * xjac
    )
    A = v * R * xjac * (
        corr_neg_dens(rho) * T + corr_neg_dens_imp(rhoimp) * alpha_imp_temperature(T)
        + (GAMMA - 1) * ionization_energy
    )
    return EvolutionEquation("model600_total_energy", v, A, B)
