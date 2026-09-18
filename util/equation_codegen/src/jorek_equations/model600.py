"""Model-600 weak equations, extending the model-199 DSL definitions."""

from .equations import ConstraintEquation, EvolutionEquation
from .external import external_function
from .model199 import FIELDS as MODEL199_FIELDS
from .model199 import j, omega, psi, rho, T, u, xjac
from .operators import R, dphi, element_bracket, dot, grad
from .symbols import coefficient, field, test_function


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


vpar_psi_equation_1 = induction_equation_1


def _parallel_norm(flux):
    return (F0**2 + grad(flux)[0]**2 + grad(flux)[1]**2) / R**2
