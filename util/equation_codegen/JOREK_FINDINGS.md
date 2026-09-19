# Model-600 Jacobian: missing and inconsistent terms

This file records the differences between the volume terms assembled by
`models/model600/mod_elt_matrix_fft.f90` and the automatically linearized
weak forms in `src/jorek_equations/model600.py`.

Every finding below was obtained by generating the two term reports and
comparing them block by block:

```bash
.venv/bin/python examples/export_model600_terms.py
.venv/bin/python examples/diff_model600_reports.py
```

`diff_model600_reports.py` compares each `rhs_ij`/`amat` assignment as a
multiset of monomials and prints the algebraic residual, so a block that only
reorders its terms is reported as identical. At the time of writing, 81 of the
100 exported blocks are identical; the remaining 19 are the five findings
below. Nothing in this list has been fixed in the Fortran.

Status of the five exported rows:

| Row | Residual | Jacobian |
|---|---|---|
| `psi` | reproduced | finding 2 |
| `u` | reproduced | findings 1 and 3 |
| `zj` | reproduced | reproduced |
| `w` | reproduced | reproduced |
| `rho` | reproduced | findings 1, 4 and 5 |

Every source term of every row is reproduced by the generator. With the single
exception of finding 3, the differences are all terms the generator produces
and the element routine does not.

The NEO branch is excluded from the default reports (`--include-neo` enables
it) and has not been audited here.

---

## Background: the one-temperature factor of one half

`construct_pressure` is written once for both temperature models and always
builds the species pressures from the species temperatures:

```fortran
Pi0 = (r0 + rimp0*alpha_i) * Ti0
Pe0 = (r0 + rimp0*alpha_e) * Te0
P0  = Pi0 + Pe0
```

The one-temperature model evolves the **total** temperature and enters that
routine with `Ti0 = Te0 = T0/2` (`mod_elt_matrix_fft.f90`, the `Ti0 = T0/2.d0`
block). Two consequences matter for the Jacobian:

* a tangent taken with respect to `T` must carry the chain factor
  `dTi0/dT0 = dTe0/dT0 = 1/2`;
* the stored one-temperature closure values are the means of the per-species
  ones, `alpha_imp = (alpha_i+alpha_e)/2`,
  `alpha_imp_bis = (alpha_i+alpha_e_bis)/2` and `alpha_imp_tri = alpha_e_tri/4`.

Most of the one-temperature Jacobian applies this correctly — it is why, for
example, `amat(var_u,var_T)` uses `tauIC` where the residual uses `tauIC*2.`,
and why `amat_n(var_psi,var_T)` uses `tauIC` where `amat_n(var_psi,var_Te)`
uses `tauIC*2.`. Finding 3 is the one place where it does not.

---

## Finding 1 — the impurity ion pressure is not differentiated in the diamagnetic terms

**Where:** `amat(var_u,var_Ti)` (two-temperature branch),
`amat(var_u,var_T)` (one-temperature branch) and `amat(var_u,var_rhoimp)` in
the `var_u` block of `mod_elt_matrix_fft.f90`, and the same three columns of
the `var_rho` block.

**What happens:** the momentum residual uses the full ion pressure,

```fortran
- v * tauIC*2. * BigR**4 * (Pi0_s * w0_t - Pi0_t * w0_s)                                  * tstep
- tauIC*2. * BigR**3 * Pi0_y * (v_x*u0_x + v_y*u0_y)                               * xjac * tstep
- v * tauIC*2. * BigR**4 * (u0_xy*(Pi0_xx - Pi0_yy) - Pi0_xy*(u0_xx - u0_yy))      * xjac * tstep
```

with `Pi0 = (r0 + rimp0*alpha_i)*Ti0`, but the corresponding tangents are
written as if `Pi0 = r0*Ti0`. For instance `amat(var_u,var_Ti)` contains

```fortran
+ v * tauIC*2. * BigR**4 * r0 * (Ti_s * w0_t - Ti_t * w0_s)                * theta * tstep
+ v * tauIC*2. * BigR**4 * Ti * (r0_s * w0_t - r0_t * w0_s)                * theta * tstep
+ tauIC*2. * BigR**3 * (r0_y*Ti + r0*Ti_y) * (v_x*u0_x + v_y*u0_y) * xjac  * theta * tstep
+ v * tauIC*2. * BigR**4 * ( (u0_xy * (Ti_xx*r0 + 2.d0*Ti_x*r0_x + Ti*r0_xx    &
                                     - Ti_yy*r0 - 2.d0*Ti_y*r0_y - Ti*r0_yy))  &
                           - (Ti_xy*r0 + Ti_x*r0_y + Ti_y*r0_x + Ti*r0_xy)     &
                             * (u0_xx - u0_yy) )                          * xjac * theta * tstep
```

and `amat(var_u,var_rhoimp)` contains no diamagnetic contribution at all.

**Suggested fix:** in the four lines above replace every background density by
the ion pressure density, i.e. `r0 -> (r0 + rimp0*alpha_i)`,
`r0_s -> (r0_s + rimp0_s*alpha_i)`, `r0_x -> (r0_x + rimp0_x*alpha_i)` and so
on (`alpha_i` carries no temperature dependence, `dalpha_i_dT = 0`, so no
`bis`/`tri` coefficients are needed). The same four lines with
`r0 -> alpha_i*rhoimp`, `Ti -> Ti0` are the contribution missing from
`amat(var_u,var_rhoimp)`. In the one-temperature branch, `amat(var_u,var_T)`
needs the same substitution with the chain factor of `1/2` that the block
already applies through its `tauIC` prefactor.

The density equation carries the same pressure through its diamagnetic drift
term,

```fortran
+ v * 2.d0 * tauIC*2. * Pi0_y * BigR * xjac * tstep * factor(var_rho,7)
```

whose tangents are again written for `Pi0 = r0*Ti0` only:

```fortran
amat(var_rho,var_Ti) = - v * 2.d0 * tauIC*2. * (Ti_y * r0 + Ti*r0_y) * BigR * xjac * theta * tstep
amat(var_rho,var_rho) = ... - v * 2.d0 * tauIC*2. * (rho_y * Ti0 + rho*Ti0_y) * BigR * xjac * theta * tstep
```

`amat(var_rho,var_rho)` is correct, because `Pi0` depends on `rho` only
through `r0*Ti0`; `amat(var_rho,var_Ti)` (and its one-temperature
`amat(var_rho,var_T)` counterpart) needs `r0 -> (r0 + rimp0*alpha_i)`, and
`amat(var_rho,var_rhoimp)` needs the new contribution

```fortran
- v * 2.d0 * tauIC*2. * alpha_i * (rhoimp_y * Ti0 + rhoimp * Ti0_y) * BigR * xjac * theta * tstep
```

**Scale:** 22 monomials missing from each of `amat(var_u,var_Ti)`,
`amat(var_u,var_T)` and `amat(var_u,var_rhoimp)`, and 2 monomials from each of
`amat(var_rho,var_Ti)`, `amat(var_rho,var_T)` and `amat(var_rho,var_rhoimp)`
(each column in both temperature branches), 96 monomials in total.

**Related, not visible in the reports:** the work variables that feed the
diamagnetic viscosity drop the same contribution and would need the same
substitution:

```fortran
Pi0_x_Ti   = r0_x  * Ti +       r0   * Ti_x
Pi0_xx_Ti  = r0_xx * Ti + 2.0 * r0_x * Ti_x  + r0 * Ti_xx
...
Pi0_x_rho  = rho_x  * Ti0 +       rho   * Ti0_x
Pi0_xx_rho = rho_xx * Ti0 + 2.0 * rho_x * Ti0_x + rho * Ti0_xx
```

`W_dia_Ti` and `W_dia_rho` are built from these, and `W_dia` itself is built
from the full `Pi0`, so the tangents are inconsistent with the value. The
generator cannot see this because it treats `W_dia` as an opaque external
quantity; the inconsistency was found by reading `construct_pressure` against
the `Pi0_*_Ti` definitions.

---

## Finding 2 — the impurity electron pressure is not differentiated in the induction equation

**Where:** `amat(var_psi,var_Te)`, `amat_n(var_psi,var_Te)` and
`amat(var_psi,var_rhoimp)` (and their one-temperature `var_T` counterparts).
`amat_n(var_psi,var_rhoimp)` is missing from the element routine entirely.

**What happens:** the induction residual uses `Pe0 = (r0 + rimp0*alpha_e)*Te0`
in its diamagnetic coupling, but the tangents are written as if
`Pe0 = r0*Te0`:

```fortran
amat(var_psi,var_Te) = ... &
   + v * tauIC*2./(r0_corr*BB2) * F0**2/BigR**2 * r0 * (ps0_s * Te_t  - ps0_t * Te_s) * theta * tstep &
   + v * tauIC*2./(r0_corr*BB2) * F0**2/BigR**2 * Te * (ps0_s * r0_t  - ps0_t * r0_s) * theta * tstep &
   - v * tauIC*2./(r0_corr*BB2) * F0**3/BigR**3 * Te * r0_p                    * xjac * theta * tstep

amat_n(var_psi,var_Te) = - v * tauIC*2./(r0_corr*BB2) * F0**3/BigR**3 * r0 * Te_p * xjac * theta * tstep
```

and `amat(var_psi,var_rhoimp)` carries only the resistivity tangent
`- deta_drimp0 * v * rhoimp * (zj0-current_source-Jb-aux_jre_ind)/BigR`.

**Suggested fix:** following `construct_pressure`, the electron-temperature
tangent needs `r0 -> (r0 + rimp0*alpha_e_bis)` wherever it multiplies a
derivative of the trial temperature, `r0_t -> (r0_t + rimp0_t*alpha_e)`
wherever it multiplies the trial temperature itself, plus the second-derivative
term `rimp0*alpha_e_tri*Te*Te0_t` that comes from differentiating
`alpha_e_bis`. The `rhoimp` column needs the whole `Pe0` impurity tangent,

```fortran
+ v * tauIC*2./(r0_corr*BB2) * F0**2/BigR**2                                    &
    * ( ps0_s * (rhoimp_t*alpha_e*Te0 + rhoimp*alpha_e_bis*Te0_t)               &
      - ps0_t * (rhoimp_s*alpha_e*Te0 + rhoimp*alpha_e_bis*Te0_s) ) * theta * tstep &
- v * tauIC*2./(r0_corr*BB2) * F0**3/BigR**3                                    &
    * (rhoimp_p*alpha_e*Te0 + rhoimp*alpha_e_bis*Te0_p)            * xjac * theta * tstep
```

plus a new `amat_n(var_psi,var_rhoimp)` assignment for the toroidal channel,

```fortran
amat_n(var_psi,var_rhoimp) = &
    + v * tauIC*2./(r0_corr*BB2) * F0**3/BigR**3 * alpha_e * Te0 * rhoimp_p * xjac * theta * tstep
```

**Scale:** 8 monomials missing from `amat(var_psi,var_Te)`, 1 from
`amat_n(var_psi,var_Te)`, 5 from `amat(var_psi,var_rhoimp)` and 1 from the
absent `amat_n(var_psi,var_rhoimp)`, in each temperature branch.

---

## Finding 3 — the one-temperature diamagnetic-viscosity tangent is a factor two too large

**Where:** the `else` (single-temperature) branch of the `var_u` block,
`amat(var_u,var_T)`:

```fortran
! --- Contributions of the diamagnetic viscosity
- dvisco_dT     * bigR * W_dia_Ti * (v_x*Ti0_x + v_y*Ti0_y)  * xjac * theta * tstep  &
- visco_T       * bigR * W_dia_Ti * (v_xx + v_x/bigR + v_yy) * xjac * theta * tstep  &
- dvisco_dT     * bigR * W_dia    * (v_x*T_x  + v_y*T_y )    * xjac * theta * tstep  &
```

**What happens:** these three lines are the tangent of the residual terms

```fortran
+ dvisco_dT * bigR * W_dia * (v_x*Ti0_x + v_y*Ti0_y) * xjac * tstep
+ visco_T   * bigR * W_dia * (v_xx + v_x/bigR + v_yy) * xjac * tstep
```

with respect to `T`. Both factors that depend on the ion temperature are
differentiated without the chain factor `dTi0/dT0 = 1/2`:

* the third line uses `(v_x*T_x + v_y*T_y)`; the variation of `Ti0_x` is
  `Ti_x = T_x/2`, so it should read `(v_x*T_x + v_y*T_y)/2`;
* `W_dia_Ti` is assembled from `Pi0_x_Ti`, `Pi0_xx_Ti` and `Pi0_yy_Ti`, which
  are built from the **trial basis function** `Ti`, `Ti_x`, `Ti_xx`. In the
  single-temperature branch that basis function is the `var_T` trial function,
  so `W_dia_Ti = 2 * dW_dia/dT0` and the first two lines are twice their
  correct value as well.

The remaining two lines of the same block,

```fortran
- d2visco_dT2*T * bigR * W_dia * (v_x*Ti0_x + v_y*Ti0_y)  * xjac * theta * tstep  &
- dvisco_dT*T   * bigR * W_dia * (v_xx + v_x/bigR + v_yy) * xjac * theta * tstep  &
```

are correct: they differentiate `visco(T0)`/`dvisco_dT(T0)`, which depend on
the evolved temperature directly and carry no chain factor.

**Suggested fix:** halve the three lines listed first, for example by writing
them with `W_dia_Ti/2.d0` and `(v_x*T_x + v_y*T_y)/2.d0`.

**Effect:** the two-temperature `amat(var_u,var_Ti)` block is unaffected — the
chain factor is one there — so this is a one-temperature-only error. It makes
the corresponding Newton block inconsistent with the residual and can slow or
prevent convergence when the diamagnetic viscosity is active
(`Wdia = .true.`) in a single-temperature run.

---

## Finding 4 — the density inward pinch is not differentiated with respect to psi

**Where:** the `var_rho` block, `amat(var_rho,var_psi)`.

**What happens:** the density residual contains the weak form of
`-div(r0 * V_pinch)` with `V_pinch = -V_prof_pinch * grad(psi)/|grad(psi)|`:

```fortran
- V_prof_pinch / sqrt(psi_grad2) * (v_x * ps0_x + v_y * ps0_y) * r0 * BigR * xjac * tstep * factor(var_rho,14)
```

with `psi_grad2 = ps0_x**2 + ps0_y**2`. The pinch direction is therefore a
function of the evolved poloidal flux, but `amat(var_rho,var_psi)` contains no
corresponding tangent. Every other flux-dependent quantity of the same
equation *is* differentiated there — `BB2` through `BB2_psi`, the parallel
gradients through `Bgrad_rho_star_psi`, `Bgrad_rho_psi` and
`Bgrad_rhoimp_psi` — which is what makes the omission stand out.

**Suggested fix:** the missing contribution is the transverse projection of
the trial flux gradient,

```fortran
+ V_prof_pinch / sqrt(psi_grad2) * r0 * BigR * xjac * theta * tstep                       &
    * ( (v_x * psi_x + v_y * psi_y)                                                       &
      - (ps0_x*v_x + ps0_y*v_y) * (ps0_x*psi_x + ps0_y*psi_y) / psi_grad2 )
```

(the first bracket from varying `grad(psi)` in the dot product, the second from
varying `1/|grad(psi)|`).

**Scale:** 6 monomials, in each temperature branch.

**Note:** `V_prof_pinch`, `D_prof`, `D_prof_imp`, `D_par_local`,
`D_par_local_imp`, `tau_sc` and `D_perp_num_psin` are profile work values that
also depend on the state through `psi_norm` (and, for `tau_sc`, through the
pressure). The element routine freezes all of them, and the generator follows
that convention, so they produce no report differences. The pinch term is
different: its flux dependence is explicit in the assembled expression rather
than hidden inside a work value.

---

## Finding 5 — the temperature dependence of alpha_e is not differentiated in the density source terms

**Where:** the `var_rho` block, `amat(var_rho,var_Te)` (two-temperature) and
`amat(var_rho,var_T)` (one-temperature).

**What happens:** the ionization and recombination sources of the density
equation carry the electron-fraction coefficient `alpha_e`, which is a
function of the electron temperature (`alpha_e = m_i_over_m_imp*Z_imp - 1`,
with `dalpha_e_dT = m_i_over_m_imp*dZ_imp_dT`):

```fortran
+ v * (r0+alpha_e*rimp0) * rn0 * BigR * Sion_T          * xjac * tstep * factor(var_rho,8) &
- v * (r0+alpha_e*rimp0) * (r0-rimp0) * BigR * Srec_T   * xjac * tstep * factor(var_rho,9)
```

but the temperature tangent differentiates only the rates:

```fortran
amat(var_rho,var_Te) = - v * BigR * (r0+alpha_e*rimp0) * rn0 * dSion_dT * Te        * xjac * theta * tstep &
                       + v * BigR * (r0+alpha_e*rimp0) * (r0-rimp0) * dSrec_dT * Te * xjac * theta * tstep
```

The momentum equation differentiates the same coefficient in the same source
terms — `amat(var_u,var_Te)` contains
`- BigR**3 * (dalpha_e_dT * rimp0 * rn0 * Sion_T * Te) * ...` and its
recombination partner — so the two equations disagree about whether `alpha_e`
is a state function.

**Suggested fix:** add, to `amat(var_rho,var_Te)` (and to
`amat(var_rho,var_T)` with the same spelling in `T`):

```fortran
- v * BigR * dalpha_e_dT * rimp0 * rn0 * Sion_T * Te          * xjac * theta * tstep &
+ v * BigR * dalpha_e_dT * rimp0 * (r0-rimp0) * Srec_T * Te   * xjac * theta * tstep
```

**Scale:** 3 monomials, in each temperature branch.

---

## Reproducing

```bash
cd util/equation_codegen
python3 -m venv .venv && .venv/bin/python -m pip install -e .
.venv/bin/python examples/export_model600_terms.py
.venv/bin/python examples/diff_model600_reports.py
```

Every line that `diff_model600_reports.py` prints as `generated` and not as
`source` is a term the linearization produces and the element routine does
not. A block reported with both `source` and `generated` lines and a non-zero
residual is a coefficient disagreement; `amat(var_u,var_T)` is currently the
only one.
