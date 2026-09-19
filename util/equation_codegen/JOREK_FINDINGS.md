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
reorders its terms is reported as identical.  The residual is evaluated after
rewriting every `_s`/`_t` derivative through the chain rule, because the
element routine spells the same poloidal bracket sometimes as
`a_s*b_t - a_t*b_s` and sometimes as `xjac*(a_x*b_y - a_y*b_x)` — even between
a residual and its own tangent.  A block whose two sides differ only by that
rewrite is reported as `SAME`. At the time of writing, 127 of the
162 exported blocks are identical; the remaining 35 are the ten findings
below. Nothing in this list has been fixed in the Fortran.

Status of the five exported rows:

| Row | Residual | Jacobian |
|---|---|---|
| `psi` | reproduced | finding 2 |
| `u` | reproduced | findings 1 and 3 |
| `zj` | reproduced | reproduced |
| `w` | reproduced | reproduced |
| `rho` | reproduced | findings 1, 4 and 5 |
| `vpar` | findings 6 and 9 | findings 4, 6, 7, 8 and 9 |
| `rhoimp` | finding 10 | reproduced |

Every source term of every row is reproduced by the generator. Apart from
findings 3, 6 and 9, the differences are all terms the generator produces and
the element routine does not.

Findings 3, 6, 9 and 10 are disagreements about a coefficient or a sign rather
than omissions, and finding 10 is the only one that sits in a residual rather
than in a tangent.

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

The parallel-velocity equation advects `v_par` with the same pinch velocity
and has the same gap: `amat(var_vpar,var_psi)` carries no pinch tangent
either.

**Scale:** 6 monomials in `amat(var_rho,var_psi)` and 6 in
`amat(var_vpar,var_psi)`, in each temperature branch.

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

## Finding 6 — the parallel-velocity time term is linearized inconsistently

**Where:** the `var_vpar` block: `rhs_ij(var_vpar)` (the `zeta` history),
`amat(var_vpar,var_vpar)`, `amat(var_vpar,var_psi)` and
`amat(var_vpar,var_rho)`.

**What happens:** the parallel momentum density is `r0_corr*vpar*BB2*BigR`,
with `BB2 = (F0**2 + ps0_x**2 + ps0_y**2)/BigR**2` — the same `BB2` the whole
equation uses for its inertia, source and Taylor-Galerkin terms. Its
linearization is implemented as three pieces that do not come from one
functional:

```fortran
amat(var_vpar,var_vpar) = v * Vpar * r0_corr * F0**2 / BigR * xjac * (1.d0 + zeta) ...
amat(var_vpar,var_psi)  = ... + v * r0_corr * vpar0 / BigR * (ps0_x*psi_x + ps0_y*psi_y) * xjac * (1.d0 + zeta) ...
amat(var_vpar,var_rho)  = ... + fact_conservative_u * v * rho * vpar0 * F0**2 / BigR * xjac * (1.d0 + zeta) ...
rhs_ij(var_vpar)        = ... + zeta * v * delta_g(mp,var_vpar,ms,mt) * r0_corr * F0**2 / BigR * xjac    &
                              + zeta * v * r0_corr * vpar0 * (ps0_x*delta_ps_x + ps0_y*delta_ps_y) / BigR * xjac ...
```

* the `vpar` column and the `vpar` history keep only the toroidal part
  `F0**2/BigR` of `BB2*BigR`; the poloidal part
  `(ps0_x**2 + ps0_y**2)/BigR` is missing;
* the `psi` column carries `(ps0_x*psi_x + ps0_y*psi_y)`, which is exactly
  **half** of the variation of `BB2*BigR`, namely
  `2*(ps0_x*psi_x + ps0_y*psi_y)/BigR`;
* the conservative `fact_conservative_u` copy of the same term repeats the
  first problem in its `rho` column and has no `psi` column at all.

There is no scalar `A` whose Gateaux derivative gives these terms together:
the mixed second derivative is not symmetric. Differentiating the implemented
`vpar` column with respect to the flux gives zero, because `F0**2/BigR`
carries no flux dependence, while differentiating the implemented `psi`
column with respect to `vpar` gives
`v*r0_corr*(ps0_x*psi_x + ps0_y*psi_y)/BigR`.

Because the `zeta` history in `rhs_ij(var_vpar)` uses the same combination as
the `(1.d0 + zeta)` mass terms, the time scheme is self-consistent; what it
discretises, however, is not `d/dt (rho*vpar*BB2*BigR)`, which is the density
the rest of the equation assumes (`BB2` appears in its inertia, source,
kinetic-coupling and Taylor-Galerkin terms).

**Not specific to model 600.** The same three expressions appear verbatim in
the `var_vpar` blocks of models 303, 305, 306, 307, 401, 501 and 502, so this
is a long-standing shared convention rather than a model-600 slip, and a fix
would have to be applied consistently across them.

**Suggested fix:** derive all four from `A = v*r0_corr*vpar*BB2*BigR*xjac`,
that is use `BB2*BigR` in place of `F0**2/BigR` in the `vpar` and `rho`
columns and in the history, and double the `psi` column.

**Scale:** 8 residual monomials in `rhs_ij(var_vpar)`, 8 in
`amat(var_vpar,var_psi)`, 4 in `amat(var_vpar,var_vpar)` and 4 in
`amat(var_vpar,var_rho)`, in each temperature branch — half of them from the
main mass term and half from its `fact_conservative_u` copy.

---

## Finding 7 — BB2 is not differentiated in the tgnum_vpar tangent

**Where:** `amat(var_vpar,var_psi)` and `amat_k(var_vpar,var_psi)`.

**What happens:** all three Taylor-Galerkin terms of the parallel-velocity
residual carry the factor `BB2`, for example

```fortran
- tgnum_vpar * 0.25d0 * r0 * Vpar0**2 * BB2                                             &
          * (-(ps0_s*vpar0_t - ps0_t*vpar0_s)/xjac + F0/BigR*vpar0_p) / BigR            &
          * (-(ps0_s*v_t     - ps0_t*v_s)    /xjac)  * xjac * tstep * tstep
```

`BB2` depends on the poloidal flux, and the element routine differentiates it
everywhere else in this equation (the `BB2_psi` terms of the inertia, source
and kinetic-coupling contributions). The `tgnum_vpar` flux tangent, however,
varies only the two poloidal brackets: every one of its six terms still
carries `BB2`, never `BB2_psi`.

**Suggested fix:** add, for each of the three residual terms, the
corresponding `BB2 -> BB2_psi` copy to `amat(var_vpar,var_psi)` and
`amat_k(var_vpar,var_psi)`.

**Scale:** 84 residual monomials in `amat(var_vpar,var_psi)` and 20 in
`amat_k(var_vpar,var_psi)`, in each temperature branch.

---

## Finding 8 — the toroidal channel of the parallel-parallel viscosity tangent is missing

**Where:** `amat_n(var_vpar,var_vpar)` and `amat_kn(var_vpar,var_vpar)`.

**What happens:** the residual term

```fortran
- visco_par_par * F0**2 / (BigR * BB2) * Bgrad_vpar * Bgrad_rho_star * xjac * tstep
```

uses the full parallel gradient
`Bgrad_vpar = (F0/BigR*vpar0_p + vpar0_x*ps0_y - vpar0_y*ps0_x)/BigR`, but its
trial-function counterpart is defined with the poloidal half only:

```fortran
Bgrad_vpar_vpar = ( vpar_x * ps0_y - vpar_y * ps0_x ) / BigR
```

so the `F0/BigR*vpar_p` part never reaches the Jacobian. Consequently
`amat_n(var_vpar,var_vpar)` and `amat_kn(var_vpar,var_vpar)` carry no
`visco_par_par` contribution at all, even though the same tangent is present
in the `p` and `k` channels.

**Suggested fix:** add a toroidal companion
`Bgrad_vpar_vpar_n = (F0/BigR*vpar_p)/BigR` and use it in the `n` and `kn`
channels, exactly as `Bgrad_rho_rho_n` is used in the density equation.

**Scale:** 2 monomials in `amat_n(var_vpar,var_vpar)` and 1 in
`amat_kn(var_vpar,var_vpar)`, in each temperature branch.

---

## Finding 9 — the parallel-velocity inward-pinch tangent has the wrong sign

**Where:** `amat(var_vpar,var_rho)` and `amat(var_vpar,var_vpar)`.

**What happens:** the residual contribution is

```fortran
rhs_ij(var_vpar) = rhs_ij(var_vpar) &
    + V_prof_pinch / sqrt(psi_grad2) * (ps0_x*vpar0_x + ps0_y*vpar0_y) &
            * r0 * v * BigR * xjac * tstep * factor(var_vpar,13)
```

and its two tangents repeat that sign:

```fortran
amat(var_vpar,var_rho)  = ... + V_prof_pinch / sqrt(psi_grad2) * (ps0_x*vpar0_x + ps0_y*vpar0_y) * rho  * v * BigR * xjac * theta * tstep
amat(var_vpar,var_vpar) = ... + V_prof_pinch / sqrt(psi_grad2) * (ps0_x*vpar_x  + ps0_y*vpar_y ) * r0   * v * BigR * xjac * theta * tstep
```

The JOREK convention is `AMAT = -theta*tstep*dB`, so both should be negative.
The density equation gets this right for the same pinch term: its residual
contribution is negative and `amat(var_rho,var_rho)` is positive.

**Suggested fix:** flip the sign of both `amat(var_vpar,...)` pinch terms.

**Note:** the residual sign itself also deserves a look. With
`V_pinch = -V_prof_pinch * grad(psi)/|grad(psi)|` as documented in the source
comment, the advection `rho * V_pinch . grad(v_par)` is
`-V_prof_pinch/|grad(psi)| * rho * grad(psi).grad(v_par)`, i.e. negative,
while the assembled term is positive. The generator follows the assembled
residual, so this note is a reading of the source comment rather than a
report difference.

**Scale:** 2 monomials in each of `amat(var_vpar,var_rho)` and
`amat(var_vpar,var_vpar)`, in each temperature branch.

---

## Finding 10 — the impurity parallel diffusivity loses its shock-capturing part in the toroidal channel

**Where:** `rhs_ij_k(var_rhoimp)`.

**What happens:** the parallel impurity diffusion term of the impurity-density
residual is split between the two FFT channels. The poloidal channel uses the
full diffusivity,

```fortran
rhs_ij(var_rhoimp) = &
    - ((D_par_local_imp+D_par_imp_sc_num*tau_sc)-D_prof_imp) * BigR / BB2 * Bgrad_rho_star   * (Bgrad_rhoimp) * xjac * tstep * factor(var_rhoimp,1) &
```

while the toroidal channel drops the shock-capturing contribution
`D_par_imp_sc_num*tau_sc`:

```fortran
rhs_ij_k(var_rhoimp) = &
    - (D_par_local_imp-D_prof_imp)                          * BigR / BB2 * Bgrad_rho_k_star * (Bgrad_rhoimp) * xjac * tstep * factor(var_rhoimp,1) &
```

Two other places in the same routine write the same physics with the full
coefficient, which is what makes this a slip rather than a convention:

* the density equation splits the identical impurity term across the same two
  channels and keeps `D_par_imp_sc_num*tau_sc` in both
  (`rhs_ij_k(var_rho)`);
* all three toroidal Jacobian blocks of the impurity equation —
  `amat_k(var_rhoimp,var_psi)`, `amat_k(var_rhoimp,var_rhoimp)` and
  `amat_kn(var_rhoimp,var_rhoimp)` — use the full coefficient, so the
  tangent does not differentiate the residual it belongs to.

**Suggested fix:** replace `(D_par_local_imp-D_prof_imp)` by
`((D_par_local_imp+D_par_imp_sc_num*tau_sc)-D_prof_imp)` in
`rhs_ij_k(var_rhoimp)`.

**Effect:** unlike findings 1, 2, 4, 5, 7 and 8, this one is in the residual,
so it changes the converged solution and not only the Newton convergence: with
shock capturing active (`tau_sc /= 0`) the toroidal part of the parallel
impurity diffusion is under-resolved relative to the poloidal part. It also
makes the impurity Jacobian inconsistent with its own residual.

**Scale:** 3 monomials, in each temperature branch.

---

## Reproducing

```bash
cd util/equation_codegen
python3 -m venv .venv && .venv/bin/python -m pip install -e .
.venv/bin/python examples/export_model600_terms.py
.venv/bin/python examples/diff_model600_reports.py
```

The **residual** printed under a `DIFF` block is the authoritative difference.
It is computed in the element basis, where the two spellings of a poloidal
bracket coincide, and displayed with `f_s -> f_x`, `f_t -> f_y`, `xjac -> 1`
so that it reads as a physical expression. Every line that reaches the
residual with a negative sign is a term the linearization produces and the
element routine does not.

In the aligned Markdown reports, a blank cell on the generated side is now a
genuine missing term: across all 162 blocks the only source lines without a
generated counterpart are the 7 of finding 3, the 6 of finding 6 and the 4 of
finding 9.

Residual composition of the parallel-velocity blocks, for orientation:

| Block | residual lines | by finding |
|---|---:|---|
| `rhs_ij(var_vpar)` | 8 | 6 |
| `amat(var_vpar,var_psi)` | 98 | 7 (84), 4/9 (6), 6 (8) |
| `amat_k(var_vpar,var_psi)` | 20 | 7 |
| `amat(var_vpar,var_rho)` | 6 | 6 (4), 9 (2) |
| `amat(var_vpar,var_vpar)` | 6 | 6 (4), 9 (2) |
| `amat_n(var_vpar,var_vpar)` | 2 | 8 |
| `amat_kn(var_vpar,var_vpar)` | 1 | 8 | A block reported with both `source` and `generated` lines and a non-zero
residual is a coefficient disagreement: `amat(var_u,var_T)` (finding 3), the
four `var_vpar` mass blocks (finding 6) and the two `var_vpar` pinch blocks
(finding 9).

A block reported as `SAME` differs only in whether a poloidal bracket is
written in element or physical coordinates; `amat(var_vpar,var_psi)` still
contains 16 such lines mixed in with its genuine differences.
