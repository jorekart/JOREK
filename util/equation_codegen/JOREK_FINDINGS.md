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
rewrite is reported as `SAME`. At the time of writing, 174 of the
233 exported blocks are identical; the remaining 59 are the seventeen findings
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
| `Ti` | reproduced | findings 11 and 12 |
| `Te` | reproduced | findings 11, 12, 13, 14 and 15 |
| `T` | finding 16 | findings 11, 12, 13, 14, 15 and 17 |

Every source term of every row is reproduced by the generator. Apart from
findings 3, 6 and 9, the differences are all terms the generator produces and
the element routine does not.

Findings 3, 6, 9, 10, 11, 13, 14 and 16 are disagreements about a
coefficient, a power of `BigR`, a sign, a missing `theta` or a wrong variable
rather than omissions.  Findings 10 and 16 are the only two that sit in a
residual rather than in a tangent, so they are the only two that change the
converged solution.

The NEO branch is excluded from the default reports (`--include-neo` enables
it) and has not been audited here.

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

## Finding 12 — three smaller omissions in the energy tangents

All three are terms the generator produces and `mod_elt_matrix_fft.f90` does
not.

**(a) The kinetic-coupling term does not differentiate `BB2`.**
`amat(var_Ti,var_psi)` differentiates the flux dependence of `BB2` for the
friction sources,

```fortran
- v * ((GAMMA - 1.) / BigR) * vpar0**2 * (psi_x * ps0_x + psi_y * ps0_y) * ((r0+alpha_e*rimp0)*rn0*Sion_T) * xjac * theta * tstep &
- v * ((GAMMA - 1.) / BigR) * vpar0**2 * (psi_x * ps0_x + psi_y * ps0_y) * (particle_source + source_pellet + source_bg_drift + source_imp_drift) * xjac * theta * tstep &
```

but the residual carries one more term of exactly the same shape,

```fortran
+ (gamma-1.d0)*0.5d0 * v * aux_rho0 * vpar0**2 * BB2 * BigR * xjac * tstep * factor(var_Ti,16) &
```

whose `BB2_psi` contribution is absent. 4 monomials, and 4 more in
`amat(var_T,var_psi)`, which repeats the omission.

**(b) The perpendicular conductivity's density derivative is missing from the
toroidal channel.** `amat(var_Ti,var_rho)` carries both

```fortran
- dZKi_prof_drho * rho * BigR / BB2 * Bgrad_T_star * Bgrad_Ti  * xjac * theta * tstep &
+ dZKi_prof_drho * rho * BigR * (v_x*Ti0_x + v_y*Ti0_y)        * xjac * theta * tstep
```

while `amat_k(var_Ti,var_rho)` contains only `tgnum_Ti` terms, although
`rhs_ij_k(var_Ti)` has the matching `Bgrad_T_k_star` and `v_p*Ti0_p`
contributions. `ZKi_prof` depends on the density
(`ZKi_prof = get_zk_iperp(psi_norm) * max(r0,zkperp_density_floor)`), so the
`k` channel needs the same two terms. 4 monomials. The electron energy
equation repeats this exactly: `dZKe_prof_drho` is in `amat(var_Te,var_rho)`
and missing from `amat_k(var_Te,var_rho)`, and `dZK_prof_drho` likewise from
`amat_k(var_T,var_rho)` — another 4 monomials each.

**(c) The recombination sink does not differentiate its rate.** The residual
ends with

```fortran
- v * Ti0 * BigR * r0_corr * r0_corr  * Srec_T * xjac * tstep * factor(var_Ti,18)
```

and `Srec_T` is a function of the electron temperature — `amat(var_rho,var_Te)`
uses `dSrec_dT` for the same rate — but `amat(var_Ti,var_Te)` carries only the
`dSion_dT` friction terms. The missing entry is

```fortran
+ v * BigR * Ti0 * r0_corr * r0_corr * dSrec_dT * Te * xjac * theta * tstep
```

1 monomial. This is the same family as finding 5.

---

## Finding 13 — four ionization-energy tangent lines are missing their theta factor

**Where:** `amat(var_Te,var_Te)` and `amat_k(var_Te,var_Te)`, and the same
two blocks of the single-temperature equation, `amat(var_T,var_T)` and
`amat_k(var_T,var_T)` — the diffusive flux of the impurity ionization
potential energy.

**What happens:** every contribution to an `amat` block in this routine is
scaled by `theta * tstep`, the implicitness factor of the time scheme. These
four lines carry `tstep` alone:

```fortran
! amat(var_Te,var_Te)
+ (GAMMA - 1.) * dE_ion_dT * Te * ((D_par_local_imp+D_par_imp_sc_num*tau_sc)-D_prof_imp) * BigR / BB2 * Bgrad_rho_star * (Bgrad_rhoimp) * xjac * tstep &
+ (GAMMA - 1.) * dE_ion_dT * Te * D_prof_imp * BigR  * (v_x*(rimp0_x) + v_y*(rimp0_y)) * xjac * tstep &

! amat_k(var_Te,var_Te)
+ (GAMMA - 1.) * dE_ion_dT * Te * ((D_par_local_imp+D_par_imp_sc_num*tau_sc)-D_prof_imp) * BigR / BB2 * Bgrad_rho_k_star * (Bgrad_rhoimp) * xjac * tstep &
+ (GAMMA - 1.) * dE_ion_dT * Te * D_prof_imp * BigR  * (                v_p*(rimp0_p) /BigR**2 ) * xjac * tstep &
```

The check is exact: all 40 source monomials of `amat(var_Te,var_Te)` and all
20 of `amat_k(var_Te,var_Te)` that are not `tgnum_Te` terms become their
generated counterparts under the single substitution
`xjac*tstep -> xjac*theta*tstep`, and nothing else differs. The matching
residual terms in `rhs_ij(var_Te)` and `rhs_ij_k(var_Te)` are reproduced
exactly, so the residual is right and only the tangent is mis-scaled: those
entries are `1/theta` times too large.

**Suggested fix:** insert `* theta` on the four lines.

**Scale:** 60 monomials in the two-temperature branch and 60 in the
single-temperature one.

---

## Finding 14 — amat(var_Te,var_Te) contradicts amat_k(var_Te,var_Te) over alpha_e

**Where:** the `tgnum_Te` bracket of `amat(var_Te,var_Te)`, and the
`tgnum_T` bracket of `amat_k(var_T,var_psi)`.

**What happens:** the two blocks differentiate the same quantity and write its
toroidal part differently, two lines apart:

```fortran
! amat(var_Te,var_Te)
+ tgnum_Te * 0.25d0 / BigR * vpar0**2 &
    * Te * ((r0_x+alpha_e_bis*rimp0_x)*ps0_y - (r0_y+alpha_e_bis*rimp0_y)*ps0_x + F0 / BigR * (r0_p+alpha_e    *rimp0_p)) &
    * ( v_x * ps0_y -  v_y * ps0_x ) * xjac * theta * tstep * tstep &

! amat_k(var_Te,var_Te)
+ tgnum_Te * 0.25d0 / BigR * vpar0**2 &
    * Te * ((r0_x+alpha_e_bis*rimp0_x)*ps0_y - (r0_y+alpha_e_bis*rimp0_y)*ps0_x + F0 / BigR * (r0_p+alpha_e_bis*rimp0_p)) &
    * (                                 + F0 / BigR * v_p) * xjac * theta * tstep * tstep &
```

`alpha_e_bis` is the correct one. Differentiating
`d_a Pe = (r0_a + alpha_e*rimp0_a)*Te0 + (r0 + alpha_e_bis*rimp0)*Te0_a` with
respect to Te gives `(r0_a + alpha_e_bis*rimp0_a)*Te + ...` for every
direction `a`, because `alpha_e_bis = alpha_e + dalpha_e_dT*Te0`; the
poloidal parts of both blocks already use `alpha_e_bis`, and only the
toroidal part of the first disagrees.

The single-temperature equation has the mirror image of the same slip, with
the roles of the two channels swapped:

```fortran
! amat(var_T,var_psi)   -- correct
* (r0+alpha_imp_bis*rimp0) * (T0_x * psi_y - T0_y * psi_x) &
! amat_k(var_T,var_psi) -- alpha_imp where the residual uses alpha_imp_bis
* (r0+alpha_imp    *rimp0) * (T0_x * psi_y - T0_y * psi_x) &
```

The bracket multiplies a temperature gradient, so `alpha_imp_bis` is right;
the residual at `rhs_ij(var_T)` and the `p` channel of the same tangent both
use it.

**Suggested fix:** `alpha_e -> alpha_e_bis` in the first bracket and
`alpha_imp -> alpha_imp_bis` in the second.

**Scale:** 2 monomials each.

---

## Finding 15 — four omissions in the electron and total energy tangents

**(a) The ionization-energy terms that do not already carry `dE_ion_dT` are
not differentiated with respect to Te.** `amat(var_Te,var_Te)` differentiates
`E_ion` in the two residual terms where it multiplies a temperature gradient
and in the two diffusive-flux terms, but not in the six where it multiplies a
density gradient, the compression or the parallel flow:

```fortran
+ (GAMMA-1.) * v * E_ion * BigR**2 * (rimp0_s * u0_t - rimp0_t * u0_s)           * tstep * factor(var_Te,17)&
- (GAMMA-1.) * v * E_ion * F0 / BigR * Vpar0 * rimp0_p                    * xjac * tstep * factor(var_Te,17)&
- (GAMMA-1.) * v * E_ion * Vpar0 * (rimp0_s * ps0_t - rimp0_t * ps0_s)           * tstep * factor(var_Te,17)&
+ (GAMMA-1.) * v * E_ion * rimp0 * 2.d0 * BigR * u0_y                     * xjac * tstep * factor(var_Te,17)&
- (GAMMA-1.) * v * E_ion * rimp0 * (vpar0_s * ps0_t - vpar0_t * ps0_s)           * tstep * factor(var_Te,17)&
- (GAMMA-1.) * v * E_ion * rimp0 * F0 / BigR * vpar0_p                    * xjac * tstep * factor(var_Te,17)&
```

Each needs its `E_ion -> dE_ion_dT*Te` copy. 18 monomials, and 18 more in
`amat(var_T,var_T)`, which repeats the omission.

Note that the element routine supplies no second derivative of `E_ion`, so the
generator follows its convention and treats `dE_ion_dT` as frozen; the
`d2E_ion_dT2` contributions are therefore not counted here.

**(b) `amat(var_Te,var_rhon)` uses uncorrected densities.** The residual's
line radiation uses the corrected ones,

```fortran
- v * BigR * (r0_corr+alpha_e*rimp0_corr) * rn0_corr * LradDrays_T * xjac * tstep * factor(var_Te,13) &
```

but its neutral-density tangent drops the corrections:

```fortran
amat(var_Te,var_rhon) = ... + v * BigR * rhon * (r0 + rimp0 * alpha_e) * LradDrays_T * xjac * theta * tstep
```

The ionization sink two lines above it is consistent, because its residual
uses `r0` and `rn0` uncorrected. 2 monomials.

**(c) `amat(var_Te,var_rhoimp)` misses three radiation and ionization
derivatives.** The impurity density enters the electron density as
`r0 + alpha_e*rimp0` in the ionization sink and in the `LradDrays` and
`LradDcont` channels. The block carries the `Lrad` and `frad_bg` derivatives
but not those three. 5 monomials.

**(d) The friction sources do not differentiate `alpha_e`.** The released
kinetic energy carries the electron density `(r0+alpha_e*rimp0)`, and
`amat(var_T,var_T)` differentiates `Sion_T` in it but not `alpha_e`, while
the ionization sink two lines above does differentiate both. The impurity
column repeats the omission. 10 monomials in `amat(var_T,var_T)` and 11 in
`amat(var_T,var_rhoimp)`.

---

## Finding 16 — the single-temperature parallel convection uses the neutral density

**Where:** `rhs_ij(var_T)`, the parallel convection of the total pressure.

**What happens:** the poloidal half of the parallel pressure convection is
written with `rn0`, the neutral density, where every analogous line uses
`rimp0`:

```fortran
! rhs_ij(var_T), line ~2013
- v * (r0 + rn0  *alpha_imp_bis) * Vpar0 * (T0_s  * ps0_t - T0_t  * ps0_s) * tstep * factor(var_T,4 ) &

! the toroidal half of the same term, one line above
- v * (r0 + rimp0*alpha_imp_bis) * F0 / BigR * Vpar0 * T0_p         * xjac * tstep * factor(var_T,4 ) &

! the same line in the two-temperature equations
- v * (r0 + rimp0*alpha_i)     * Vpar0 * (Ti0_s * ps0_t - Ti0_t * ps0_s)   * tstep * factor(var_Ti,4) &
- v * (r0 + rimp0*alpha_e_bis) * Vpar0 * (Te0_s * ps0_t - Te0_t * ps0_s)   * tstep * factor(var_Te,4) &
```

`alpha_imp_bis` is the impurity closure coefficient, so pairing it with the
neutral density has no physical reading; the toroidal half of the very same
term, and both two-temperature equations, use `rimp0`.

**Suggested fix:** `rn0 -> rimp0`.

**Effect:** this is in the residual, so like finding 10 it changes the
converged solution, not only the Newton convergence — whenever neutrals and
impurities are both present and their densities differ.

**Scale:** 4 monomials.

---

## Finding 17 — amat_n(var_T,var_vpar) drops the impurity pressure

**Where:** `amat_n(var_T,var_vpar)`.

**What happens:** the residual convects the full pressure along the field,

```fortran
- v * (r0 + rimp0*alpha_imp) * T0 * GAMMA * F0 / BigR * vpar0_p * xjac * tstep * factor(var_T,3 ) &
```

but its toroidal parallel-velocity tangent keeps only the main-ion part:

```fortran
amat_n(var_T,var_vpar) = + v * r0 * GAMMA * T0 * F0 / BigR * vpar_p * xjac * theta * tstep &
```

Both two-temperature equations get this right —
`amat_n(var_Ti,var_vpar)` uses `(r0+rimp0*alpha_i)` and
`amat_n(var_Te,var_vpar)` uses `(r0+rimp0*alpha_e)` — which is what makes the
single-temperature one a slip rather than a convention.

**Suggested fix:** `r0 -> (r0 + rimp0*alpha_imp)`.

**Scale:** 2 monomials.

---

## Complete inventory

Every report line that is blank on one side, across all 233 blocks, belongs to
one of the findings above.  There is nothing unexplained left.

| Finding | Report lines | Blocks |
|---:|---:|---|
| 1 — impurity ion pressure in the diamagnetic tangents | 96 | `amat(var_rho,var_rhoimp)` 4, `amat(var_rho,var_t)` 2, `amat(var_rho,var_ti)` 2, `amat(var_u,var_rhoimp)` 44, `amat(var_u,var_t)` 22, `amat(var_u,var_ti)` 22 |
| 2 — impurity electron pressure in the induction tangents | 30 | `amat(var_psi,var_rhoimp)` 10, `amat(var_psi,var_t)` 8, `amat(var_psi,var_te)` 8, `amat_n(var_psi,var_rhoimp)` 2, `amat_n(var_psi,var_t)` 1, `amat_n(var_psi,var_te)` 1 |
| 4 — pinch not differentiated with respect to psi | 24 | `amat(var_rho,var_psi)` 12, `amat(var_vpar,var_psi)` 12 |
| 5 — alpha_e(T) in the density sources | 6 | `amat(var_rho,var_t)` 3, `amat(var_rho,var_te)` 3 |
| 6 — parallel-velocity time term | 36 | `amat(var_vpar,var_psi)` 8, `amat(var_vpar,var_rho)` 8, `amat(var_vpar,var_vpar)` 8, `rhs_ij(var_vpar)` 12 |
| 7 — BB2 in the tgnum_vpar tangent | 240 | `amat(var_vpar,var_psi)` 192, `amat_k(var_vpar,var_psi)` 48 |
| 8 — toroidal channel of the parallel-parallel viscosity | 6 | `amat_kn(var_vpar,var_vpar)` 2, `amat_n(var_vpar,var_vpar)` 4 |
| 10 — impurity parallel diffusivity, toroidal channel | 6 | `rhs_ij_k(var_rhoimp)` 6 |
| 12 — omissions in the energy tangents (a,b) | 21 | `amat(var_t,var_psi)` 4, `amat(var_ti,var_psi)` 4, `amat(var_ti,var_te)` 1, `amat_k(var_t,var_rho)` 4, `amat_k(var_te,var_rho)` 4, `amat_k(var_ti,var_rho)` 4 |
| 14 — alpha_e/alpha_imp versus their "bis" forms | 8 | `amat(var_te,var_te)` 4, `amat_k(var_t,var_psi)` 4 |
| 15 — omissions in the electron and total energy tangents | 70 | `amat(var_t,var_rhoimp)` 15, `amat(var_t,var_t)` 28, `amat(var_te,var_rhoimp)` 5, `amat(var_te,var_rhon)` 4, `amat(var_te,var_te)` 18 |
| 16 — neutral density in the single-T parallel convection | 8 | `rhs_ij(var_t)` 8 |
| 17 — impurity pressure in amat_n(var_T,var_vpar) | 2 | `amat_n(var_t,var_vpar)` 2 |
| **total** | **553** | |

Findings 3, 9, 11 and 13 do not appear in the table: their terms exist on both
sides and differ only by a coefficient, a sign, a power of `BigR` or a missing
`theta`, so the report puts them on one line rather than leaving a blank.  The
counts above are report lines and are not the residual monomial counts quoted
under each finding; several generated lines cancel against each other inside a
residual.

Only four source lines in the whole export have no generated partner at all:
the two `alpha_e` lines of finding 14 and the two uncorrected-density lines of
finding 15(b).  Both are symbol substitutions rather than scalings, so the
matcher deliberately does not pair them; guessing symbol equivalences would
hide exactly the kind of error they represent.

---

## Regression benchmark

The inventory above is frozen in
[`reference/model600_discrepancies.json`](reference/model600_discrepancies.json)
and checked by [`final_test.py`](final_test.py).  Each block of the reference
names the findings it belongs to, so a failure points straight back at this
document:

```bash
./final_test.py
```

Fixing any finding in the Fortran will make the benchmark fail with `-` lines
for the monomials that stopped disagreeing.  That is the intended signal:
re-run with `--update`, and strike the finding from this file.

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

The tool also verifies its own alignment: the blank cells of every block must
account for that block's multiset difference, and a surplus is reported as
`MISALIGNED`.  In the aligned Markdown reports, a blank cell on the generated
side is therefore a genuine missing term: across all 233 blocks every source line has a generated
counterpart on the same line, including the ones that differ only by a
coefficient (findings 3, 6 and 9) or by a power of `BigR` (finding 11).

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
