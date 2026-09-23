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
