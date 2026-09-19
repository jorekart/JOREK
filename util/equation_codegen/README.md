# JOREK equation code generator

This package implements the symbolic DSL described in
[`SPECIFICATION.md`](SPECIFICATION.md). The initial version provides fields,
test functions, structured spatial derivatives, external-function dependency
policies, directional linearization, and the JOREK RHS/AMAT sign convention.

Model-199 equations, FFT-channel classification, and Fortran emission will be
added after the symbolic core has been validated.

Equation definitions use bare fields (`psi`, `u`, `rho`, ...). The
linearization stage decides when those fields become current/background values
and when the selected field becomes a trial function or increment. Explicit
forms such as `psi.current` remain available for low-level expressions and
tests.

## Setup and tests

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -e .
.venv/bin/python -m unittest discover -s tests -v
```

To print the current symbolic and JOREK-Fortran-style model-199 equation-1
linearization:

```bash
.venv/bin/python examples/print_model199_equation1.py
```

To compare generated equation-1 terms with the current model-199 Fortran
source:

```bash
.venv/bin/python examples/compare_model199_equation1.py
```

The comparison is symbolic rather than textual. A mismatch raises an error
that includes the source expression, generated expression, and their
algebraic difference.

## Model 600 reports

The integrated model-600 checker validates the currently implemented `psi`,
`u`, `zj`, `w` and `rho` rows and writes two line-aligned Markdown reports:

```bash
.venv/bin/python examples/check_model600_integrated.py
meld reports/model600_fortran_terms.md reports/model600_generated_terms.md
```

The integrated script accepts the same selector, for example
`--equation psi`; with no selector it exports all five rows.

To generate only the reports, without running the pass/fail checks:

```bash
.venv/bin/python examples/export_model600_terms.py
```

For a faster report while developing one row, select it explicitly:

```bash
.venv/bin/python examples/export_model600_terms.py --equation psi
```

To inspect only the factored perpendicular-momentum RHS:

```bash
.venv/bin/python examples/export_model600_terms.py \
  --equation u --assignment 'rhs_ij(var_u)'
```

This focused mode preserves the Fortran outer-term structure and does not
build the expensive `u` AMAT columns.

The available selections are `psi`, `u`, `zj`, `w` and `rho`. Repeat
`--equation` to select more than one row. If no selection is supplied, all
five rows are exported.

To compare the two reports block by block instead of line by line:

```bash
.venv/bin/python examples/diff_model600_reports.py
.venv/bin/python examples/diff_model600_reports.py --quiet   # verdicts only
```

Each assignment block is compared as a multiset of monomials, and a block that
differs is printed together with the algebraic residual `source-generated`.
Reordering alone therefore does not show up as a difference.

## Temperature conventions

The element routine builds every pressure in `construct_pressure`, which is
written once for both temperature models and always uses the species
temperatures `Ti0`/`Te0` and the per-species impurity coefficients `alpha_i`,
`alpha_e` and `alpha_e_bis`. The one-temperature model evolves the *total*
temperature and enters that routine with

```
Ti0 = Te0 = T0/2
```

so a single-species temperature is half of the evolved field, and the stored
one-temperature closure values are the means

```
alpha_imp     = (alpha_i + alpha_e)/2
alpha_imp_bis = (alpha_i + alpha_e_bis)/2
alpha_imp_tri = alpha_e_tri/4
```

The comparison therefore expands both reports in the two-species basis. The
DSL follows the same convention: `_diamagnetic_pressure` and the induction
equation's electron pressure use `T/2` in the one-temperature branch. This
replaced two earlier "legacy tangent" overrides that reproduced JOREK's
factor of one half by hand; that half is physical, not a legacy quirk.

## Known JOREK differences

With the conventions above, every `psi`, `u`, `zj`, `w` and `rho` source term
is reproduced by the generator, and 81 of the 100 exported assignment blocks
are identical. The remaining 19 are terms the generator produces and the element
routine does not, or coefficients that disagree. They are recorded, with their
suggested Fortran fixes, in [`JOREK_FINDINGS.md`](JOREK_FINDINGS.md):

1. the impurity part of the ion pressure is not differentiated in the
   diamagnetic terms of the momentum and density equations
   (`amat(var_u,var_Ti)`, `amat(var_u,var_T)`, `amat(var_u,var_rhoimp)` and
   the same three columns of `var_rho`);
2. the impurity part of the electron pressure is not differentiated in the
   induction equation (`amat(var_psi,var_Te)`, `amat_n(var_psi,var_Te)`,
   `amat(var_psi,var_rhoimp)`, and an `amat_n(var_psi,var_rhoimp)` that the
   element routine does not have at all);
3. the one-temperature `amat(var_u,var_T)` diamagnetic-viscosity tangent is a
   factor two too large;
4. the density inward-pinch term is not differentiated with respect to `psi`
   (`amat(var_rho,var_psi)`);
5. the temperature dependence of `alpha_e` is not differentiated in the
   density ionization/recombination sources (`amat(var_rho,var_Te)`,
   `amat(var_rho,var_T)`), although the momentum equation does differentiate
   it in the same sources.

Assignments and terms follow their order in
`models/model600/mod_elt_matrix_fft.f90`. A source term that is not available
from the equation generator has an empty line at the same position in the
generated report. Generated-only terms are appended to their corresponding
assignment block with an empty source line. NEO AMAT terms remain explicitly
source-only until their symbolic linearization is enabled.
