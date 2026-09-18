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
`u`, `zj`, and `w` rows and writes two line-aligned Markdown reports:

```bash
.venv/bin/python examples/check_model600_integrated.py
meld reports/model600_fortran_terms.md reports/model600_generated_terms.md
```

The integrated script accepts the same selector, for example
`--equation psi`; with no selector it exports all four rows.

To generate only the reports, without running the pass/fail checks:

```bash
.venv/bin/python examples/export_model600_terms.py
```

For a faster report while developing one row, select it explicitly:

```bash
.venv/bin/python examples/export_model600_terms.py --equation psi
```

The available selections are `psi`, `u`, `zj`, and `w`. Repeat
`--equation` to select more than one row. If no selection is supplied, all
four rows are exported.

Assignments and terms follow their order in
`models/model600/mod_elt_matrix_fft.f90`. A source term that is not available
from the equation generator has an empty line at the same position in the
generated report. Generated-only terms are appended to their corresponding
assignment block with an empty source line. NEO AMAT terms remain explicitly
source-only until their symbolic linearization is enabled.
