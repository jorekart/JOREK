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
