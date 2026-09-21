#!/bin/bash
#
# Regression test for the equation_codegen tool.
#
# Unlike the JOREK non-regression tests under reg_tests/ (which compile
# binaries, run an MPI job and diff HDF5 restart files), this test is a pure
# Python check: it recomputes every model-600 discrepancy between the
# hand-written element routine (models/model600/mod_elt_matrix_fft.f90) and
# the symbolic linearization (util/equation_codegen/src/jorek_equations),
# and compares it against the frozen reference in
# reference/model600_discrepancies.json.  A change to either side that shifts
# which terms disagree — a real regression, not just a reordering — makes
# this test fail and print exactly which monomials are no longer accounted
# for.
#
# Follows the same pass/fail convention as reg_tests/run_test.sh: exit 0 and
# a line "Test 'equation_codegen' passed." on success, non-zero and a
# detailed failure report otherwise.  Meant to be run directly, or as its own
# step in CI — it does not go through reg_tests/run_test.sh.
#
# Usage:
#   util/equation_codegen/run_test.sh              # check the current tree
#   util/equation_codegen/run_test.sh --update      # rewrite the reference
#                                                    # after an intended change

set -u

TESTNAME="equation_codegen"
STARTDIR=$(readlink -f "$(dirname "$0")")
cd "$STARTDIR" || exit 1

# --- Set up the virtual environment the first time this runs.
if [ ! -x ".venv/bin/python" ]; then
    echo "Setting up .venv for ${TESTNAME} ..."
    python3 -m venv .venv                          || exit 1
    .venv/bin/python -m pip install --quiet -e .    || exit 1
fi

# --- Run the regression check.  final_test.py prints, on failure, exactly
#     which monomials each disagreeing block is missing or has gained.
.venv/bin/python final_test.py "$@"
status=$?

if [ $status -eq 0 ]; then
    echo "Test '${TESTNAME}' passed."
else
    echo "Test '${TESTNAME}' FAILED."
fi

exit $status
