#!/usr/bin/env python3
"""Compare the two model-600 term reports block by block.

Meld shows every line that differs, including the many lines that differ only
because source and generated monomials are listed in a different order.  This
script instead compares each assignment block as a multiset of monomials and,
where the multisets differ, prints the algebraic residual ``source-generated``.
A block whose residual is zero contains the same physics in a different
arrangement; a block with a non-zero residual is a real disagreement.
"""

from collections import Counter, OrderedDict
from pathlib import Path
import argparse
import re
import sys

import sympy as sp

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def read_blocks(path):
    """Return {(section, assignment, repeat): [monomial, ...]}."""

    blocks = OrderedDict()
    key = None
    section = ""
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith("## "):
            section = line[3:].strip()
            continue
        heading = re.match(r"^#### `(.*)`", line)
        if heading:
            repeats = sum(
                1 for existing in blocks
                if existing[0] == section and existing[1] == heading.group(1)
            )
            key = (section, heading.group(1), repeats)
            blocks[key] = []
            continue
        if key is not None and line.strip():
            blocks[key].append(line.strip())
    return blocks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    reports = PROJECT_ROOT / "reports"
    parser.add_argument("source", nargs="?", default=reports / "model600_fortran_terms.md")
    parser.add_argument("generated", nargs="?", default=reports / "model600_generated_terms.md")
    parser.add_argument(
        "--quiet", action="store_true", help="print only the per-block verdicts",
    )
    args = parser.parse_args()

    source = read_blocks(args.source)
    generated = read_blocks(args.generated)
    differing = 0
    keys = list(OrderedDict.fromkeys(list(source) + list(generated)))
    for key in keys:
        source_only = Counter(source.get(key, ())) - Counter(generated.get(key, ()))
        generated_only = Counter(generated.get(key, ())) - Counter(source.get(key, ()))
        label = "{} {} #{}".format(*key)
        if not source_only and not generated_only:
            print("OK   {}  ({} terms)".format(label, len(source.get(key, ()))))
            continue
        differing += 1
        residual = sp.expand(
            sum(sp.sympify(term) for term in source_only.elements())
            - sum(sp.sympify(term) for term in generated_only.elements())
        )
        print("DIFF {}: source-only={} generated-only={} residual terms={}".format(
            label, sum(source_only.values()), sum(generated_only.values()),
            0 if residual == 0 else len(sp.Add.make_args(residual)),
        ))
        if args.quiet:
            continue
        for term in sorted(source_only.elements()):
            print("   source    {}".format(term))
        for term in sorted(generated_only.elements()):
            print("   generated {}".format(term))
        if residual == 0:
            print("   residual: 0 (the block agrees; only the ordering differs)")
        else:
            print("   residual (source-generated):")
            for term in sp.Add.make_args(residual):
                print("     {}".format(term))
    print("\nBlocks with differences: {} / {}".format(differing, len(keys)))
    return 1 if differing else 0


if __name__ == "__main__":
    sys.exit(main())
