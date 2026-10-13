#!/usr/bin/env python3
"""
Plan 08 section 3.5 (ii), FUSS standalone: the coupled ablating wall BC 307 with its exchange values from the TEST ONLY
coupled-env-file must reproduce BC 306 with lambda 0 and the same values as its env table, byte for byte -- one SEB
routine, two input paths; equal table rows interpolate to the row value exactly.

    python3 verify.py RUN307_DIR RUN306_DIR

Compares every file of OUTPUT/ (field files, seb.txt, wall files, bc306-state.dat, energy-budget.txt, darcy.txt, the
probes); the two runs must also write the same set of files. Exit 0 on PASS, 1 on FAIL. FUSS_GATE_SELFTEST=1 is applied
by the gate to the run307 input (h_r 1.5000001e6 in env307.dat), which must FAIL. Python 3.6.
"""
import filecmp
import os
import sys


def main():
    a, b = (os.path.join(os.path.abspath(d), "OUTPUT") for d in sys.argv[1:3])
    fa, fb = sorted(os.listdir(a)), sorted(os.listdir(b))
    bad = []
    if fa != fb:
        bad.append("file sets differ: only 307 %s, only 306 %s" % (sorted(set(fa) - set(fb)), sorted(set(fb) - set(fa))))
    same = [f for f in fa if f in fb and filecmp.cmp(os.path.join(a, f), os.path.join(b, f), shallow=False)]
    diff = [f for f in fa if f in fb and f not in same]
    if diff:
        bad.append("differ: " + " ".join(diff))
    for need in ("seb.txt", "field.tec", "bc306-state.dat"):
        if need not in fa:
            bad.append("missing " + need)
    print("coupled307: %d of %d files identical%s" % (len(same), len(fa), "" if not bad else " -- " + "; ".join(bad)))
    print("PASS" if not bad else "FAIL")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
