#!/usr/bin/env python3
"""
R5's type-2 variant (plan 08 section 2.9 (3)): the two sides of every block-connection (101) face put opposite
conduction, gas mass and gas enthalpy flows into their cells -- one flux, conserved across the unequal-spacing interface.

    python3 verify.py [RUN_DIR]

From OUTPUT/interface-flux.txt (rows `time b i j k f bs is js ks fs q_if qm_if qh_if`, every solution output): for each
face and its partner at the same time, |q_A + q_B| over the run, per quantity, relative to the largest |q| of that
quantity over the run, <= 1e-12. Measured 2026-10-13: q_if 1.6e-13 (the conduction kernel's round-off), qm_if and qh_if
exactly 0 (both sides evaluate one expression on the same data). Each quantity must also be nonzero, so the pairing is
not vacuous. Before plan 08 the Darcy row used the cell's own metric at a full-cell distance: the two sides'
conductances differed by the spacing ratio, 4 here. FUSS_GATE_SELFTEST=1 sets the bound to 0. Exit 0 on PASS.
"""
import os
import sys

TOL = 0.0 if os.environ.get("FUSS_GATE_SELFTEST") == "1" else 1.0e-12


def main():
    run = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
    by = {}
    for line in open(os.path.join(run, "OUTPUT", "interface-flux.txt")):
        if line.startswith("#") or not line.strip():
            continue
        r = line.split()
        by[(r[0], tuple(int(x) for x in r[1:6]))] = (tuple(int(x) for x in r[6:11]), [float(x) for x in r[11:14]])
    ok = True
    npair = sum(1 for (t, key), (part, v) in by.items() if (t, part) in by)
    for q, name in enumerate(("q_if", "qm_if", "qh_if")):
        worst, mx = 0.0, 0.0
        for (t, key), (part, v) in by.items():
            if (t, part) not in by:
                continue
            worst = max(worst, abs(v[q] + by[(t, part)][1][q]))
            mx = max(mx, abs(v[q]))
        rel = worst / mx if mx > 0 else float("inf")
        good = mx > 0 and rel <= TOL
        ok &= good
        print("(pair)  %-5s max|q_A + q_B| %.3e of max|q| %.3e: %.3e (<= %.0e)  %s"
              % (name, worst, mx, rel, TOL, "PASS" if good else "FAIL"))
    print("%d face rows paired; %s" % (npair, "PASS" if ok else "FAIL"))
    return 0 if ok and npair > 0 else 1


if __name__ == "__main__":
    sys.exit(main())
