#!/usr/bin/env python3
"""
Conjugate two-layer slab, unequal cell spacing across the interface.

WHAT THIS GATES
---------------
That FUSS's block-to-block interface flux uses the NEIGHBOUR's geometry, not
just the local side's. `BC_Connection` builds the interface temperature as

    T_int = (k1 T1 D1 + k2 T2 D2) / (k1 D1 + k2 D2)

with D2 taken from bc%Mg -- the only place in the solver that field is read.
For a two-layer slab that expression is exactly the flux-continuity interface
temperature between cell centres at dx1/2 and dx2/2, so a correct scheme
reproduces the analytic profile to round-off.

WHY THE TOLERANCE IS SO TIGHT
-----------------------------
The steady solution of a two-layer slab is piecewise LINEAR, and a linear field
is exactly representable on this grid: the interior Laplacian of a linear
function is zero, and the half-cell Dirichlet closure is exact for it too. So
there is no discretisation error to hide in. The only thing that can put a kink
in the profile is the interface treatment, which makes the expected result
machine precision rather than a budget.

The tolerance below is therefore set by the RESIDUAL threshold of the run
(res-threshold = 1e-10), not by any discretisation argument. Observed on the
correct code: 4.5e-09 K against a 700 K span.

WHY dx2/dx1 = 4 MATTERS, AND WHY THIS CASE HAD TO BE WRITTEN
------------------------------------------------------------
With equal spacing either side, D1 = D2 and the neighbour's metric cancels out
of the weighting entirely -- the formula collapses to a pure conductivity
blend. Every other multi-block case in this suite has matched spacing, so all
of them pass whether or not bc%Mg is correct. Halving every bc%Mg on
multimat_Plate leaves it bit-identical over 214 steps; doing the same here
moves the interface temperature by 28 K.

That is the whole point of the case: it is the configuration in which the
defect is visible, and it is the configuration real TPS stacks are in.
"""

import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from generate import (K1, K2, L1, NX1, NX2, Q, T_COLD, T_HOT,   # noqa: E402
                      T_INT, DX1, DX2, t_exact)

SOLUTION = os.path.join(HERE, "OUTPUT", "field.tec")

# Set by the run's residual threshold (1e-10), not by discretisation: the exact
# solution is representable on this grid. Two orders above what the correct
# code achieves, and nine below what a wrong interface metric produces.
TOL_K = 1.0e-6


def read_zones(path):
    lines = open(path).read().split("\n")
    names = [t.strip('"') for t in lines[0].split("=", 1)[-1].split() if t.strip('"')]  # VARIABLES names, quoted (ORION <= v1.6) or bare (ORION >= v1.7.0 writes `VARIABLES = x y z T ...`)
    starts = [i for i, l in enumerate(lines) if l.strip().startswith("ZONE")]

    zones = []
    for zi, i in enumerate(starts):
        m = re.search(r"I=(\d+),\s*J=(\d+),\s*K=(\d+)", lines[i])
        ni, nj, nk = (int(x) for x in m.groups())
        end = starts[zi + 1] if zi + 1 < len(starts) else len(lines)

        vals = []
        for line in lines[i + 1:end]:
            s = line.strip()
            if not s:
                continue
            try:
                vals.append(float(s))
            except ValueError:
                pass

        n_nodal = ni * nj * nk
        n_cell = (ni - 1) * (nj - 1) * max(nk - 1, 1)
        d, off = {}, 0
        for nm in names:
            n = n_nodal if nm in ("x", "y", "z") else n_cell
            d[nm] = vals[off:off + n]
            off += n
        zones.append((ni, d))
    return zones


def main():
    try:
        zones = read_zones(SOLUTION)
    except FileNotFoundError:
        print("ERROR: file not found - %s" % SOLUTION)
        return 1

    if len(zones) != 2:
        print("ERROR: expected 2 blocks, found %d" % len(zones))
        return 1

    print()
    print("=" * 70)
    print(" CONJUGATE TWO-LAYER SLAB, UNEQUAL SPACING ACROSS THE INTERFACE")
    print("=" * 70)
    print("  k1 / k2            %.1f / %.1f  W/m/K" % (K1, K2))
    print("  dx1 / dx2          %.3e / %.3e m   (ratio %g)" % (DX1, DX2, DX2 / DX1))
    print("  interface weights  k1/dx1 = %.0f  vs  k2/dx2 = %.0f   (%.0f:1)"
          % (K1 / DX1, K2 / DX2, (K1 / DX1) / (K2 / DX2)))
    print("  analytic q         %.6f W/m2" % Q)
    print("  analytic T_int     %.6f K" % T_INT)
    print()

    worst = 0.0
    edge = {}
    for zi, (ni, d) in enumerate(zones):
        xs = d["x"][:ni]                       # first i-line of nodes
        temp = d["T"][:ni - 1]                 # first j row, first k plane
        xc = [0.5 * (xs[i] + xs[i + 1]) for i in range(ni - 1)]
        err = [t - t_exact(x) for t, x in zip(temp, xc)]
        e = max(abs(v) for v in err)
        worst = max(worst, e)
        print("  block %d   %2d cells   max |T - T_exact| = %.4e K"
              % (zi + 1, ni - 1, e))
        edge[zi] = (temp[-1], temp[0])

    # Reconstruct the interface temperature the scheme must have produced, from
    # the two cells either side. Reported separately because it is the specific
    # quantity bc%Mg controls -- a whole-field norm can dilute it.
    t1 = edge[0][0]      # last cell of block 1
    t2 = edge[1][1]      # first cell of block 2
    t_int_num = (K1 * t1 / DX1 + K2 * t2 / DX2) / (K1 / DX1 + K2 / DX2)

    print()
    print("  T last cell blk1   %.9f K" % t1)
    print("  T first cell blk2  %.9f K" % t2)
    print("  implied T_int      %.9f K   (analytic %.9f)" % (t_int_num, T_INT))
    print("  interface error    %+.4e K" % (t_int_num - T_INT))
    print()
    print("  tolerance          %.1e K  (set by res-threshold, not by dx)" % TOL_K)

    ok = worst < TOL_K and abs(t_int_num - T_INT) < TOL_K
    if ok:
        print()
        print("  Result: PASS")
    else:
        print()
        print("  Result: FAIL")
        print("  The exact solution is piecewise linear and representable on this")
        print("  grid, so this is not a resolution problem. An error of order 10 K")
        print("  here means the interface temperature is being built from the")
        print("  wrong side's metric -- check bc%Mg reaches BC_Connection.")
    print()
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
