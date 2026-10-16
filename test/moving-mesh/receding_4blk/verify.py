#!/usr/bin/env python3
"""
Receding-surface conduction on a translating mesh, four blocks  (plan 05 gate
5.3 on the 4 x 1 split of ../receding/N080; plan 11 stage 2).

    usage: verify.py [case_dir]        (default: the current directory)

Run it from the case (or a copy of it) after the solve; it reads
OUTPUT/field.tec and INPUT/bc.txt there. The problem constants (alpha, the
recession rate v, the decay length d = alpha/v, the slab thickness L, Tw,
Tinf) and the exact cell average are imported from ../receding/generate.py,
which wrote N080.

Checks:

  * The physics checks of ../receding/verify.py that one grid can make, on
    the first j-row of EVERY zone of OUTPUT/field.tec (that script reads zone
    1 only; here the slab runs through four blocks along x, and the decay fit
    window 1 <= xi/d <= 5 spans blocks 1 and 2). The row is assembled block
    after block, in order of x, each shared node taken once:
      - DISCRETE DECAY LENGTH, fitted over 1 <= xi/d <= 5, against the closed
        form d_num/d = Pe/ln(1+Pe) of the steady three-term recurrence:
        relative difference <= 2 %.
      - POINTWISE ERROR against the exact cell averages, within the predicted
        budget 1.5 (Tw - Tinf) (d_num/d - 1) / e.
    The CONVERGENCE ORDER needs the grid family N040..N320 and stays with
    ../receding/verify.py; this case is one grid (N080's).

  * SEAMS. Every 101 record of INPUT/bc.txt pairs a cell face with its
    neighbour block's; the four corner nodes of the two faces must be
    identical in field.tec. The translation law moves every node by the same
    v dt, so the blocks agree on every shared node. Compared as written
    (15 significant digits, so equal text is necessary, not proof of bit
    equality).
"""

import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.dont_write_bytecode = True     # no __pycache__ beside ../receding/generate.py
sys.path.insert(0, os.path.join(HERE, "..", "receding"))
from generate import (ALPHA, DECAY, LX, TINF, TW, VEL,   # noqa: E402
                      t_exact_cell_average)

# Acceptance thresholds (../receding/verify.py).
DECAY_REL_TOL = 0.02              # measured vs closed-form discrete decay rate
POINTWISE_SLACK = 1.5             # budget = SLACK * the perturbation estimate
FIT_LO, FIT_HI = 1.0, 5.0         # fit window in units of the decay length


def read_tec_zones(path):
    """Every zone of a BLOCK-format Tecplot file -> (names, [(dims, {name: [values]})])."""
    with open(path) as fh:
        lines = fh.read().split("\n")
    names = [t.strip('"') for t in lines[0].split("=", 1)[-1].split() if t.strip('"')]
    starts = [n for n, line in enumerate(lines) if line.lstrip().startswith("ZONE")]
    if not starts:
        raise RuntimeError("no ZONE in %s" % path)
    zones = []
    for z, first in enumerate(starts):
        m = re.search(r"I=(\d+),\s*J=(\d+),\s*K=(\d+)", lines[first])
        if not m:
            raise RuntimeError("could not parse the dimensions of zone %d of %s" % (z + 1, path))
        I, J, K = (int(x) for x in m.groups())
        last = starts[z + 1] if z + 1 < len(starts) else len(lines)
        vals = [float(line) for line in lines[first + 1:last] if line.strip()]
        n_nodal = I * J * K
        n_cell = (I - 1) * (J - 1) * max(K - 1, 1)
        data, off = {}, 0
        for name in names:
            n = n_nodal if name in ("x", "y", "z") else n_cell
            data[name] = vals[off:off + n]
            off += n
        if off != len(vals):
            raise RuntimeError("zone %d of %s: %d values, %d expected" % (z + 1, path, len(vals), off))
        zones.append(((I, J, K), data))
    return names, zones


def first_row(zones):
    """Node x and cell T along the first j-row of the first k-plane, through
    every zone in order of x, each shared node once."""
    order = sorted(range(len(zones)), key=lambda z: zones[z][1]["x"][0])
    xs, temp = [], []
    for n, z in enumerate(order):
        (I, J, K), data = zones[z]
        row_x = data["x"][:I]
        if n > 0:
            if row_x[0] != xs[-1]:
                raise RuntimeError("block %d starts at x = %r, the previous one ends at %r"
                                   % (z + 1, row_x[0], xs[-1]))
            row_x = row_x[1:]
        xs += row_x
        temp += data["T"][:I - 1]
    return xs, temp


def fit_decay(xi_c, temp):
    """Least-squares slope of ln(T - T_base) against xi over the fit window
    (../receding/verify.py: T_base is the asymptote of the FINITE domain)."""
    e_l = math.exp(-LX / DECAY)
    t_base = TINF - (TW - TINF) * e_l / (1.0 - e_l)

    xs, ys = [], []
    for x, t in zip(xi_c, temp):
        if not (FIT_LO * DECAY <= x <= FIT_HI * DECAY):
            continue
        d = t - t_base
        if d <= 0.0:
            continue
        xs.append(x)
        ys.append(math.log(d))

    if len(xs) < 4:
        raise RuntimeError("too few points in the decay fit window")

    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den = sum((x - mx) ** 2 for x in xs)
    return -1.0 / (num / den), n


def face_nodes(i, j, k, f):
    """The four corner nodes (0-based) of face f of cell (i, j, k) (1-based),
    ordered by (first, second) tangential offset."""
    out = []
    for a in (0, 1):
        for c in (0, 1):
            if f in (1, 2):
                out.append((i - 1 if f == 1 else i, j - 1 + a, k - 1 + c))
            elif f in (3, 4):
                out.append((i - 1 + a, j - 1 if f == 3 else j, k - 1 + c))
            else:
                out.append((i - 1 + a, j - 1 + c, k - 1 if f == 5 else k))
    return out


def seam_check(zones, bc_path):
    """(pairs of nodes compared, mismatches, max |difference|) over every 101 record."""
    with open(bc_path) as fh:
        lines = fh.read().split("\n")
    two = (101, 103, 301, 302, 303, 304, 305, 306, 307)
    n, compared, bad, worst = 0, 0, 0, 0.0
    while n < len(lines) and lines[n].strip():
        b, i, j, k, f, ty = (int(t) for t in lines[n].split())
        if ty == 101:
            t = [int(x) for x in lines[n + 1].split()]
            bs, i_s, js, ks, fs = t[:5]
            if t[5:] != [1, 0, 0, 1]:
                raise RuntimeError("bc.txt record at line %d: orientation %s, only 1 0 0 1 is checked"
                                   % (n + 1, t[5:]))
            (I, J, K), A = zones[b - 1]
            (Is, Js, Ks), B = zones[bs - 1]
            for (p, q) in zip(face_nodes(i, j, k, f), face_nodes(i_s, js, ks, fs)):
                la = p[0] + I * (p[1] + J * p[2])
                lb = q[0] + Is * (q[1] + Js * q[2])
                for c in ("x", "y", "z"):
                    compared += 1
                    if A[c][la] != B[c][lb]:
                        bad += 1
                        worst = max(worst, abs(A[c][la] - B[c][lb]))
        n += 2 if ty in two else 1
    return compared, bad, worst


def main():
    case = sys.argv[1] if len(sys.argv) > 1 else "."
    solution = os.path.join(case, "OUTPUT", "field.tec")
    try:
        _, zones = read_tec_zones(solution)
    except (IOError, OSError):
        print("ERROR: file not found - %s" % solution)
        return 1

    try:
        xs, temp = first_row(zones)
    except RuntimeError as e:
        print("  FAIL: the blocks do not join along the row: %s" % e)
        print("  Result: FAIL")
        return 1
    nx = len(temp)
    if len(xs) != nx + 1:
        print("ERROR: %d nodes for %d cells along the row" % (len(xs), nx))
        return 1

    # Mesh-frame coordinate: the whole mesh translated rigidly, so subtracting
    # the first node recovers xi exactly, whatever time the run ended at.
    xi_n = [x - xs[0] for x in xs]
    xi_c = [0.5 * (xi_n[i] + xi_n[i + 1]) for i in range(nx)]
    exact = [t_exact_cell_average(xi_n[i], xi_n[i + 1]) for i in range(nx)]
    err = [t - e for t, e in zip(temp, exact)]
    linf = max(abs(e) for e in err)
    l2 = math.sqrt(sum(e * e for e in err) / nx)
    xi_max = xi_c[max(range(nx), key=lambda i: abs(err[i]))]

    d_num, n_fit = fit_decay(xi_c, temp)
    dx = LX / nx
    pe = VEL * dx / ALPHA
    d_pred = DECAY * pe / math.log(1.0 + pe)
    decay_rel = abs(d_num - d_pred) / d_pred
    eps = d_pred / DECAY - 1.0
    budget = POINTWISE_SLACK * (TW - TINF) * eps / math.e

    compared, bad, worst = seam_check(zones, os.path.join(case, "INPUT", "bc.txt"))

    print()
    print("=" * 74)
    print(" RECEDING-SURFACE CONDUCTION ON A TRANSLATING MESH, %d BLOCKS (gate 5.3)" % len(zones))
    print("=" * 74)
    print("  alpha            %.6e m2/s" % ALPHA)
    print("  recession rate   %.6e m/s" % VEL)
    print("  decay length d   %.6g m   (slab is %g d thick)" % (DECAY, LX / DECAY))
    print("  surface / far    %.2f / %.2f K" % (TW, TINF))
    print("  blocks / cells along x   %d / %d" % (len(zones), nx))
    print()
    print("  DISCRETE DECAY LENGTH  -- isolates the interior scheme")
    print("  %-6s %-9s %-13s %-13s %-13s %s" % ("N", "Pe", "d_measured", "d_predicted", "rel.diff", "pts"))
    print("  %-6d %-9.5f %-13.6e %-13.6e %-13.2e %d%s"
          % (nx, pe, d_num, d_pred, decay_rel, n_fit, "" if decay_rel <= DECAY_REL_TOL else "   <-- OFF"))
    print()
    print("  POINTWISE ERROR vs the analytic profile (cell averages)")
    print("  %-6s %-12s %-12s %-12s %-10s %s" % ("N", "max|dT| K", "L2 K", "budget K", "meas/pred", "at xi/d"))
    print("  %-6d %-12.4f %-12.4f %-12.4f %-10.3f %.2f%s"
          % (nx, linf, l2, budget, linf / (budget / POINTWISE_SLACK), xi_max / DECAY,
             "" if linf <= budget else "   <-- OVER"))
    print("  (the convergence order needs the grid family: ../receding/verify.py)")
    print()
    print("  seam node coordinates   %d compared, %d differ (max %.3e m)" % (compared, bad, worst))
    print()

    ok = True
    if decay_rel <= DECAY_REL_TOL:
        print("  PASS: discrete decay length matches Pe/ln(1+Pe) (rel. diff %.2e <= %.2f)"
              % (decay_rel, DECAY_REL_TOL))
    else:
        print("  FAIL: the discrete decay length does not match the closed-form")
        print("        prediction. The interior ALE/conduction balance is wrong.")
        ok = False
    if linf <= budget:
        print("  PASS: pointwise error within the predicted budget")
    else:
        print("  FAIL: pointwise error exceeds the predicted budget")
        ok = False
    if compared > 0 and bad == 0:
        print("  PASS: seams -- every node shared by two blocks is written alike")
    elif compared == 0:
        print("  FAIL: seams -- INPUT/bc.txt holds no 101 record: not a split case")
        ok = False
    else:
        print("  FAIL: seams -- blocks disagree on the position of a shared node")
        ok = False

    print()
    print("  Result: %s" % ("PASS" if ok else "FAIL"))
    print()
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
