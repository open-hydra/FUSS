#!/usr/bin/env python3
"""Two-material slab in one block under pyrolysis (test/pyrolysis/tools/make_two_layer.py) against the exact
composite-slab solution (test/pyrolysis/tools/two_layer_exact.py).

    verify.py <case>                        the error of every OUTPUT/field<n>.tec: max over the cells of |T - T_exact|
                                            at the cell centres, and its time; exit 0 when the run reached its end
                                            time with T inside [T0 - 1, TW + 1]
    verify.py --decay <case>                the same table; exit 0 when the run reached its end time with a finite
                                            field and the error of its last output is below that of its first (the
                                            check of the runs at the nominal vnn limits, below)
    verify.py --ladder <c1> <c2> [<c3> ...] the errors of N, 2N, 4N... cells at the output nearest T_SAMPLE (each at
                                            its own time), coarsest first; exit 0 by the convergence rule of
                                            tools/conv_rule.py with p = 2: the error falls on every rung and the
                                            observed order on the finest pair lies in [1.6, 2.4]

WHAT THE ORDER BAND IS. The scheme is the second-order finite-volume diffusion operator; the jump lies on a face, and
inside each layer the solution is smooth at T_SAMPLE (8 s: 46 cell diffusion times h^2/alpha of the coarsest mesh,
0.175 s for LOW at N = 20, after the step at t = 0; the heated front, sqrt(alpha_LOW t) = 3.4 mm, has reached the
interface 5 mm from the heated face). With the two half cells in series at the jump the interface flux is
consistent to O(h^2) and the error falls at second order: p = 2, and the band is the convergence rule's
[0.8 p, 1.2 p] = [1.6, 2.4] (tools/conv_rule.py, the rule of every exact-solution gate). The mean of the two k_s at
the jump (the scheme before 2026-11-26) gives the interface an O(1) relative conductance error at every mesh (10.5
against 1.9 times k_LOW), a contact resistance that shrinks with h only: the error falls at FIRST order, outside the
band. Neither bound is taken from a FUSS run. At a fixed vnn dt ~ h^2, so the time error is O(h^6) for
RK3, O(h^4) for RK2 and O(h^2) for euler: the band holds for every scheme at any vnn below the limit (for euler the
error constant grows with vnn).

AT THE NOMINAL vnn LIMIT (0.5 euler and RK2, 0.628 RK3) the shortest mode of the scheme is not damped: its
amplification factor has modulus 1 (euler 1 - z at z = 2, RK3 1 - z + z^2/2 - z^3/6 at z = 2.51), for one material as
for two. The step at t = 0 leaves an oscillation from cell to cell near the heated face that decays only through the
slightly smaller modes (+-8 K for RK3, +-15 K for euler at 8 s on 40 cells, the same on one material), and at the RK3
limit the field leaves [T0, TW] by a few K early on (no maximum principle there). So a nominal-limit run is graded on
stability only: it reaches its end and its error falls between the first output and the last (an unstable mode
would grow). Accuracy is graded below the limit (--ladder).

Every number read is printed. FUSS_GATE_SELFTEST=1 empties the order window of --ladder (the proven-red hook of
conv_rule.py). Python 3.6, no numpy.
"""
import glob
import json
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tools"))
import two_layer_exact  # noqa: E402
import conv_rule as cr  # noqa: E402

T_SAMPLE = 8.0
P_EXPECTED = 2


def read_field(path):
    """[(xc list, T list)] per zone of a FUSS field*.tec (x of the cell centres along i, j = k = 1 row first),
    and the SOLUTIONTIME"""
    lines = open(path).read().split("\n")
    names = [t.strip('"') for t in re.split(r"[\s,]+", lines[0].split("=", 1)[-1]) if t.strip('"')]
    starts = [i for i, l in enumerate(lines) if l.strip().startswith("ZONE")]
    zones, time = [], None
    for zi, s in enumerate(starts):
        m = re.search(r"I=(\d+),\s*J=(\d+),\s*K=(\d+)", lines[s])
        ni, nj, nk = (int(v) for v in m.groups())
        mt = re.search(r"SOLUTIONTIME=\s*([-+0-9.eEdD]+)", lines[s])
        time = float(mt.group(1).replace("D", "E"))
        end = starts[zi + 1] if zi + 1 < len(starts) else len(lines)
        vals = []
        for l in lines[s + 1:end]:
            for tok in l.split():
                try:
                    vals.append(float(tok))
                except ValueError:
                    pass
        nn, nc = ni * nj * nk, (ni - 1) * (nj - 1) * max(nk - 1, 1)
        d, off = {}, 0
        for nm in names:
            n = nn if nm in ("x", "y", "z") else nc
            d[nm] = vals[off:off + n]
            off += n
        x = d["x"]
        xc = [0.5 * (x[i] + x[i + 1]) for i in range(ni - 1)]
        rows = [d["T"][r * (ni - 1):(r + 1) * (ni - 1)] for r in range((nj - 1) * max(nk - 1, 1))]
        zones.append((xc, rows))
    return zones, time


def fields(case):
    out = []
    for f in glob.glob(os.path.join(case, "OUTPUT", "field*.tec")):
        z, t = read_field(f)
        out.append((t, f, z))
    return sorted(out)


def error(case, ex, z, t):
    worst, where = 0.0, None
    for xc, rows in z:
        for row in rows:
            for x, T in zip(xc, row):
                if not math.isfinite(T):
                    return float("inf"), x
                e = abs(T - ex.T(x, t))
                if e > worst:
                    worst, where = e, x
    return worst, where


def one(case):
    cfg = json.load(open(os.path.join(case, "case.json")))
    ex = two_layer_exact.from_case(cfg)
    fl = fields(case)
    if not fl:
        print("no OUTPUT/field*.tec in %s" % case)
        return 1
    print("%s: %d cells, %s vnn %g (%s)" % (case, cfg["ncells"], cfg["time_scheme"], cfg["vnn"],
                                            "type2" if cfg["type2"] else "type1"))
    ok = True
    for t, f, z in fl:
        e, x = error(case, ex, z, t)
        Tmax = max(max(max(r) for r in rows) for _xc, rows in z)
        Tmin = min(min(min(r) for r in rows) for _xc, rows in z)
        print("  t = %10.6f s  max |T - T_exact| = %.6e K at x = %s   T in [%.3f, %.3f] K  (%s)"
              % (t, e, "%.6e" % x if x is not None else "-", Tmin, Tmax, os.path.basename(f)))
        if not math.isfinite(e) or Tmin < cfg["T0"] - 1.0 or Tmax > cfg["TW"] + 1.0:
            ok = False
    tend = None
    for l in open(os.path.join(case, "input.ini")):
        if l.startswith("time-threshold"):
            tend = float(l.split("=")[1])
    reached = abs(fl[-1][0] - tend) <= 1e-9 * tend or fl[-1][0] >= tend
    print("  last field at t = %.6f s of %g s: %s; bounded in [T0 - 1, TW + 1]: %s"
          % (fl[-1][0], tend, "reached" if reached else "NOT reached", "yes" if ok else "NO"))
    return 0 if (ok and reached) else 1


def decay(case):
    """--decay: reached the end, finite, last error below the first"""
    cfg = json.load(open(os.path.join(case, "case.json")))
    ex = two_layer_exact.from_case(cfg)
    fl = [q for q in fields(case) if q[0] > 0.0]
    if len(fl) < 2:
        print("%s: fewer than two field files after t = 0" % case)
        return 1
    errs = []
    for t, f, z in fl:
        e, x = error(case, ex, z, t)
        errs.append(e)
        print("  t = %10.6f s  max |T - T_exact| = %.6e K" % (t, e))
    tend = None
    for l in open(os.path.join(case, "input.ini")):
        if l.startswith("time-threshold"):
            tend = float(l.split("=")[1])
    reached = abs(fl[-1][0] - tend) <= 1e-9 * tend or fl[-1][0] >= tend
    ok = reached and all(math.isfinite(e) for e in errs) and errs[-1] < errs[0]
    print("  %s: end %s; error %.4e K at t = %.3f s -> %.4e K at t = %.3f s: %s"
          % (os.path.basename(case.rstrip("/")), "reached" if reached else "NOT reached", errs[0], fl[0][0], errs[-1],
             fl[-1][0], "decays" if ok else "DOES NOT DECAY"))
    return 0 if ok else 1


def ladder(cases):
    errs = []
    for c in cases:
        cfg = json.load(open(os.path.join(c, "case.json")))
        ex = two_layer_exact.from_case(cfg)
        # the field file nearest T_SAMPLE (a run writes at the first step past each multiple of sol-dtime): the
        # exact solution is evaluated at that file's own time
        hit = sorted((abs(t - T_SAMPLE), t, z) for t, _f, z in fields(c))
        if not hit or hit[0][0] > 0.01 * T_SAMPLE:
            print("%s: no field file within 1 %% of t = %g s" % (c, T_SAMPLE))
            return 1
        e, x = error(c, ex, hit[0][2], hit[0][1])
        errs.append((cfg["ncells"], e, x))
        print("  N = %4d  max |T - T_exact| at t = %.9f s: %.6e K (x = %.6e m)" % (cfg["ncells"], hit[0][1], e, x))
    if len(errs) < 2:
        print("two_layer ladder: FAIL  the convergence rule needs a ladder of 2 or more runs")
        return 1
    for (n1, e1, _), (n2, e2, _) in zip(errs, errs[1:]):
        print("  order %4d -> %4d: %.3f  (printed; the rule grades the finest pair)" % (n1, n2, cr.order(e1, e2, n1, n2)))
    ok, text = cr.check("max |T - T_exact| at 8 s [K]", [e for _n, e, _x in errs], [n for n, _e, _x in errs],
                        P_EXPECTED, fmt="%.6e")
    print("  %s  %s" % (text, "PASS" if ok else "FAIL"))
    print("two_layer ladder: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--ladder":
        sys.exit(ladder(sys.argv[2:]))
    if len(sys.argv) > 2 and sys.argv[1] == "--decay":
        sys.exit(decay(sys.argv[2]))
    sys.exit(one(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "N020")))
