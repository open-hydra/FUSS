#!/usr/bin/env python3
"""Gate 4.2 of plan 07: the analytic steady pressure profile and its order of convergence.

Sealed at x = 0, p = p_w at x = L, constant pi, T, mu, W and K: with Gamma proportional to p the
steady equation is linear in p^2, p^2(x) = p_w^2 + (pi R T mu/(W K)) (L^2 - x^2). The solid cp of
1e12 freezes T; each material step is one Picard iteration on the lagged Gamma(p) (plan 07
section 3.3 lags the coefficient, so a single step would give the linear parabola, not this field).

Checks over the grids N020 N040 N080 N160 (OUTPUT/field.tec, OUTPUT/darcy.txt of each):
  (a) stationarity: the last step's dp_last <= TOL_STAT (1e-10: the run stops on the convergence
      residual sqrt(sum((p - p^n)/p-ref)^2) <= 1e-10, i.e. a last Picard increment of at most ~1e-10
      relative; the line solve's own round-off floor is ~5e-13 at 160 cells, measured) and the run
      stopped before its step cap
  (b) max |T - T0| <= 1e-6 K on every grid
  (c) the L-inf error of p against the exact profile decreases monotonically
  (d) the observed order on the finest pair >= ORDER_MIN (1.9: the central flux of p^2 and the
      half-cell Dirichlet closure are second order; measured ~2.0)

Usage: verify.py [parent_dir]   Exit 0 on PASS, 1 on FAIL. FUSS_GATE_SELFTEST=1 raises ORDER_MIN
to 1e9 and sets the other tolerances to 0 (the proven-red hook). Python 3.6, no numpy.
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import darcy_tools as dt  # noqa: E402

GRIDS = ["N020", "N040", "N080", "N160"]
TOL_STAT = 1.0e-10
TOL_T = 1.0e-6
ORDER_MIN = 1.9


def main():
    parent = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else HERE
    tol_stat, tol_t, order_min = TOL_STAT, TOL_T, ORDER_MIN
    if dt.selftest():
        tol_stat, tol_t, order_min = 0.0, 0.0, 1.0e9
    ok = True
    errs, dxs = [], []
    print("sealed_profile: grids %s under %s" % (" ".join(GRIDS), parent))
    for g in GRIDS:
        case = os.path.join(parent, g)
        meta = dt.case_meta(case)
        L, pw, pi, K, T0 = meta["L"], meta["p0"], meta["pi"], meta["K"], meta["T0"]
        W, mu, hg = dt.gas_constants(case)
        dims, t, d = dt.read_tec(os.path.join(case, "OUTPUT", "field.tec"))
        nx = dims[0] - 1
        xc = dt.cell_centres_x(dims, d)
        p, T = d["p_pg"], d["T"]
        err = max(abs(p[i] - dt.exact_profile(xc[i], L, pw, pi, T0, W, mu, K)) for i in range(nx))
        dT = max(abs(v - T0) for v in T)
        drows = dt.read_darcy(os.path.join(case, "OUTPUT", "darcy.txt"))
        cap = int(dt.ini_value(case, "iter-threshold"))
        stat = drows[-1][5]
        stopped = len(drows) < cap
        good = stat <= tol_stat and stopped and dT <= tol_t
        ok &= good
        errs.append(err)
        dxs.append(L / nx)
        print("  %s: %3d cells  steps %3d (cap %d)  last dp/p-ref %.2e  max|T-T0| %.1e K  max|p - exact| %.4e Pa  p(0) %.3f  %s"
              % (g, nx, len(drows), cap, stat, dT, err, p[0], "PASS" if good else "FAIL"))
    for a, b in zip(errs, errs[1:]):
        if not b < a:
            ok = False
            print("  (c) FAIL: the error does not decrease (%.4e -> %.4e)" % (a, b))
    orders = [math.log(errs[i] / errs[i + 1]) / math.log(dxs[i] / dxs[i + 1]) for i in range(len(errs) - 1)]
    print("  observed orders: %s" % "  ".join("%.3f" % o for o in orders))
    good = orders[-1] >= order_min
    ok &= good
    print("  (d) order on the finest pair %.3f (min %g)  %s" % (orders[-1], order_min, "PASS" if good else "FAIL"))
    print("sealed_profile: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
