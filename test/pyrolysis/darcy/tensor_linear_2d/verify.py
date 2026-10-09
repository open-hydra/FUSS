#!/usr/bin/env python3
"""Gate 4.3 of plan 07: orthotropic permeability on a 2-D square.

K = R diag(K1, K2) R^T with R the 45-degree rotation in x-y, p_w on every face from the field
p^2 = p0^2 + a x. Because Gamma is proportional to p, the steady equation is div(K grad p^2) = 0
and the p^2-linear field is its exact solution; the cell velocity is u = -(K/(mu eps)) grad p, so
|u_y| / |u_x| = K_yx / K_xx in every cell. The red control (K diagonal, --control) must give u_y = 0 to the
solver tolerance: max |u_y| / max |u_x| <= TOL_RATIO. Plan 07 asked for exactly 0; that is unattainable here --
the y-faces are Dirichlet with x-dependent p_w and the damped Jacobi iterates pass through y-dependent states, so
the converged field keeps a round-off y-variation (measured 3e-13 relative, 2026-10-08).

Checks (OUTPUT/field.tec, OUTPUT/darcy.txt):
  (a) max_cells | |u_y| / |u_x| - K_yx/K_xx | <= TOL_RATIO (1e-10, plan 07 section 4.3; --control: max |u_y| / max |u_x|)
  (b) max |p - exact| / p0 printed; gated at TOL_P (1e-9: the discrete field is exact to the sweep tolerance)
  (c) stationarity: the last step's dp_last <= TOL_STAT (1e-10, the convergence stop at a residual
      root of 1e-10 as in sealed_profile) and the run stopped before its step cap
  (d) max |T - T0| <= 1e-6 K

Usage: verify.py [case_dir] [--control]   Exit 0 on PASS, 1 on FAIL. FUSS_GATE_SELFTEST=1 makes every
tolerance negative (the proven-red hook). Python 3.6, no numpy.
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import darcy_tools as dt  # noqa: E402

TOL_RATIO = 1.0e-10
TOL_P = 1.0e-9
TOL_STAT = 1.0e-10
TOL_T = 1.0e-6


def main():
    args = [a for a in sys.argv[1:] if a != "--control"]
    control = "--control" in sys.argv[1:]
    case = os.path.abspath(args[0]) if args else HERE
    tol_ratio, tol_p, tol_stat, tol_t, tol0 = TOL_RATIO, TOL_P, TOL_STAT, TOL_T, 0.0
    if dt.selftest():
        tol_ratio, tol_p, tol_stat, tol_t, tol0 = -1.0, -1.0, -1.0, -1.0, -1.0
    meta = dt.case_meta(case)
    Kt = meta["K"]
    Kxx, Kyx = Kt[0], Kt[3]
    p0, a, T0 = meta["p0"], meta["a"], meta["T0"]
    dims, t, d = dt.read_tec(os.path.join(case, "OUTPUT", "field.tec"))
    nx, ny = dims[0] - 1, dims[1] - 1
    xc = dt.cell_centres_x(dims, d)
    p, T, ux, uy = d["p_pg"], d["T"], d["u_pg_1"], d["u_pg_2"]
    ok = True

    def check(label, val, tol):
        nonlocal ok
        good = val <= tol
        ok &= good
        print("  %-48s %.3e  (tol %g)  %s" % (label, val, tol, "PASS" if good else "FAIL"))

    print("tensor_linear_2d%s: %d x %d cells, t = %g, K_xx = %g, K_yx = %g, a = %g Pa^2/m"
          % (" (red control)" if control else "", nx, ny, t, Kxx, Kyx, a))
    if control:
        check("(a) max |u_y| / max |u_x|", max(abs(v) for v in uy) / max(abs(v) for v in ux), tol_ratio)
    else:
        target = Kyx / Kxx
        worst = max(abs(abs(uy[c]) / abs(ux[c]) - target) for c in range(len(p)))
        print("  target |u_y|/|u_x| = %.8f" % target)
        check("(a) max | |u_y|/|u_x| - K_yx/K_xx |", worst, tol_ratio)
    err_p = max(abs(p[c] - math.sqrt(p0 * p0 + a * xc[c])) for c in range(len(p))) / p0
    check("(b) max |p - exact| / p0", err_p, tol_p)
    drows = dt.read_darcy(os.path.join(case, "OUTPUT", "darcy.txt"))
    cap = int(dt.ini_value(case, "iter-threshold"))
    check("(c) last dp_last / p-ref", drows[-1][5], tol_stat)
    check("(c) steps - (cap - 1), stopped before the cap", float(len(drows) - (cap - 1)), tol0)
    check("(d) max |T - T0| (K)", max(abs(v - T0) for v in T), tol_t)
    print("  info: steps %d, sweeps max %d median %d, |u_x| range %.3e..%.3e m/s"
          % (len(drows), max(r[3] for r in drows), sorted(r[3] for r in drows)[len(drows) // 2],
             min(abs(v) for v in ux), max(abs(v) for v in ux)))
    print("tensor_linear_2d: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
