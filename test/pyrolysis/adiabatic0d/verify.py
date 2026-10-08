#!/usr/bin/env python3
"""Gate 4.6 (verification record 2026-10-08): adiabatic 0-D decomposition against the exact type-1 ODE.

A 4-cell slab of real TACOT with every face adiabatic (BC 301, q = 0 on faces 1 and 2, 300 on 3/4) started at
T0 = 1000 K: the state stays uniform and follows the ODE of tools/ode0d.py (type1 variant) exactly; the only
FUSS error is the kinetics/energy splitting inside the RK3 stages (frozen stage temperature), second order in
dt here (0.051 K at vnn 0.5, 0.014 at 0.25, 0.0022 at 0.1, 0.0005 at 0.05; record 2026-10-08).

Checks over every OUTPUT/field<i>.tec:
  (a) max |T - T_or1|      <= TOL_T    (0.10 K; measured 0.051 K at vnn 0.5)
  (b) max |chi_r - chi_or1| <= TOL_CHI  (2e-3; measured 1.22e-3)
  (c) spatial non-uniformity of T and chi <= TOL_UNIF (1e-6; measured 2e-8)
  Reported, not gated: the exact type-2 (gas storage) minus type-1 difference (0.07 K here) and the FUSS error
  at the last file.

Usage: verify.py [case_dir]   Exit 0 on PASS, 1 on FAIL. Python 3.6, no numpy (imports ../tools/ode0d.py).
"""
import os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tools"))
import ode0d  # noqa: E402

TOL_T, TOL_CHI, TOL_UNIF = 0.10, 2.0e-3, 1.0e-6


def main():
    case = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "T1000")
    c = ode0d.Case(case, pato_tables=None)
    s = c.fuss_series()
    ts = [x[0] for x in s]
    T1, chi1 = c.integrate("type1", ts[-1], 2.0e-4, ts)
    T2, _ = c.integrate("type2", ts[-1], 2.0e-4, ts)
    eT = [s[i][1] - T1[i] for i in range(len(s))]
    ech = [max(abs(s[i][2][k] - chi1[i][k]) for k in range(c.nsub)) for i in range(len(s))]
    unif = max(x[3] for x in s)
    d2 = max((T2[i] - T1[i] for i in range(len(s))), key=abs)
    ok = True
    print("adiabatic0d: T0 = %g K, q = %g, %d files to %.3f s" % (c.T0, c.q, len(s), ts[-1]))
    print("  max |T - T_or1|   = %.4f K  (tol %g)  %s" % (max(abs(e) for e in eT), TOL_T, "PASS" if max(abs(e) for e in eT) <= TOL_T else "FAIL"))
    ok &= max(abs(e) for e in eT) <= TOL_T
    print("  max |chi - or1|   = %.2e    (tol %g)  %s" % (max(ech), TOL_CHI, "PASS" if max(ech) <= TOL_CHI else "FAIL"))
    ok &= max(ech) <= TOL_CHI
    print("  non-uniformity    = %.1e    (tol %g)  %s" % (unif, TOL_UNIF, "PASS" if unif <= TOL_UNIF else "FAIL"))
    ok &= unif <= TOL_UNIF
    print("  info: FUSS - or1 at t = %.2f s: %.4f K; type2 - type1 (exact): %.4f K; T_end FUSS %.3f or1 %.3f" % (ts[-1], eT[-1], d2, s[-1][1], T1[-1]))
    print("adiabatic0d: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
