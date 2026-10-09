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

Under gas-closure = type2 (plan 07 gate 4.1(d), test/pyrolysis/darcy/adiabatic0d_type2: both 301 faces open
at p-wall through open-bc-types, K = 1e-7 m2 so that the pressure excursion stays below 0.1 Pa and the
fixed-p oracle applies) the oracle is the type2 variant, with the same TOL_T; the exact type-2 minus
type-1 shift (+0.070 K) is printed as the quantity the type-2 run must reproduce. TOL_UNIF is 1e-5 K
there: the gas storage -eps p inside E makes T follow the sub-pascal pressure excursion (measured 1e-7 K).

Usage: verify.py [case_dir]   Exit 0 on PASS, 1 on FAIL. Python 3.6, no numpy (imports ../tools/ode0d.py).
FUSS_GATE_SELFTEST=1 sets every tolerance to 0 (the proven-red hook).
"""
import os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tools"))
import ode0d  # noqa: E402

TOL_T, TOL_CHI, TOL_UNIF = 0.10, 2.0e-3, 1.0e-6
TOL_UNIF_TYPE2 = 1.0e-5


def gas_closure(case):
    for line in open(os.path.join(case, "input.ini")):
        s = line.split(";")[0].strip()
        if s.startswith("gas-closure") and "=" in s:
            return s.split("=", 1)[1].strip()
    return "type1"


def main():
    case = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "T1000")
    closure = gas_closure(case)
    variant = "type2" if closure == "type2" else "type1"
    tol_T, tol_chi, tol_unif = TOL_T, TOL_CHI, (TOL_UNIF_TYPE2 if variant == "type2" else TOL_UNIF)
    if os.environ.get("FUSS_GATE_SELFTEST", "0") == "1":
        tol_T, tol_chi, tol_unif = 0.0, 0.0, 0.0
    c = ode0d.Case(case, pato_tables=None)
    s = c.fuss_series()
    ts = [x[0] for x in s]
    T1, chi1 = c.integrate("type1", ts[-1], 2.0e-4, ts)
    T2, chi2 = c.integrate("type2", ts[-1], 2.0e-4, ts)
    To, chio = (T2, chi2) if variant == "type2" else (T1, chi1)
    eT = [s[i][1] - To[i] for i in range(len(s))]
    ech = [max(abs(s[i][2][k] - chio[i][k]) for k in range(c.nsub)) for i in range(len(s))]
    unif = max(x[3] for x in s)
    d2 = max((T2[i] - T1[i] for i in range(len(s))), key=abs)
    ok = True
    print("adiabatic0d (%s, oracle %s): T0 = %g K, q = %g, %d files to %.3f s" % (closure, variant, c.T0, c.q, len(s), ts[-1]))
    print("  max |T - T_or%s|   = %.4f K  (tol %g)  %s" % (variant[-1], max(abs(e) for e in eT), tol_T, "PASS" if max(abs(e) for e in eT) <= tol_T else "FAIL"))
    ok &= max(abs(e) for e in eT) <= tol_T
    print("  max |chi - or%s|   = %.2e    (tol %g)  %s" % (variant[-1], max(ech), tol_chi, "PASS" if max(ech) <= tol_chi else "FAIL"))
    ok &= max(ech) <= tol_chi
    print("  non-uniformity    = %.1e    (tol %g)  %s" % (unif, tol_unif, "PASS" if unif <= tol_unif else "FAIL"))
    ok &= unif <= tol_unif
    print("  info: FUSS - or%s at t = %.2f s: %.4f K; type2 - type1 (exact): %.4f K; T_end FUSS %.3f or1 %.3f or2 %.3f"
          % (variant[-1], ts[-1], eT[-1], d2, s[-1][1], T1[-1], T2[-1]))
    print("adiabatic0d: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
