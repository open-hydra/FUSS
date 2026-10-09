#!/usr/bin/env python3
"""Gate 4.7 (verification record 2026-10-08): heated 0-D decomposition against the exact type-1 ODE.

The adiabatic0d slab with a uniform volumetric source q = 3.75e8 W/m3 (INPUT/st.tec) from 300 K, 3 s: the
state heats at ~650 K/s (the rate of ATC1's 1 mm probe in its first second) through the whole decomposition
and follows tools/ode0d.py (type1) exactly. Run at vnn 0.00625 (dt ~ 2.1 ms = ATC1's step): the splitting
error is 0.0016 K (second order: 15.8 K at vnn 0.5, 3.8 at 0.25, 0.59 at 0.1, 0.14 at 0.05, 0.034 at 0.025,
0.008 at 0.0125; record 2026-10-08).

Checks over every OUTPUT/field<i>.tec:
  (a) max |T - T_or1|       <= TOL_T    (0.02 K; measured 0.002 K, oracle agreement with scipy RK45 1e-4 K)
  (b) max |chi_r - chi_or1| <= TOL_CHI  (1e-5; measured 1e-6)
  (c) spatial non-uniformity <= TOL_UNIF (1e-6; measured 1e-11)
  Reported, not gated: type2 - type1 (exact, -3.9 K at the end) and type1 with PATO's table hbar minus type1
  (the data difference, -0.15 K) when PATO's TACOT tables are on this machine.

Under gas-closure = type2 (plan 07 gate 4.1(d), test/pyrolysis/darcy/heated0d_type2: both 301 faces open at
p-wall, K = 1e-7 m2 so that the pressure excursion stays below 0.5 Pa and the fixed-p oracle applies) the
oracle is the type2 variant with the same TOL_T; the exact type-2 minus type-1 shift (-3.9 K at the end) is
the quantity the type-2 run must reproduce. TOL_UNIF is 1e-5 K there (the gas storage -eps p inside E follows
the sub-pascal pressure excursion; measured 6.5e-7 K) and TOL_CHI 2e-5 (measured 1.04e-5 at this step, 2.07e-6 at
half of it: second order; the RK4 oracle at dt 1e-4 is itself uncertain by 2.1e-6 in chi and 1.2e-3 K in T against
dt 5e-5; 2026-10-09).

Usage: verify.py [case_dir]   Exit 0 on PASS, 1 on FAIL. Python 3.6, no numpy (imports ../tools/ode0d.py).
FUSS_GATE_SELFTEST=1 sets every tolerance to 0 (the proven-red hook).
"""
import os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tools"))
import ode0d  # noqa: E402

TOL_T, TOL_CHI, TOL_UNIF = 0.02, 1.0e-5, 1.0e-6
TOL_UNIF_TYPE2 = 1.0e-5
TOL_CHI_TYPE2 = 2.0e-5


def gas_closure(case):
    for line in open(os.path.join(case, "input.ini")):
        s = line.split(";")[0].strip()
        if s.startswith("gas-closure") and "=" in s:
            return s.split("=", 1)[1].strip()
    return "type1"


def main():
    case = sys.argv[1] if len(sys.argv) > 1 else HERE
    closure = gas_closure(case)
    variant = "type2" if closure == "type2" else "type1"
    tol_T, tol_chi, tol_unif = TOL_T, (TOL_CHI_TYPE2 if variant == "type2" else TOL_CHI), (TOL_UNIF_TYPE2 if variant == "type2" else TOL_UNIF)
    if os.environ.get("FUSS_GATE_SELFTEST", "0") == "1":
        tol_T, tol_chi, tol_unif = 0.0, 0.0, 0.0
    c = ode0d.Case(case)
    if c.q <= 0.0:
        print("heated0d: FAIL  INPUT/st.tec missing or q <= 0 (q = %g)" % c.q)
        return 1
    s = c.fuss_series()
    ts = [x[0] for x in s]
    T1, chi1 = c.integrate("type1", ts[-1], 1.0e-4, ts)
    T2, chi2 = c.integrate("type2", ts[-1], 1.0e-4, ts)
    To, chio = (T2, chi2) if variant == "type2" else (T1, chi1)
    eT = [s[i][1] - To[i] for i in range(len(s))]
    ech = [max(abs(s[i][2][k] - chio[i][k]) for k in range(c.nsub)) for i in range(len(s))]
    unif = max(x[3] for x in s)
    ok = True
    print("heated0d (%s, oracle %s): T0 = %g K, q = %g W/m3, %d files to %.3f s, T_end FUSS %.3f or1 %.3f or2 %.3f"
          % (closure, variant, c.T0, c.q, len(s), ts[-1], s[-1][1], T1[-1], T2[-1]))
    print("  max |T - T_or%s|   = %.4f K  (tol %g)  %s" % (variant[-1], max(abs(e) for e in eT), tol_T, "PASS" if max(abs(e) for e in eT) <= tol_T else "FAIL"))
    ok &= max(abs(e) for e in eT) <= tol_T
    print("  max |chi - or%s|   = %.2e    (tol %g)  %s" % (variant[-1], max(ech), tol_chi, "PASS" if max(ech) <= tol_chi else "FAIL"))
    ok &= max(ech) <= tol_chi
    print("  non-uniformity    = %.1e    (tol %g)  %s" % (unif, tol_unif, "PASS" if unif <= tol_unif else "FAIL"))
    ok &= unif <= tol_unif
    d2 = max((T2[i] - T1[i] for i in range(len(s))), key=abs)
    line = "  info: type2 - type1 (exact) max %.3f K" % d2
    if c.pTv is not None:
        Tt, _ = c.integrate("type1-tablehbar", ts[-1], 1.0e-4, ts)
        line += "; PATO-table hbar minus type1 max %.3f K" % max((Tt[i] - T1[i] for i in range(len(s))), key=abs)
    print(line)
    print("heated0d: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
