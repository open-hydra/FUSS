#!/usr/bin/env python3
"""Gate 4.1(a) of plan 07: a sealed, source-free slab keeps its pressure.

Every face impermeable, no gas production, uniform T: the pressure equation has the solution
p = p-init exactly, every face flux is zero and the solver must not move anything. The run ends
at its first material step: nothing moves, so the convergence residual (the T, chi and p
increments) is exactly zero and the solver stops as converged -- one step is the whole test.

Checks (OUTPUT/field.tec of the last step, OUTPUT/darcy.txt):
  (a) max |p_pg - p0| / p0 <= TOL_P   (1e-14: round-off of a zero-residual sweep; measured 0)
  (b) max |p_i - p_j| over cells == 0 exactly (every face difference, hence every face flux, is 0)
  (c) max |mdot_pg| == 0 and max |u_pg| == 0 exactly (the cell flux of the solved field)
  (d) darcy.txt: every step converged with resid == 0 and at most 1 sweep (the mandatory one)
  (e) max |T - T0| <= 1e-9 (no gas work: p constant)

Usage: verify.py [case_dir]   Exit 0 on PASS, 1 on FAIL. FUSS_GATE_SELFTEST=1 makes every tolerance
negative so that the exact zeros fail (the proven-red hook). Python 3.6, no numpy.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import darcy_tools as dt  # noqa: E402

TOL_P = 1.0e-14
TOL_T = 1.0e-9


def main():
    case = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else HERE
    tol_p, tol_t, tol0 = TOL_P, TOL_T, 0.0
    if dt.selftest():
        tol_p, tol_t, tol0 = -1.0, -1.0, -1.0
    meta = dt.case_meta(case)
    p0, T0 = meta["p0"], meta["T0"]
    dims, t, d = dt.read_tec(os.path.join(case, "OUTPUT", "field.tec"))
    p = d["p_pg"]
    ok = True

    def check(label, val, tol):
        nonlocal ok
        good = val <= tol
        ok &= good
        print("  %-44s %.3e  (tol %g)  %s" % (label, val, tol, "PASS" if good else "FAIL"))

    print("sealed_nosource: %d cells, t = %g s, p0 = %g Pa, T0 = %g K" % (len(p), t, p0, T0))
    check("(a) max |p - p0| / p0", max(abs(v - p0) for v in p) / p0, tol_p)
    check("(b) max |p_i - p_j| (uniformity, Pa)", max(p) - min(p), tol0)
    check("(c) max |mdot_pg| (kg/m2/s)", max(abs(v) for v in d["mdot_pg"]), tol0)
    check("(c) max |u_pg| (m/s)", max(abs(v) for c in ("u_pg_1", "u_pg_2", "u_pg_3") for v in d[c]), tol0)
    rows = dt.read_darcy(os.path.join(case, "OUTPUT", "darcy.txt"))
    check("(d) max resid over %d steps" % len(rows), max(r[4] for r in rows), tol0)
    check("(d) max sweeps per step - 1", max(r[3] for r in rows) - 1, tol0)
    check("(e) max |T - T0| (K)", max(abs(v - T0) for v in d["T"]), tol_t)
    print("sealed_nosource: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
