#!/usr/bin/env python3
"""Gate 4.1(b) of plan 07: the steady gas mass balance on a 1-D line.

Face 1 sealed, face 2 open at p_w, a constant gas production pi in every cell: at steady state the
mass flux through the face at x is pi x whatever T(x) is (dm/dx = pi). The case runs 2 s
time-accurate (b2, quasi-steady = false; the gas relaxation time L^2/D_p is ~0.09 s) or as
Picard iterations on the lagged Gamma (b1, quasi-steady = true, the gate's sedded copy).

The face fluxes are rebuilt here from OUTPUT/field.tec with the solver's own closure: cell
Gamma_c = p_c W K/(R T_c mu), face Gamma_f = (Gamma_L + Gamma_R)/2, flux -Gamma_f (p_R - p_L)/dx;
open face Gamma_w = (Gamma_N + Gamma(p_w, T_N))/2, flux -Gamma_w (p_w - p_N)/(dx/2).

Checks:
  (i)   max_f |mdot_f - pi x_f| / (pi L) <= TOL_FLUX  (b1: 1e-9, round-off of the exact line solve;
        b2: 1e-6 = N tol, the residual bound of plan 07 section 4.1 (b))
  (ii)  the open-face flux from the energy budget's M_gas_out over the last write interval,
        |rate / A - pi L| / (pi L) <= TOL_FLUX
  (iii) energy identity of the budget at every row: |E(t) - E(0) + Q_out - Q_src| <= 1e-10 max(|E0|, |Q_out|)
        (the stage-1 E^n rule of plan 07 section 3.5: a wrong base state drifts by eps dp per step)
  (iv)  darcy.txt: every step converged (resid <= tol of the case) with exactly 1 sweep (1-D tridiag)

Usage: verify.py [case_dir]   Exit 0 on PASS, 1 on FAIL. FUSS_GATE_SELFTEST=1 sets every tolerance
to 0 (the proven-red hook). Python 3.6, no numpy.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import darcy_tools as dt  # noqa: E402

TOL_FLUX_QS = 1.0e-9
TOL_FLUX_TR = 1.0e-6
TOL_BAL = 1.0e-10


def main():
    case = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else HERE
    meta = dt.case_meta(case)
    quasi = (dt.ini_value(case, "quasi-steady") or "false").lower() == "true"
    tol_flux = TOL_FLUX_QS if quasi else TOL_FLUX_TR
    tol_bal = TOL_BAL
    tol_solver = float(dt.ini_value(case, "tol"))
    if dt.selftest():
        tol_flux, tol_bal, tol_solver = 0.0, 0.0, 0.0
    L, pw, pi, K = meta["L"], meta["p0"], meta["pi"], meta["K"]
    W, mu, hg = dt.gas_constants(case)
    dims, t, d = dt.read_tec(os.path.join(case, "OUTPUT", "field.tec"))
    nx, ny = dims[0] - 1, dims[1] - 1
    xc = dt.cell_centres_x(dims, d)
    p, T = d["p_pg"], d["T"]
    dx = L / nx
    ok = True

    def check(label, val, tol):
        nonlocal ok
        good = val <= tol
        ok &= good
        print("  %-52s %.3e  (tol %g)  %s" % (label, val, tol, "PASS" if good else "FAIL"))

    print("mass_balance_1d (%s): %d cells, t = %g s, pi = %g kg/m3/s, K = %g m2, p_w = %g Pa"
          % ("b1 quasi-steady" if quasi else "b2 transient", nx, t, pi, K, pw))
    # (i) face fluxes on the first j-row of cells (every row is identical in 1-D)
    worst, x_worst = 0.0, None
    for i in range(nx - 1):
        gL = dt.gamma(p[i], T[i], W, mu, K)
        gR = dt.gamma(p[i + 1], T[i + 1], W, mu, K)
        mdot = -0.5 * (gL + gR) * (p[i + 1] - p[i]) / dx
        xf = 0.5 * (xc[i] + xc[i + 1])
        err = abs(mdot - pi * xf)
        if err > worst:
            worst, x_worst = err, xf
    gN = dt.gamma(p[nx - 1], T[nx - 1], W, mu, K)
    gw = dt.gamma(pw, T[nx - 1], W, mu, K)
    mdot_w = -0.5 * (gN + gw) * (pw - p[nx - 1]) / (0.5 * dx)
    err_w = abs(mdot_w - pi * L)
    if err_w > worst:
        worst, x_worst = err_w, L
    print("  face fluxes: open face %.9e vs pi L %.9e kg/m2/s; worst face at x = %s" % (mdot_w, pi * L, x_worst))
    check("(i)   max_f |mdot_f - pi x_f| / (pi L)", worst / (pi * L), tol_flux)
    # (ii) budget rate of the open face
    rows = dt.read_budget(os.path.join(case, "OUTPUT", "energy-budget.txt"))
    A = 0.0
    I, J, K3 = dims
    # open-face area: the i = nx face of every (j, k) cell, from the node coordinates
    for k in range(max(K3 - 1, 1)):
        for j in range(J - 1):
            n = lambda jj, kk: (nx) + jj * I + kk * I * J
            y0, y1 = d["y"][n(j, k)], d["y"][n(j + 1, k)]
            z0, z1 = d["z"][n(j, k)], d["z"][n(j, k + 1)]
            A += abs(y1 - y0) * abs(z1 - z0)
    if len(rows) >= 2:
        rate = (rows[-1][5] - rows[-2][5]) / (rows[-1][0] - rows[-2][0])
        print("  budget: M_gas_out rate %.9e kg/s over %.3f..%.3f s, A = %g m2 -> %.9e kg/m2/s"
              % (rate, rows[-2][0], rows[-1][0], A, rate / A))
        check("(ii)  |M_gas_out rate / A - pi L| / (pi L)", abs(rate / A - pi * L) / (pi * L), tol_flux)
    else:
        print("  (ii)  SKIP: fewer than two budget rows")
    # (iii) energy identity
    E0 = rows[0][1]
    worst_bal = 0.0
    for r in rows:
        worst_bal = max(worst_bal, abs(r[1] - E0 + r[3] - r[6]) / max(abs(E0), abs(r[3]), 1e-300))
    check("(iii) energy identity, max over %d rows" % len(rows), worst_bal, tol_bal)
    # (iv) solver health
    drows = dt.read_darcy(os.path.join(case, "OUTPUT", "darcy.txt"))
    check("(iv)  max resid over %d steps (solver tol)" % len(drows), max(r[4] for r in drows), tol_solver)
    check("(iv)  max |sweeps - 1| per step", max(abs(r[3] - 1) for r in drows), 0.0 if not dt.selftest() else -1.0)
    print("  info: p(0) = %.3f Pa, p(L) cell = %.3f Pa, T range %.6f..%.6f K" % (p[0], p[nx - 1], min(T), max(T)))
    print("mass_balance_1d: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
