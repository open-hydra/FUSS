#!/usr/bin/env python3
"""The Klinkenberg oracle of plan 07 (closed 2026-10-13): the sealed_profile cases with [FUSS-Darcy] klinkenberg = true.
FUSS's form is thesis eq. 3.15, K_eff = K + (mu/p) B (B in m2/s; the literature's K (1 + beta/p) with beta in Pa is the
same law with beta = mu B / K at constant mu). With Gamma = rho_g K_eff/mu the mass flux is -(W K/(R T mu)) (p + b) dp/dx,
b = mu B / K, so p + b obeys the plain equation exactly:
    (p + b)^2 = (p_w + b)^2 + (pi R T mu/(W K)) (L^2 - x^2)
(Wu, Pruess & Persoff 1998 eq. 3.2 for the steady form; derived for FUSS's tensor form in the 2026-10-10 benchmark survey).
Checks over the grids N020 N040 N080 N160 staged by gate_darcy_klinkenberg (each with B_v = B_c = B I and the key on):
  (a) the L-inf error against the Klinkenberg profile decreases monotonically, order on the finest pair >= ORDER_MIN (1.9)
  (b) the run is not the plain profile: max|p_exact,b - p_exact,0| on the finest grid >= 100 x its error (the term acts)
Usage: verify_klinkenberg.py <parent>   (b is read from <parent>/klinkenberg.json, written by the gate).
FUSS_GATE_SELFTEST=1 raises ORDER_MIN to 1e9. Python 3.6, no numpy.
"""
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import darcy_tools as dt  # noqa: E402

GRIDS = ["N020", "N040", "N080", "N160"]
ORDER_MIN = 1.9


def exact(x, L, pw, pi, T, W, mu, K, b):
    a = pi * dt.R_UNIV * T * mu / (W * K)
    return math.sqrt((pw + b) ** 2 + a * (L * L - x * x)) - b


def main():
    parent = os.path.abspath(sys.argv[1])
    b = json.load(open(os.path.join(parent, "klinkenberg.json")))["b"]
    order_min = 1.0e9 if dt.selftest() else ORDER_MIN
    ok = True
    errs, dxs, sep = [], [], 0.0
    print("klinkenberg: b = mu B/K = %.6g Pa, grids %s under %s" % (b, " ".join(GRIDS), parent))
    for g in GRIDS:
        case = os.path.join(parent, g)
        meta = dt.case_meta(case)
        L, pw, pi, K, T0 = meta["L"], meta["p0"], meta["pi"], meta["K"], meta["T0"]
        W, mu, hg = dt.gas_constants(case)
        dims, t, d = dt.read_tec(os.path.join(case, "OUTPUT", "field.tec"))
        nx = dims[0] - 1
        xc = dt.cell_centres_x(dims, d)
        p = d["p_pg"]
        err = max(abs(p[i] - exact(xc[i], L, pw, pi, T0, W, mu, K, b)) for i in range(nx))
        sep = max(abs(exact(x, L, pw, pi, T0, W, mu, K, b) - dt.exact_profile(x, L, pw, pi, T0, W, mu, K)) for x in xc)
        errs.append(err)
        dxs.append(L / nx)
        print("  %s: %3d cells  max|p - exact_b| %.4e Pa  p(0) %.3f (exact %.3f; without Klinkenberg %.3f)"
              % (g, nx, err, p[0], exact(xc[0], L, pw, pi, T0, W, mu, K, b), dt.exact_profile(xc[0], L, pw, pi, T0, W, mu, K)))
    for a, c in zip(errs, errs[1:]):
        if not c < a:
            ok = False
            print("  (a) FAIL: the error does not decrease (%.4e -> %.4e)" % (a, c))
    orders = [math.log(errs[i] / errs[i + 1]) / math.log(dxs[i] / dxs[i + 1]) for i in range(len(errs) - 1)]
    good = orders[-1] >= order_min
    ok &= good
    print("  (a) observed orders %s; finest pair %.3f (min %g)  %s" % (" ".join("%.3f" % o for o in orders), orders[-1], order_min,
                                                                     "PASS" if good else "FAIL"))
    good = sep >= 100.0 * errs[-1]
    ok &= good
    print("  (b) Klinkenberg moves the profile by %.4e Pa, %.0f x the finest error  %s" % (sep, sep / errs[-1], "PASS" if good else "FAIL"))
    print("klinkenberg: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
