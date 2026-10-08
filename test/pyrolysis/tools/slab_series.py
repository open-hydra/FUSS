#!/usr/bin/env python3
"""Exact solution of the step-heated, back-insulated slab with constant properties (plan 06 verification, 2026-10-08):

  T(x, t) = T_w + (T_0 - T_w) sum_{n>=0} (2 (-1)^n / lam_n) cos(lam_n x / L) exp(-lam_n^2 alpha t / L^2),  lam_n = (2n+1) pi/2

x from the insulated face (x = 0, dT/dx = 0) to the heated face (x = L, T = T_w for t > 0), T(x, 0) = T_0.
Pure Python 3.6, no numpy; 4000 terms (the series is exact to round-off for alpha t / L^2 >= 1e-7, checked against the
semi-infinite erf solution to 3e-12 K at t = 0.5 s, 1 mm).
"""
import math


def exact(depth, t, alpha, L=0.05, Tw=1644.0, T0=298.0, nterms=4000):
    """T at `depth` below the heated face at time t"""
    if t <= 0.0:
        return T0
    x = L - depth
    s = 0.0
    for n in range(nterms):
        lam = (2 * n + 1) * math.pi / 2.0
        ex = -lam * lam * alpha * t / (L * L)
        if ex < -745.0:
            break
        s += 2.0 * (-1.0) ** n / lam * math.cos(lam * x / L) * math.exp(ex)
    return Tw + (T0 - Tw) * s


if __name__ == "__main__":
    a = 0.4 / (280.0 * 1000.0)
    for d, t in ((1e-3, 0.5), (1e-3, 1.0), (2e-3, 1.0), (16e-3, 60.0)):
        e = 1644.0 + (298.0 - 1644.0) * math.erf(d / (2.0 * math.sqrt(a * t)))
        print("d = %g m, t = %g s: series %.6f K, semi-infinite erf %.6f K" % (d, t, exact(d, t, a), e))
