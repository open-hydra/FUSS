#!/usr/bin/env python3
"""Exact solution of the step-heated, back-insulated composite slab (two layers, constant properties).

    layer 1: 0 <= x <= A (k1, w1 = rho1 cp1, alpha1 = k1 / w1), insulated at x = 0
    layer 2: A <= x <= L (k2, w2, alpha2), T = TW at x = L for t > 0; T(x, 0) = T0; T and k dT/dx continuous at A

    T(x, t) = TW + (T0 - TW) sum_n c_n X_n(x) exp(-lam_n^2 t)
    X_n = cos(b1 x) on layer 1, B sin(b2 (L - x)) on layer 2,   b_i = lam / sqrt(alpha_i)
    eigenvalues: f(lam) = k1 b1 sin(b1 A) sin(b2 C) - k2 b2 cos(b1 A) cos(b2 C) = 0,   C = L - A
    B = cos(b1 A) / sin(b2 C) (or k1 b1 sin(b1 A) / (k2 b2 cos(b2 C)) where sin(b2 C) is small)
    c_n = int w X_n dx / int w X_n^2 dx   (the eigenfunctions are orthogonal with the weight w = rho cp)

The roots are bracketed on a uniform lam grid fine enough to separate them (step = pi / (200 (A/sqrt(alpha1) +
C/sqrt(alpha2))), against a smallest spacing of about pi / (A/sqrt(alpha1) + C/sqrt(alpha2))) and refined by
bisection to the last bit. Terms are summed until exp(-lam^2 t) < 1e-18 (the remainder is below 1e-15 relative for
the sampled times; checked by the self-test below). Python 3.6, no numpy.

Self-test (python3 two_layer_exact.py): with two equal layers it reproduces slab_series.exact (the one-material
series) to round-off, and the series at t -> 0+ returns T0 inside the slab.
"""
import math


class TwoLayer:
    def __init__(self, k1, w1, k2, w2, A, L, T0, TW):
        self.k1, self.w1, self.k2, self.w2, self.A, self.L, self.T0, self.TW = k1, w1, k2, w2, A, L, T0, TW
        self.s1 = math.sqrt(k1 / w1)          # sqrt(alpha1)
        self.s2 = math.sqrt(k2 / w2)
        self.C = L - A
        self.roots = []
        self._lam_done = 0.0
        self._step = math.pi / (200.0 * (A / self.s1 + self.C / self.s2))

    def f(self, lam):
        b1, b2 = lam / self.s1, lam / self.s2
        return (self.k1 * b1 * math.sin(b1 * self.A) * math.sin(b2 * self.C)
                - self.k2 * b2 * math.cos(b1 * self.A) * math.cos(b2 * self.C))

    def _extend(self, lam_max):
        lo = self._lam_done if self._lam_done > 0.0 else 0.5 * self._step
        flo = self.f(lo)
        while lo < lam_max:
            hi = lo + self._step
            fhi = self.f(hi)
            if flo == 0.0:
                self.roots.append(lo)
            elif flo * fhi < 0.0:
                a, b, fa = lo, hi, flo
                for _ in range(200):
                    m = 0.5 * (a + b)
                    if m == a or m == b:
                        break
                    fm = self.f(m)
                    if fa * fm <= 0.0:
                        b = m
                    else:
                        a, fa = m, fm
                self.roots.append(0.5 * (a + b))
            lo, flo = hi, fhi
        self._lam_done = lo

    def term(self, lam):
        """(c_n, B) of the eigenvalue lam"""
        b1, b2 = lam / self.s1, lam / self.s2
        sa, ca = math.sin(b1 * self.A), math.cos(b1 * self.A)
        sc, cc = math.sin(b2 * self.C), math.cos(b2 * self.C)
        if abs(sc) > abs(cc):
            B = ca / sc
        else:
            B = self.k1 * b1 * sa / (self.k2 * b2 * cc)
        num = self.w1 * sa / b1 + self.w2 * B * (1.0 - cc) / b2
        den = (self.w1 * (0.5 * self.A + math.sin(2.0 * b1 * self.A) / (4.0 * b1))
               + self.w2 * B * B * (0.5 * self.C - math.sin(2.0 * b2 * self.C) / (4.0 * b2)))
        return num / den, B

    def T(self, x, t):
        if t <= 0.0:
            return self.T0
        lam_max = math.sqrt(18.0 * math.log(10.0) / t)       # exp(-lam^2 t) < 1e-18 beyond
        if self._lam_done < lam_max:
            self._extend(lam_max)
        s = 0.0
        for lam in self.roots:
            if lam > lam_max:
                break
            c, B = self.term(lam)
            if x <= self.A:
                X = math.cos(lam / self.s1 * x)
            else:
                X = B * math.sin(lam / self.s2 * (self.L - x))
            s += c * X * math.exp(-lam * lam * t)
        return self.TW + (self.T0 - self.TW) * s


def from_case(cfg):
    """TwoLayer of a make_two_layer.py case.json"""
    m1, m2 = cfg["materials"]
    return TwoLayer(m1["k"], m1["rho"] * m1["cp"], m2["k"], m2["rho"] * m2["cp"], cfg["A"], cfg["L"],
                    cfg["T0"], cfg["TW"])


if __name__ == "__main__":
    import os
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import slab_series
    a = 0.4 / 2.8e5
    ex = TwoLayer(0.4, 2.8e5, 0.4, 2.8e5, 0.025, 0.05, 298.0, 1644.0)
    worst = 0.0
    for t in (0.5, 5.0, 60.0):
        for x in (0.0, 0.01, 0.025, 0.04, 0.049):
            worst = max(worst, abs(ex.T(x, t) - slab_series.exact(0.05 - x, t, a)))
    print("equal layers vs slab_series: max |dT| = %.3e K" % worst)
    ex2 = TwoLayer(8.0, 1.12e7, 0.4, 2.8e5, 0.005, 0.01, 300.0, 1000.0)
    print("t = 1e-4 s, x = 0.002 / 0.007: %.6f %.6f K (T0 = 300)" % (ex2.T(0.002, 1e-4), ex2.T(0.007, 1e-4)))
    eps = 1e-9
    for t in (0.5, 4.0):
        jump = ex2.T(0.005 - eps, t) - ex2.T(0.005 + eps, t)
        q1 = -8.0 * (ex2.T(0.005 - eps, t) - ex2.T(0.005 - 2 * eps - 1e-7, t)) / (eps + 1e-7)
        q2 = -0.4 * (ex2.T(0.005 + 2 * eps + 1e-7, t) - ex2.T(0.005 + eps, t)) / (eps + 1e-7)
        print("t = %g s: T jump at A %.2e K, flux left/right of A %.6e %.6e W/m2" % (t, jump, q1, q2))
