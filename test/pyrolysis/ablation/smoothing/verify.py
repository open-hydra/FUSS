#!/usr/bin/env python3
"""
Plan 08 section 3.6: the wall smoothing (smoothing-width, thesis Eq. 5.1) damps the period-4 node mode of the recession
by its transfer factor, 5 for the 5-node stencil.

    python3 verify.py W1_DIR W5_DIR [W5_2BLK_DIR]

From each run's OUTPUT/field.tec (the last output, the moved mesh): the recession of the top surface nodes,
s_k = y_k(0) - y_k(t) at k = 0..40 (y_k(0) = 0.05 m, j = 20, the k = 0 layer), fitted at the interior nodes 3..37 (the
truncated stencil excluded) by least squares to m + a cos(pi k/2) + b sin(pi k/2); the mode's amplitude is
sqrt(a^2 + b^2). Checks:
  (mode)   the width-1 run has the mode: amplitude > 10 % of the 0.1 * mean recession the pattern imposes
  (ratio)  amplitude(width 1) / amplitude(width 5) >= 4.5 (theory 5.0: (1 + 2 cos(pi/2) + 2 cos(pi))/5 = -1/5;
           the margin covers the nonlinear SEB response and the geometry feedback over 20 s)
  (sign)   the width-5 mode is reversed: a(w5) a(w1) + b(w5) b(w1) < 0 (the transfer is negative)
  (seam)   with W5_2BLK_DIR (w5 cut into two blocks joined by 101 faces at x = 0.1 m): the windows of the nodes beside
           the seam cross it, so the surface equals w5's to 1e-10 m at every node, and T to 1e-7 of its maximum (the
           connection kernels, plan 08 section 2.9). The two runs solve the pressure along different sweep paths (the
           line-uniform correction spans one block) and stop within the Darcy tol (1e-8) apart: measured 2026-10-13
           7.9e-12 m and 4.0e-9. The old connection rows give 5.6e-8 m and 1.1e-4 here; a window truncated at the
           seam would move the nodes beside it by ~1e-5 m
Exit 0 on PASS, 1 on FAIL. FUSS_GATE_SELFTEST=1 asks a ratio of 1e9.
"""
import math
import os
import re
import sys

LY, NY, NX = 0.05, 20, 40
SELFTEST = os.environ.get("FUSS_GATE_SELFTEST") == "1"
RATIO_MIN = 1.0e9 if SELFTEST else 4.5
TOL_SEAM_S, TOL_SEAM_T = (0.0, 0.0) if SELFTEST else (1.0e-10, 1.0e-7)


def zones(run):
    """[(I, J, y nodes, T cells)] of every zone of OUTPUT/field.tec."""
    txt = open(os.path.join(run, "OUTPUT", "field.tec")).read()
    names = txt.split("\n", 1)[0].split("=", 1)[1].split()
    out = []
    for part in re.split(r"\n\s*ZONE", txt)[1:]:
        I, J, K = (int(v) for v in re.search(r"I=(\d+), J=(\d+), K=(\d+)", part.split("\n", 1)[0]).groups())
        body = [float(v) for v in part.split("\n", 1)[1].split()]
        nn, nc = I * J * K, (I - 1) * (J - 1) * (K - 1)
        it = names.index("T")
        out.append((I, J, body[nn:2 * nn], body[3 * nn + (it - 3) * nc:3 * nn + (it - 2) * nc]))
    return out


def top_recession(run):
    """s_k at the top nodes k = 0..40 (the k = 0 layer, index (k J + j) I + i), blocks side by side along x"""
    s = []
    for n, (I, J, y, T) in enumerate(zones(run)):
        row = [LY - y[(J - 1) * I + i] for i in range(I)]
        s.extend(row if n == 0 else row[1:])                # a seam node once
    assert len(s) == NX + 1
    return s


def cell_T(run):
    """T of the cells in the global order (j, i), blocks side by side along x"""
    zs = zones(run)
    nxb = [I - 1 for (I, J, y, T) in zs]
    out = []
    for j in range(NY):
        for (I, J, y, T), nx in zip(zs, nxb):
            out.extend(T[j * nx:(j + 1) * nx])
    return out


def fit(s):
    ks = list(range(3, NX - 2))                              # 3..37
    rows = [(1.0, math.cos(math.pi * k / 2), math.sin(math.pi * k / 2)) for k in ks]
    # normal equations (3 x 3)
    A = [[sum(r[a] * r[b] for r in rows) for b in range(3)] for a in range(3)]
    y = [sum(r[a] * s[k] for r, k in zip(rows, ks)) for a in range(3)]
    # Gaussian elimination
    M = [A[i][:] + [y[i]] for i in range(3)]
    for c in range(3):
        p = max(range(c, 3), key=lambda r: abs(M[r][c]))
        M[c], M[p] = M[p], M[c]
        for r in range(3):
            if r != c:
                f = M[r][c] / M[c][c]
                M[r] = [x - f * z for x, z in zip(M[r], M[c])]
    return [M[i][3] / M[i][i] for i in range(3)]


def main():
    w1, w5 = sys.argv[1], sys.argv[2]
    m1, a1, b1 = fit(top_recession(w1))
    m5, a5, b5 = fit(top_recession(w5))
    A1, A5 = math.hypot(a1, b1), math.hypot(a5, b5)
    ratio = A1 / A5 if A5 > 0 else float("inf")
    ok = True
    good = A1 > 0.1 * 0.1 * m1
    ok &= good
    print("(mode)   width 1: mean recession %.4e m, period-4 amplitude %.4e m (%.2f %% of the mean)  %s"
          % (m1, A1, 100 * A1 / m1, "PASS" if good else "FAIL"))
    good = ratio >= RATIO_MIN
    ok &= good
    print("(ratio)  amplitude w1 / w5 = %.4f (>= %g; theory 5)  width 5: mean %.4e m, amplitude %.4e m  %s"
          % (ratio, RATIO_MIN, m5, A5, "PASS" if good else "FAIL"))
    good = a1 * a5 + b1 * b5 < 0.0
    ok &= good
    print("(sign)   mode w1 (%.3e, %.3e), w5 (%.3e, %.3e): reversed  %s" % (a1, b1, a5, b5, "PASS" if good else "FAIL"))
    if len(sys.argv) > 3:
        s2, s5 = top_recession(sys.argv[3]), top_recession(w5)
        ds = max(abs(a - b) for a, b in zip(s2, s5))
        T2, T5 = cell_T(sys.argv[3]), cell_T(w5)
        dT = max(abs(a - b) for a, b in zip(T2, T5)) / max(abs(v) for v in T5)
        good = ds <= TOL_SEAM_S and dT <= TOL_SEAM_T
        ok &= good
        print("(seam)   two blocks vs one: surface max %.3e m (<= %.0e), T %.3e of max (<= %.0e)  %s"
              % (ds, TOL_SEAM_S, dT, TOL_SEAM_T, "PASS" if good else "FAIL"))
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
