#!/usr/bin/env python3
"""Gate (plan 10, D14): BC 303's qrad must reach BOTH the flux and the wall output.

The steady 1-D slab of generate.py has a closed form (T_w = 700 K, q_w = 4.0e4 W/m2,
T linear between the faces). Two checks, one per dispatch site of BC 303:

  (1) interior: max |T_i - T(x_i)| over every cell of OUTPUT/field.tec <= TOL_K.
      The flux site (Mod_BC_Fluxes) decides the interior: with it passing
      anything but qrad the interior is wrong by tens of K.
  (2) wall: on every cell of zone B1F2 of OUTPUT/wall.tec, |Tw - T_w| <= TOL_K and
      |qw - q_w| <= TOL_QREL q_w. The wall output (IO_Wall) recomputes Tw and qw
      with a second call of the BC routine, so this sees the output site.

Tolerance. The run stops on residuotot = sqrt(sum dT^2) <= res-threshold = 1e-10
K per iteration. The slowest error mode (Dirichlet-Robin, mu cot mu = -Bi, Bi =
hL/k = 1, mu = 2.029) decays per iteration by dt lambda = 0.5 vnn (mu/L)^2 /
(1/dx^2 + 1/dy^2) = 5.1e-3 (IRS leaves the slowest mode untouched), so the error
left is about 1e-10 / 5.1e-3 = 2e-8 K; round-off and the 15-digit output are
1e-12 K at 700 K. TOL_K = 1e-6 K keeps a 50x margin (the budget
conjugate_two_layer uses for the same reason). The wall flux: 1e-8 relative.

Usage: verify.py [case_dir]   (default: this directory). Exit 0 PASS, 1 FAIL.
Python 3.6, no numpy.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from generate import t_exact, T_W, Q_W   # noqa: E402

TOL_K = 1.0e-6
TOL_QREL = 1.0e-8


def read_zones(path):
    """{zone title: {var: [values]}} of a block-packed Tecplot file, every zone."""
    text = open(path).read().split("\n")
    names = []
    zones, cur, raw = {}, None, []

    def close():
        if cur is None:
            return
        title, ni, nj, nk, cc = cur
        n_node = ni * nj * nk
        n_cell = max(ni - 1, 1) * max(nj - 1, 1) * max(nk - 1, 1)
        out, idx = {}, 0
        for v, name in enumerate(names, 1):
            n = n_cell if cc[0] <= v <= cc[1] else n_node
            out[name] = raw[idx:idx + n]
            idx += n
        if idx != len(raw):
            sys.exit("FAIL: %s zone %s holds %d values, the header implies %d" % (path, title, len(raw), idx))
        zones[title] = out

    for line in text:
        s = line.strip()
        if re.match(r"(?i)variables", s):
            names = [t.strip('",') for t in s.split("=", 1)[1].split() if t.strip('",')]
            continue
        if re.match(r"(?i)zone", s):
            close()
            g = lambda p: re.search(p, s, re.I)
            m = g(r"\[(\d+)-(\d+)\]\s*=\s*CELLCENTERED")
            cc = tuple(map(int, m.groups())) if m else (0, -1)
            cur = (g(r"T\s*=\s*([\w-]+)").group(1), int(g(r"\bI\s*=\s*(\d+)").group(1)),
                   int(g(r"\bJ\s*=\s*(\d+)").group(1)), int(g(r"\bK\s*=\s*(\d+)").group(1)), cc)
            raw = []
            continue
        for tok in s.split():
            try:
                raw.append(float(tok))
            except ValueError:
                pass
    close()
    return zones


def main():
    case = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else HERE
    ok = True

    field = read_zones(os.path.join(case, "OUTPUT", "field.tec"))
    if len(field) != 1:
        print("FAIL: expected one block in field.tec, found %d" % len(field))
        return 1
    z = list(field.values())[0]
    xs = z["x"]
    T = z["T"]
    # cell centres from the nodes (nodal x is i-fastest, cells too)
    from generate import NX, NY
    err = 0.0
    for j in range(NY):
        for i in range(NX):
            xc = 0.5 * (xs[i] + xs[i + 1])
            e = abs(T[j * NX + i] - t_exact(xc))
            err = max(err, e)
    good = err <= TOL_K
    ok &= good
    print("  %s  (1) interior  max |T - T(x)| = %.3e K   (tol %.1e)" % ("PASS" if good else "FAIL", err, TOL_K))

    wall = read_zones(os.path.join(case, "OUTPUT", "wall.tec"))
    if "B1F2" not in wall:
        print("FAIL: wall.tec has no zone B1F2 (zones: %s)" % ", ".join(sorted(wall)))
        return 1
    w = wall["B1F2"]
    etw = max(abs(v - T_W) for v in w["Tw"])
    eqw = max(abs(v - Q_W) for v in w["qw"])
    good = etw <= TOL_K and eqw <= TOL_QREL * Q_W
    ok &= good
    print("  %s  (2) wall      max |Tw - %.1f| = %.3e K, max |qw - %.1f| = %.3e W/m2   (tol %.1e K, %.1e W/m2)"
          % ("PASS" if good else "FAIL", T_W, etw, Q_W, eqw, TOL_K, TOL_QREL * Q_W))

    print("bc303_qrad_1d: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
