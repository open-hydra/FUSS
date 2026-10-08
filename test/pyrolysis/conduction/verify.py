#!/usr/bin/env python3
"""Gate 4.8 (verification record 2026-10-08): pure conduction on the pyrolysis path against the exact slab series.

ATC1's geometry and graded mesh (100 cells, last/first 0.1, L = 0.05 m) with constant-property tables
(cp 1000 J/kg/K, k 0.4 W/m/K, rho_v 280 kg/m3 -> alpha = k / (rho cp)), no reactions, the heated face 2 stepped
to 1644 K at t = 0 (BC 302 time file `0 1644 / 60 1644`), face 1 insulated, T0 = 298 K, 60 s. Every probe is a
cell value compared with tools/slab_series.py at that cell's centre (no interpolation), at the probe times.

Checks (measured on 2026-10-08: 1.10 K at 1 mm / 0.5 s, 0.50 K at 1 mm / 1 s, 0.49 K at 2 mm / 1 s, < 0.14 K at >= 4 mm):
  (a) max |T - exact| over every probe and t in {0.5, 1, 2, 4, 8, 16, 30} s <= TOL_ALL   (1.3 K)
  (b) |T - exact| at the 1 mm pair, t = 1 s                                 <= TOL_1MM_1S (0.7 K)
  (c) max |T - exact| over the pairs at >= 2 mm and every time               <= TOL_DEEP  (0.6 K; measured 0.49 K at 2 mm, 1 s)
  The uniform-mesh companions (not gated) show second order: 14.8 K (100 cells) -> 3.8 K (200) at 1 mm / 0.5 s.

Usage: verify.py [case_dir]   Exit 0 on PASS, 1 on FAIL. Python 3.6, no numpy (imports ../tools/slab_series.py).
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tools"))
from slab_series import exact  # noqa: E402

TOL_ALL, TOL_1MM_1S, TOL_DEEP = 1.3, 0.7, 0.6
TIMES = [0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 30.0]
RHO, CP, K = 280.0, 1000.0, 0.4


def read_series(path):
    out = []
    for line in open(path):
        p = line.split()
        if len(p) >= 2:
            try:
                out.append((float(p[0]), float(p[1])))
            except ValueError:
                pass
    return out


def at(series, t):
    lo, hi = 0, len(series) - 1
    if t <= series[0][0]:
        return series[0][1]
    if t >= series[-1][0]:
        return series[-1][1]
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if series[mid][0] <= t:
            lo = mid
        else:
            hi = mid
    t0, v0 = series[lo]
    t1, v1 = series[hi]
    return v0 + (t - t0) * (v1 - v0) / (t1 - t0)


def main():
    case = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "graded100")
    pj = json.load(open(os.path.join(case, "probes.json")))
    L = pj["L"]
    alpha = K / (RHO * CP)
    Tw, T0 = 1644.0, 298.0
    ok = True
    worst_all = worst_deep = 0.0
    e1mm1s = 0.0
    print("conduction: %d cells, alpha = %.4e m2/s; FUSS cell value minus exact at the cell centre [K]" % (pj["n_cells"], alpha))
    print("  probe     depth[mm] " + " ".join("t=%-5g" % t for t in TIMES))
    for p in pj["probes"]:
        if p["role"] == "back":
            continue
        ser = read_series(os.path.join(case, "OUTPUT", p["name"] + ".txt"))
        depth = L - p["x"]
        row = []
        for t in TIMES:
            tt = min(ser, key=lambda s: abs(s[0] - t))[0]
            e = at(ser, tt) - exact(depth, tt, alpha, L, Tw, T0)
            row.append(e)
            worst_all = max(worst_all, abs(e))
            if p["depth_mm"] >= 2:
                worst_deep = max(worst_deep, abs(e))
            if p["depth_mm"] == 1 and t == 1.0:
                e1mm1s = max(e1mm1s, abs(e))
        print("  %-8s %8.4f  " % (p["name"], depth * 1e3) + " ".join("%7.3f" % v for v in row))
    for label, val, tol in (("max |T - exact|, all probes/times", worst_all, TOL_ALL), ("|T - exact| at 1 mm, t = 1 s", e1mm1s, TOL_1MM_1S), ("max |T - exact| at >= 2 mm", worst_deep, TOL_DEEP)):
        print("  %-36s = %.3f K (tol %g)  %s" % (label, val, tol, "PASS" if val <= tol else "FAIL"))
        ok &= val <= tol
    print("conduction: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
