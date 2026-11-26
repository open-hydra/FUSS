#!/usr/bin/env python3
"""Gate 4.8 (verification record 2026-10-08): pure conduction on the pyrolysis path against the exact slab series.

ATC1's geometry and graded mesh (N cells, last/first 0.1, L = 0.05 m; the committed case graded100 has 100) with
constant-property tables (cp 1000 J/kg/K, k 0.4 W/m/K, rho_v 280 kg/m3 -> alpha = k / (rho cp)), no reactions, the
heated face 2 stepped to 1644 K at t = 0 (BC 302 time file `0 1644 / 60 1644`), face 1 insulated, T0 = 298 K, 60 s.
The exact solution is tools/slab_series.py, at the probe sample nearest each listed time.

Three quantities per run, the error taken at the exact depths 1, 2, 4, 8, 16 mm (the same depths on every rung): the
error T - exact of each of the two cell centres that bracket the depth, interpolated linearly to it. The error, not T,
is interpolated because the exact solution exists at the cell centres, where FUSS has values; interpolating T would add
the exact profile's own interpolation error, which depends on where the centres fall on each rung. (a) max |error| over
every depth and t in {0.5, 1, 2, 4, 8, 16, 30} s; (b) |error| at 1 mm, t = 1 s; (c) max |error| at >= 2 mm and every
time. Check (plan 17 P17-D1, the
convergence rule of tools/conv_rule.py): on a ladder of runs of make_case.py --kind conduction (coarsest first; the gate
runs 50, 100 and 200 cells, the same 0.1 grading) each quantity falls on every rung and its order on the finest pair
lies in [0.8 p, 1.2 p], p = 2 (second order in space: test/pyrolysis/README.md; the uniform-mesh companion went
14.8 K (100 cells) -> 3.8 K (200) at 1 mm / 0.5 s; the step size, RK3 at fixed vnn, scales as dx^2).

Usage: verify.py case_dir [case_dir ...]   (one dir: the table only, exit 1). Exit 0 on PASS, 1 on FAIL.
FUSS_GATE_SELFTEST=1 empties the order window. Python 3.6, no numpy (imports ../tools/slab_series.py, conv_rule.py).
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tools"))
from slab_series import exact  # noqa: E402
import conv_rule as cr  # noqa: E402

P_EXPECTED = 2
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


def errors(case):
    """(n_cells, (a), (b), (c)) of one run, printing its table"""
    pj = json.load(open(os.path.join(case, "probes.json")))
    L = pj["L"]
    alpha = K / (RHO * CP)
    Tw, T0 = 1644.0, 298.0
    worst_all = worst_deep = 0.0
    e1mm1s = 0.0
    print("conduction: %d cells, alpha = %.4e m2/s; T - exact at the two bracketing cell centres, interpolated to the "
          "exact depth [K]" % (pj["n_cells"], alpha))
    print("  depth[mm] " + " ".join("t=%-5g" % t for t in TIMES))
    px = {p["name"]: p["x"] for p in pj["probes"]}
    for dmm, (lo, hi) in sorted(pj["pairs"].items(), key=lambda kv: int(kv[0])):
        slo = read_series(os.path.join(case, "OUTPUT", lo + ".txt"))
        shi = read_series(os.path.join(case, "OUTPUT", hi + ".txt"))
        depth = int(dmm) * 1.0e-3
        w = ((L - depth) - px[lo]) / (px[hi] - px[lo])
        if not -1.0e-9 <= w <= 1.0 + 1.0e-9:
            raise SystemExit("conduction: FAIL  probes %s/%s do not bracket %s mm (w = %g)" % (lo, hi, dmm, w))
        row = []
        for t in TIMES:
            tt = min(slo, key=lambda s: abs(s[0] - t))[0]
            e = ((1.0 - w) * (at(slo, tt) - exact(L - px[lo], tt, alpha, L, Tw, T0))
                 + w * (at(shi, tt) - exact(L - px[hi], tt, alpha, L, Tw, T0)))
            row.append(e)
            worst_all = max(worst_all, abs(e))
            if int(dmm) >= 2:
                worst_deep = max(worst_deep, abs(e))
            if int(dmm) == 1 and t == 1.0:
                e1mm1s = max(e1mm1s, abs(e))
        print("  %8d   " % int(dmm) + " ".join("%7.3f" % v for v in row))
    return pj["n_cells"], worst_all, e1mm1s, worst_deep


def main():
    cases = sys.argv[1:]
    if not cases:
        sys.exit(__doc__)
    runs = [errors(c) for c in cases]
    if len(runs) < 2:
        print("conduction: FAIL  the convergence rule needs a ladder of 2 or more runs")
        return 1
    ok = True
    ns = [r[0] for r in runs]
    for k, label in ((1, "(a) max |T - exact|, all depths/times [K]"), (2, "(b) |T - exact| at 1 mm, t = 1 s [K]"),
                     (3, "(c) max |T - exact| at >= 2 mm [K]")):
        good, text = cr.check(label, [r[k] for r in runs], ns, P_EXPECTED, fmt="%.4f")
        ok &= good
        print("  %s  %s" % (text, "PASS" if good else "FAIL"))
    print("conduction: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
