#!/usr/bin/env python3
"""Gate 4.5 of plan 06: Ablation Test Case 2 chemistry-off against a local PATO run.

The same slab, mesh and tables as ATC1 with T(0) = 300 K and the heated face
under BC 303 `hconv = 300 W/m2K, qrad = 0, Tref = 2800 K` for 60 s. The oracle
is PATO 3.1's tutorials/1D/AblationTestCase_2.x_chemistryOff run with the
emissivity column of its TACOT virgin/char tables set to 0, every row of
constant/porousMat/BoundaryConditions at hconv = 300, Tedge = 2800,
chemistryOn = 0 (ramp and cool-down removed) and endTime 60 -- edit_pato_case.py
applies the edits -- so that both codes evaluate the same closed-form surface
energy balance hconv (Tedge - T_w) with no re-radiation. PATO's probe output
(output/porousMat/scalar/Ta_plot: header line, then rows `t v1..v7` every
0.1 s, v1..v7 at 1, 2, 4, 8, 12, 16, 24 mm below the heated face, interpolated
by PATO at the points) is read from $FUSS_PATO_ATC2_OUT when set, else from the
bundled ref/ copy of the 2026-10-07 run.

FUSS samples each depth with two bracketing index-position probes (probes.json),
linear in x then linear in time onto PATO's 0.1 s grid.

Pass: max_{t = 0.1..60 s, 7 depths} |T_FUSS - T_PATO| <= budget_K (expected.json,
recorded on the first run as measured + 2 K with --record; provisional 20 K).

Usage: verify.py [case_dir] [--record]   Exit 0 on PASS, 1 on FAIL. Python 3.6, no numpy.
"""
import json, os, sys

DEPTHS_MM = [1, 2, 4, 8, 12, 16, 24]
BUDGET_DEFAULT = 20.0


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


def interp_t(series, t):
    if not series or t < series[0][0] or t > series[-1][0]:
        return None
    lo, hi = 0, len(series) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if series[mid][0] <= t:
            lo = mid
        else:
            hi = mid
    t0, v0 = series[lo]
    t1, v1 = series[hi]
    return v0 if t1 == t0 else v0 + (v1 - v0) * (t - t0) / (t1 - t0)


def read_pato_plot(path):
    rows = []
    for line in open(path):
        s = line.strip()
        if not s or s.startswith("//") or s.startswith("#"):
            continue
        try:
            rows.append([float(x) for x in s.split()])
        except ValueError:
            pass
    return rows


def probe_x(pj):
    out = {}
    probes = pj["probes"]
    if isinstance(probes, dict):
        for name, v in probes.items():
            out[name] = v["x"] if isinstance(v, dict) else v
    else:
        for v in probes:
            out[v["name"]] = v["x"]
    return out


def main():
    args = [a for a in sys.argv[1:] if a != "--record"]
    record = "--record" in sys.argv[1:]
    case = os.path.abspath(args[0]) if args else os.path.dirname(os.path.abspath(__file__))
    here = os.path.dirname(os.path.abspath(__file__))
    pj = json.load(open(os.path.join(case, "probes.json")))
    L = pj["L"]
    px = probe_x(pj)
    oracle_dir = os.environ.get("FUSS_PATO_ATC2_OUT") or os.path.join(here, "ref")
    plot = os.path.join(oracle_dir, "Ta_plot")
    if not os.path.isfile(plot):
        plot = os.path.join(oracle_dir, "scalar", "Ta_plot")
    if not os.path.isfile(plot):
        print("SKIP: no PATO Ta_plot under %s (set FUSS_PATO_ATC2_OUT to the PATO output/porousMat directory)" % oracle_dir)
        return 0
    pato = read_pato_plot(plot)
    expected_path = os.path.join(here, "expected.json")
    expected = json.load(open(expected_path)) if os.path.isfile(expected_path) else {}
    budget = expected.get("budget_K", BUDGET_DEFAULT)

    series = {}
    for d in DEPTHS_MM:
        lo, hi = pj["pairs"][str(d)]
        s_lo = read_series(os.path.join(case, "OUTPUT", lo + ".txt"))
        s_hi = read_series(os.path.join(case, "OUTPUT", hi + ".txt"))
        x_t = L - d * 1.0e-3
        w = (x_t - px[lo]) / (px[hi] - px[lo])
        if not (-1.0e-9 <= w <= 1.0 + 1.0e-9):
            print("FAIL: probes %s/%s do not bracket x = %g (w = %g)" % (lo, hi, x_t, w)); return 1
        series[d] = [(t0, (1.0 - w) * v0 + w * v1) for (t0, v0), (t1, v1) in zip(s_lo, s_hi)]

    print("ATC2-chemistryOff gate: FUSS (bracketing probes) vs PATO (%s), t = 0.1..60 s" % plot)
    worst_all = 0.0
    for col, d in enumerate(DEPTHS_MM, 1):
        worst, t_worst, n = 0.0, None, 0
        for row in pato:
            t = row[0]
            if t < 0.1 - 1e-9 or t > 60.0 + 1e-9:
                continue
            v = interp_t(series[d], t)
            if v is None:
                continue
            n += 1
            dT = abs(v - row[col])
            if dT > worst:
                worst, t_worst = dT, t
        if n == 0:
            print("FAIL: no FUSS sample inside PATO's window for depth %d mm" % d); return 1
        worst_all = max(worst_all, worst)
        print("  %2d mm: max |dT| = %6.2f K at t = %5.1f s  (%d samples)" % (d, worst, t_worst, n))
    if record:
        expected["measured_K"] = worst_all
        expected["budget_K"] = worst_all + 2.0
        expected["oracle"] = plot
        json.dump(expected, open(expected_path, "w"), indent=1, sort_keys=True)
        budget = worst_all + 2.0
        print("  recorded %s" % expected_path)
    print("  max over the seven depths: %.2f K  (budget %.1f K)" % (worst_all, budget))
    ok = worst_all <= budget
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
