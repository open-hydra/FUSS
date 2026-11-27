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
linear in x, then onto PATO's 0.1 s grid by the cubic through the four probe rows
around each PATO time (tools/abl_tools.interp4: a probe row lands up to one step
after its 0.1 s multiple, and linear in time that offset moves the maximum with dt).

Pass: max_{t = 0.1..60 s, 7 depths} |T_FUSS - T_PATO| <= the T spread of Ablation Test Case 2.1 (plan 17
P17-D4: one external code here, so the tolerance is the PATO - Amaryllis spread of ATC 2.1, same material and heating,
computed from test/pyrolysis/ablation/ATC21/ref by tools/abl_tools.spread_atc21).

Usage: verify.py [case_dir]   Exit 0 on PASS, 1 on FAIL. FUSS_GATE_SELFTEST=1 sets the tolerance to 0.
Python 3.6, no numpy.
"""
import json, os, sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools"))
import abl_tools  # noqa: E402

DEPTHS_MM = [1, 2, 4, 8, 12, 16, 24]


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
    return abl_tools.interp4(series, t)


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
    args = sys.argv[1:]
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
    s21 = abl_tools.spread_atc21()[0]["T"]
    budget = 0.0 if os.environ.get("FUSS_GATE_SELFTEST", "0") == "1" else s21[0]

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
    print("  max over the seven depths: %.2f K  (ATC 2.1 T spread %.2f K, at %g s %s)" % (worst_all, budget, s21[1], s21[2]))
    ok = worst_all <= budget
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
