#!/usr/bin/env python3
"""Gate 4.4 of plan 07: Ablation Test Case 2 chemistry-off under the type-2 (Darcy) closure against
the same local PATO run as plan 06's gate 4.5, now with the pressure probed too.

The case is test/pyrolysis/ATC2-chemistryOff with [FUSS-Pyrolysis] gas-closure = type2 and
[FUSS-Darcy] p-init = p-wall = 101325, solver tridiag, TACOT's permeability (K_v 1.6e-11, K_c
2e-11 m2); the heated 303 face is open, the back face impermeable. PATO solves the same type-2
model (DarcyLawMassModel + PyrolysisEnergyModel). The oracle run of 2026-10-08 added `p` to the
probed fields of plotDict and is otherwise the 2026-10-07 run byte for byte (ref/Ta_plot identical):
ref/Ta_plot and ref/p_plot, rows `t v1..v7` every 0.1 s at 1, 2, 4, 8, 12, 16, 24 mm below the
heated face. $FUSS_PATO_ATC2_OUT (a PATO output/porousMat/scalar directory) overrides ref/.

FUSS samples each depth with two bracketing index-position probes (probes.json; columns t T p_pg),
linear in x then linear in time onto PATO's 0.1 s grid.

Checks (expected.json):
  (T)  max_{t = 0.1..60 s, 7 depths} |T_FUSS - T_PATO| <= budget_K, the regression budget = this case's measured
       maximum + 2 K (plan 06 section 4.5's rule, the same for both closures; 10.52 K measured 2026-10-13, after the open
       face took the gas enthalpy at the wall temperature; 15.79 K before)
  (acc) the same maximum <= acceptance_K = B_T = the measured type-1 maximum, 22.75 K at 1 mm / 2.6 s (plan 07 D4 and
       section 5: removing the type-1 approximation must not make the answer worse). Resolution of 2026-10-13: the
       24.75 K of the 2026-10-08 decision was the type-1 REGRESSION budget, not this acceptance; with 22.75 K the
       derivation of the 0.15 below holds (14.6 %; the plan's 14.9 % left dmu/mu unscaled, and 24.75 K would give 15.4 %)
  (dir) at the probe and time where the type-1 run's |T_FUSS - T_PATO| peaked (1 mm, 2.6 s), the
       type-2 value is smaller than that peak
  (p)  max_{t,probe} |dp_FUSS - dp_PATO| <= dp_tol_rel * max_{t,probe} dp_PATO, dp = p - p_w (0.15: the
       derivation of plan 07 section 4.4)
  (sign, with --type1 DIR = the type-1 run of test/pyrolysis/ATC2-chemistryOff on the same mesh)
       T_type2 - T_type1 <= SIGN_SLACK wherever T_PATO >= T_HOT (700 K): the gas storage makes type 2
       cooler at every hot probe (the exact 0-D shift is -1 K at 700 K, -2.2 K at 1000 K, -3 K at
       1500 K); the magnitude is printed, not gated

Usage: verify.py [case_dir] [--type1 DIR] [--record]   Exit 0 on PASS, 1 on FAIL. Python 3.6, no numpy.
FUSS_GATE_SELFTEST=1 sets every tolerance to 0 (the proven-red hook).
"""
import json
import os
import sys

DEPTHS_MM = [1, 2, 4, 8, 12, 16, 24]
P_WALL = 101325.0
T_HOT = 700.0
SIGN_SLACK = 0.05      # K: the -eps p gas-work term heats a cell by ~eps dp/(rho c) ~ 6e-3 K at dp = 2 kPa


def read_series(path, col):
    out = []
    for line in open(path):
        p = line.split()
        if len(p) > col:
            try:
                out.append((float(p[0]), float(p[col])))
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
    return {v["name"]: v["x"] for v in pj["probes"]}


def bracket_series(case, pj, px, L, col):
    series = {}
    for d in DEPTHS_MM:
        lo, hi = pj["pairs"][str(d)]
        s_lo = read_series(os.path.join(case, "OUTPUT", lo + ".txt"), col)
        s_hi = read_series(os.path.join(case, "OUTPUT", hi + ".txt"), col)
        x_t = L - d * 1.0e-3
        w = (x_t - px[lo]) / (px[hi] - px[lo])
        if not (-1.0e-9 <= w <= 1.0 + 1.0e-9):
            raise SystemExit("FAIL: probes %s/%s do not bracket x = %g (w = %g)" % (lo, hi, x_t, w))
        series[d] = [(t0, (1.0 - w) * v0 + w * v1) for (t0, v0), (t1, v1) in zip(s_lo, s_hi)]
    return series


def main():
    argv = sys.argv[1:]
    record = "--record" in argv
    type1 = None
    if "--type1" in argv:
        i = argv.index("--type1")
        type1 = os.path.abspath(argv[i + 1])
        argv = argv[:i] + argv[i + 2:]
    args = [a for a in argv if a != "--record"]
    case = os.path.abspath(args[0]) if args else os.path.dirname(os.path.abspath(__file__))
    here = os.path.dirname(os.path.abspath(__file__))
    pj = json.load(open(os.path.join(case, "probes.json")))
    L = pj["L"]
    px = probe_x(pj)
    oracle_dir = os.environ.get("FUSS_PATO_ATC2_OUT") or os.path.join(here, "ref")
    ta_plot = os.path.join(oracle_dir, "Ta_plot")
    p_plot = os.path.join(oracle_dir, "p_plot")
    for f in (ta_plot, p_plot):
        if not os.path.isfile(f):
            print("SKIP: no %s (set FUSS_PATO_ATC2_OUT to a PATO output/porousMat/scalar directory probing Ta and p)" % f)
            return 0
    pato_T = read_pato_plot(ta_plot)
    pato_p = read_pato_plot(p_plot)
    expected_path = os.path.join(here, "expected.json")
    expected = json.load(open(expected_path)) if os.path.isfile(expected_path) else {}
    budget = expected.get("budget_K", 12.52)
    accept = expected.get("acceptance_K", 22.7483856120341)
    peak1 = expected.get("type1_peak_K", 22.7483856120341)
    peak_d, peak_t = expected.get("type1_peak_depth_mm", 1), expected.get("type1_peak_t", 2.6)
    dp_tol = expected.get("dp_tol_rel", 0.15)
    sign_slack = SIGN_SLACK
    if os.environ.get("FUSS_GATE_SELFTEST", "0") == "1":
        budget, accept, peak1, dp_tol, sign_slack = 0.0, 0.0, 0.0, 0.0, -1.0

    sT = bracket_series(case, pj, px, L, 1)
    sp = bracket_series(case, pj, px, L, 2)
    s1 = bracket_series(type1, pj, px, L, 1) if type1 else None
    ok = True
    print("ATC2-chemistryOff-type2 gate: FUSS type 2 (bracketing probes) vs PATO (%s), t = 0.1..60 s" % oracle_dir)

    # (T)
    worst_all, worst_at = 0.0, None
    dT_at_peak1 = None
    for col, d in enumerate(DEPTHS_MM, 1):
        worst, t_worst, n = 0.0, None, 0
        for row in pato_T:
            t = row[0]
            if t < 0.1 - 1e-9 or t > 60.0 + 1e-9:
                continue
            v = interp_t(sT[d], t)
            if v is None:
                continue
            n += 1
            dT = abs(v - row[col])
            if dT > worst:
                worst, t_worst = dT, t
            if d == peak_d and abs(t - peak_t) < 1e-6:
                dT_at_peak1 = dT
        if n == 0:
            print("FAIL: no FUSS sample inside PATO's window for depth %d mm" % d)
            return 1
        if worst > worst_all:
            worst_all, worst_at = worst, (d, t_worst)
        print("  %2d mm: max |dT| = %6.2f K at t = %5.1f s  (%d samples)" % (d, worst, t_worst, n))
    good = worst_all <= budget
    ok &= good
    print("  (T)   max over the seven depths: %.2f K at %s  (budget %.2f K)  %s" % (worst_all, worst_at, budget, "PASS" if good else "FAIL"))
    good = worst_all <= accept
    ok &= good
    print("  (acc) the same maximum against the measured type-1 maximum %.2f K  %s" % (accept, "PASS" if good else "FAIL"))
    if dT_at_peak1 is None:
        print("  (dir) FAIL: no sample at the type-1 peak (%g mm, %g s)" % (peak_d, peak_t))
        ok = False
    else:
        good = dT_at_peak1 < peak1
        ok &= good
        print("  (dir) |dT| at the type-1 peak (%g mm, %g s): %.2f K vs type 1's %.2f K  %s"
              % (peak_d, peak_t, dT_at_peak1, peak1, "PASS" if good else "FAIL"))

    # (p)
    worst_dp, dp_max, dp_at = 0.0, 0.0, None
    fuss_dp_max = 0.0
    for col, d in enumerate(DEPTHS_MM, 1):
        for row in pato_p:
            t = row[0]
            if t < 0.1 - 1e-9 or t > 60.0 + 1e-9:
                continue
            v = interp_t(sp[d], t)
            if v is None:
                continue
            dpp = row[col] - P_WALL
            dpf = v - P_WALL
            dp_max = max(dp_max, dpp)
            fuss_dp_max = max(fuss_dp_max, dpf)
            if abs(dpf - dpp) > worst_dp:
                worst_dp, dp_at = abs(dpf - dpp), (d, t, dpf, dpp)
    good = worst_dp <= dp_tol * dp_max
    ok &= good
    print("  (p)   max |dp_FUSS - dp_PATO| = %.1f Pa at %d mm, t = %.1f s (FUSS %.1f, PATO %.1f); PATO peak dp %.1f Pa, FUSS peak %.1f Pa"
          % (worst_dp, dp_at[0], dp_at[1], dp_at[2], dp_at[3], dp_max, fuss_dp_max))
    print("        ratio %.3f  (tol %.2f)  %s" % (worst_dp / dp_max if dp_max > 0 else float("inf"), dp_tol, "PASS" if good else "FAIL"))

    # (sign) against the type-1 run on the same mesh
    shifts = {}
    if s1 is not None:
        worst_sign = -1e30
        for col, d in enumerate(DEPTHS_MM, 1):
            dmax, dmin = -1e30, 1e30
            for row in pato_T:
                t = row[0]
                if t < 0.1 - 1e-9 or t > 60.0 + 1e-9 or row[col] < T_HOT:
                    continue
                v2, v1 = interp_t(sT[d], t), interp_t(s1[d], t)
                if v2 is None or v1 is None:
                    continue
                dmax, dmin = max(dmax, v2 - v1), min(dmin, v2 - v1)
            if dmax > -1e30:
                shifts[str(d)] = [dmin, dmax]
                worst_sign = max(worst_sign, dmax)
                print("  %2d mm: T_type2 - T_type1 over T_PATO >= %g K: %.2f .. %.2f K" % (d, T_HOT, dmin, dmax))
        good = worst_sign <= sign_slack
        ok &= good
        print("  (sign) max (T_type2 - T_type1) at the hot samples: %.3f K  (slack %g K)  %s" % (worst_sign, sign_slack, "PASS" if good else "FAIL"))
    else:
        print("  (sign) SKIP: no --type1 run given")

    if record:
        expected.update({"measured_K": worst_all, "measured_at": list(worst_at), "measured_dp_rel": worst_dp / dp_max,
                         "measured_dT_at_type1_peak_K": dT_at_peak1, "pato_dp_peak_Pa": dp_max, "fuss_dp_peak_Pa": fuss_dp_max,
                         "type2_minus_type1_K": shifts})
        json.dump(expected, open(expected_path, "w"), indent=1, sort_keys=True)
        print("  recorded %s" % expected_path)
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
