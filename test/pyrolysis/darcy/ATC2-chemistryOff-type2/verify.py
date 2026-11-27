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
linear in x, then onto PATO's 0.1 s grid by the cubic through the four probe rows around each PATO time
(tools/abl_tools.interp4). A probe row lands up to one step after its 0.1 s multiple; linear in time, that
offset moved (T) by 0.05 K between the time steps vnn 0.25 / 0.125 / 0.0625, against 3e-4 K of the solution.

Checks (plan 17: one external code, so the T tolerance is the PATO - Amaryllis T spread of Ablation Test Case 2.1,
P17-D4, from test/pyrolysis/ablation/ATC21/ref by tools/abl_tools.spread_atc21):
  (T)  max_{t = 0.1..60 s, 7 depths} |T_FUSS - T_PATO| <= the ATC 2.1 T spread
  (acc, with --type1) the same maximum <= the type-1 run's own maximum, from the run given (plan 07 D4 and section 5:
       removing the type-1 approximation must not make the answer worse)
  (dir, with --type1) at the probe and time where that type-1 run's |T_FUSS - T_PATO| peaks, the type-2 value is
       smaller than that peak
  (p)  max_{t,probe} |dp_FUSS - dp_PATO| <= dp_tol_rel * max_{t,probe} dp_PATO, dp = p - p_w; dp_tol_rel = the T terms
       of plan 07 section 4.4 (delta T/T + delta mu/mu + delta W/W = 8.4 % at delta T = 20 K, 850 K, scaled to delta T =
       the ATC 2.1 T spread) + the ATC 2.1 mdot_g spread over PATO's mdot_g peak (P17-D13, option of the plan)
  (sign, with --type1 DIR = the type-1 run of test/pyrolysis/ATC2-chemistryOff on the same mesh)
       T_type2 - T_type1 <= SIGN_SLACK wherever T_PATO >= T_HOT (700 K): the gas storage makes type 2
       cooler at every hot probe (the exact 0-D shift is -1 K at 700 K, -2.2 K at 1000 K, -3 K at
       1500 K); the magnitude is printed, not gated

Usage: verify.py [case_dir] [--type1 DIR]   Exit 0 on PASS, 1 on FAIL. Python 3.6, no numpy.
FUSS_GATE_SELFTEST=1 sets every tolerance to 0 (the proven-red hook).
"""
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tools"))
import abl_tools  # noqa: E402

DEPTHS_MM = [1, 2, 4, 8, 12, 16, 24]
DP_T_TERMS_PER_K = (0.024 + 0.021 + 0.039) / 20.0   # plan 07 section 4.4: dT/T + dmu/mu + dW/W at dT = 20 K, 850 K
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
    type1 = None
    if "--type1" in argv:
        i = argv.index("--type1")
        type1 = os.path.abspath(argv[i + 1])
        argv = argv[:i] + argv[i + 2:]
    case = os.path.abspath(argv[0]) if argv else os.path.dirname(os.path.abspath(__file__))
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
    sp21, pk21 = abl_tools.spread_atc21()
    budget = sp21["T"][0]
    dp_tol = DP_T_TERMS_PER_K * budget + sp21["mdot_g"][0] / pk21["mdot_g"]
    sign_slack, slack1 = SIGN_SLACK, 0.0
    if os.environ.get("FUSS_GATE_SELFTEST", "0") == "1":
        budget, dp_tol, sign_slack, slack1 = 0.0, 0.0, -1.0, -1.0e9

    sT = bracket_series(case, pj, px, L, 1)
    sp = bracket_series(case, pj, px, L, 2)
    s1 = bracket_series(type1, pj, px, L, 1) if type1 else None
    ok = True
    print("ATC2-chemistryOff-type2 gate: FUSS type 2 (bracketing probes) vs PATO (%s), t = 0.1..60 s" % oracle_dir)

    # the type-1 run's own maximum and where it lies (for (acc) and (dir))
    peak1, peak_d, peak_t = None, None, None
    if s1 is not None:
        peak1 = 0.0
        for col, d in enumerate(DEPTHS_MM, 1):
            for row in pato_T:
                t = row[0]
                if t < 0.1 - 1e-9 or t > 60.0 + 1e-9:
                    continue
                v = interp_t(s1[d], t)
                if v is not None and abs(v - row[col]) > peak1:
                    peak1, peak_d, peak_t = abs(v - row[col]), d, t

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
            if peak_d is not None and d == peak_d and abs(t - peak_t) < 1e-6:
                dT_at_peak1 = dT
        if n == 0:
            print("FAIL: no FUSS sample inside PATO's window for depth %d mm" % d)
            return 1
        if worst > worst_all:
            worst_all, worst_at = worst, (d, t_worst)
        print("  %2d mm: max |dT| = %6.2f K at t = %5.1f s  (%d samples)" % (d, worst, t_worst, n))
    good = worst_all <= budget
    ok &= good
    print("  (T)   max over the seven depths: %.2f K at %s  (ATC 2.1 T spread %.2f K)  %s" % (worst_all, worst_at, budget, "PASS" if good else "FAIL"))
    if peak1 is None:
        print("  (acc) SKIP: no --type1 run given")
        print("  (dir) SKIP: no --type1 run given")
    else:
        good = worst_all <= peak1 + slack1
        ok &= good
        print("  (acc) the same maximum against the type-1 run's maximum %.2f K (%g mm, %g s)  %s"
              % (peak1, peak_d, peak_t, "PASS" if good else "FAIL"))
        if dT_at_peak1 is None:
            print("  (dir) FAIL: no sample at the type-1 peak (%g mm, %g s)" % (peak_d, peak_t))
            ok = False
        else:
            good = dT_at_peak1 < peak1 + slack1
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
    print("        ratio %.3f  (tol %.4f)  %s" % (worst_dp / dp_max if dp_max > 0 else float("inf"), dp_tol, "PASS" if good else "FAIL"))

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

    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
