#!/usr/bin/env python3
"""
verify.py -- gate of the Zhang (2010) three-code decomposing slab (make_zhang2010.py case), 2026-10-13.

Usage: verify.py <case> [--record]
Reference: ref/zhang_fig2_2_T.csv (temperatures, K) and ref/zhang_fig2_3_mloss.csv (gas mass flux, kg/m2/s) of the three
codes of the thesis (MATLAB, UEL, UMATHT), extracted from the vector figures (ref/README.md: 0.16 K, 1.7e-6 kg/m2/s).
FUSS: exposed and unexposed surface temperatures from wall*.tec (faces 1 and 2), the mid-plane from the two cells around
it (probes midlo, midhi), the gas mass flux at the exposed face from the first cell's mdot_pg (type 1 outflow, face 1).
Checks, for t > 20 s:
  (T)    max |T_FUSS - T_MATLAB| at each of the three planes <= dT_max (the three codes themselves spread 4.75-5.26 K;
         red before 2026-10-13, when BC 305's wall update dropped h_conv: 7.4 K at the mid-plane)
  (mdot) the FUSS peak gas flux within [mdot_lo, mdot_hi] x the lowest/highest reference peak, its time within dt_peak
--record writes the measured values into expected.json. FUSS_GATE_SELFTEST=1 zeroes the budgets. Python 3.6, no numpy.
"""
import csv
import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import darcy_tools as dt   # noqa: E402
import fuss_out as fo      # noqa: E402

CODES = ("MATLAB", "UEL", "UMATHT")


def read_ref(name):
    rows = list(csv.DictReader(open(os.path.join(HERE, "ref", name))))
    return {k: [(float(r["time_s"]), float(r[k])) for r in rows if r[k] not in ("", None)] for k in rows[0] if k != "time_s"}


def interp(series, t):
    for (t0, v0), (t1, v1) in zip(series[:-1], series[1:]):
        if t0 <= t <= t1:
            return v0 + (v1 - v0) * (t - t0) / (t1 - t0)
    return None


def fuss_series(case):
    ex = [(t, T) for t, T, q in fo.wall_series(case, "B1F1")]
    un = [(t, T) for t, T, q in fo.wall_series(case, "B1F2")]
    lo, hi = fo.probe(case, "midlo"), fo.probe(case, "midhi")
    mid = [(a[0], 0.5 * (a[1] + b[1])) for a, b in zip(lo, hi)]
    md = []
    for f in glob.glob(os.path.join(case, "OUTPUT", "field*.tec")):
        dims, t, d = dt.read_tec(f)
        md.append((t, d["mdot_pg"][0]))
    return {"exposed": ex, "middle": mid, "unexposed": un, "mdot": sorted(md)}


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if not args:
        sys.exit(__doc__)
    case = os.path.abspath(args[0])
    path = os.path.join(HERE, "expected.json")
    exp = json.load(open(path))
    sel = dt.selftest()
    rT, rM = read_ref("zhang_fig2_2_T.csv"), read_ref("zhang_fig2_3_mloss.csv")
    s = fuss_series(case)
    ok = True
    meas = {}
    print("Zhang (2010) decomposing slab: %s against the MATLAB code (UEL, UMATHT for the spread)" % case)
    dT_max = 0.0 if sel else exp["dT_max"]
    for key in ("exposed", "middle", "unexposed"):
        d = [v - interp(rT["MATLAB_%s_K" % key], t) for t, v in s[key] if t > 20.0 and interp(rT["MATLAB_%s_K" % key], t) is not None]
        spread = max(abs(v - interp(rT["%s_%s_K" % (c, key)], t)) for c in ("UEL", "UMATHT") for t, v in rT["MATLAB_%s_K" % key]
                     if t > 20.0 and rT.get("%s_%s_K" % (c, key)) and interp(rT["%s_%s_K" % (c, key)], t) is not None)
        worst = max(abs(v) for v in d)
        good = worst <= dT_max
        ok &= good
        meas[key] = {"max": worst, "mean": sum(d) / len(d)}
        print("  (T) %-9s max|T - T_MATLAB| %.2f K (mean %+.2f K; codes' spread %.2f K; budget %.2f K)  %s"
              % (key, worst, sum(d) / len(d), spread, dT_max, "PASS" if good else "FAIL"))
    pk = max(s["mdot"], key=lambda r: r[1])
    refpk = [max(rM["%s_kg_s_m2" % c], key=lambda r: r[1]) for c in CODES]
    lo = (1.0e9 if sel else exp["mdot_lo"]) * min(r[1] for r in refpk)
    hi = (0.0 if sel else exp["mdot_hi"]) * max(r[1] for r in refpk)
    tmid = sum(r[0] for r in refpk) / len(refpk)
    dt_ok = (-1.0 if sel else exp["dt_peak"])
    good = lo <= pk[1] <= hi and abs(pk[0] - tmid) <= dt_ok
    ok &= good
    meas["mdot_peak"] = pk
    print("  (mdot) FUSS peak %.4e kg/m2/s at %.0f s; the codes %s; window [%.4e, %.4e], time within %g s of %.0f s  %s"
          % (pk[1], pk[0], ", ".join("%s %.4e at %.0f s" % (c, r[1], r[0]) for c, r in zip(CODES, refpk)), lo, hi, dt_ok, tmid,
             "PASS" if good else "FAIL"))
    if "--record" in sys.argv:
        exp["measured"] = meas
        with open(path, "w") as fh:
            json.dump(exp, fh, indent=1, sort_keys=True)
            fh.write("\n")
        print("  recorded the measured values in %s" % path)
    print("zhang2010: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
