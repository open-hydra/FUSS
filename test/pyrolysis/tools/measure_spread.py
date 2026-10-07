#!/usr/bin/env python3
"""
measure_spread.py -- PATO-vs-FIAT spread on AblationTestCase_1.0, the budget of plan 06 gate 4.4.

Implements plan 06 section 4.4 ("Budget rule"): the budget of the FUSS-vs-FIAT comparison is
1.5 x the maximum |T_PATO - T_FIAT| over 0..60 s and the six non-surface probes, measured on a
local PATO 3.1 run of tutorials/1D/AblationTestCase_1.0 (the surface column is the Dirichlet
input and is not gated).

Usage:
    measure_spread.py <Ta_plot> <Ta_surfacePatch> <FIAT/T> [--write-expected <json>]

PATO output layout (libPATOx IOModel/SampleFunction/SampleFunction.C:122-131, verified 2026-10-07):
one header line `//t(s) probe0(x,y,z) probe1(x,y,z) ...`, then one row per output time,
`t v0 v1 ...` (time first). With the tutorial's dictionaries:
    Ta_plot          data columns 1..7 = plotDict points y = 0.049, 0.048, 0.046, 0.042, 0.038,
                     0.034, 0.026 m = depths 1, 2, 4, 8, 12, 16, 24 mm below the heated face
    Ta_surfacePatch  data columns 1, 2 = surfacePatchDict points y = 0.05 (top) and y = 0 (back)
FIAT/T (data/ref/FIAT/T): no header, 61 rows t = 0..60 s, 8 columns
    t, T(0 mm surface), T(1), T(2), T(4), T(8), T(16), T(50 mm back face).
"Column n" below means the n-th DATA column after the time column.

Pairing (plan 06 section 4.4, established empirically 2026-09-15 and confirmed by the tutorial's
plot.py): FIAT col 1 <-> Ta_surfacePatch col 1 (surface); FIAT cols 2-5 <-> Ta_plot cols 1-4
(1, 2, 4, 8 mm); FIAT col 6 <-> Ta_plot col 6 (16 mm); FIAT col 7 <-> Ta_surfacePatch col 2
(back face). PATO is interpolated linearly in time onto FIAT's 1 s grid.

Prints, per probe, max |dT| over 0-60 s and its time, max |dT| over 5-60 s and dT at 60 s, then
the worst-probe sweep for t_cut = 0 / 5 / 10 / 20 / 30 s. With --write-expected writes
    {"budget_K": ceil(1.5 x max over the non-surface probes), "measured_spread_K": that max,
     "pato_run": {"date": ..., "ic_K": ...}, ...}
Python 3.6, no third-party modules.
"""

import argparse
import datetime
import json
import math
import os
import sys

PAIRS = [
    # label,            FIAT col, source file, source col
    ("0 mm (surface)",  1, "surf", 1),
    ("1 mm",            2, "plot", 1),
    ("2 mm",            3, "plot", 2),
    ("4 mm",            4, "plot", 3),
    ("8 mm",            5, "plot", 4),
    ("16 mm",           6, "plot", 6),
    ("50 mm (back)",    7, "surf", 2),
]
T_CUTS = [0.0, 5.0, 10.0, 20.0, 30.0]
T_END = 60.0


def die(msg):
    sys.stderr.write("[ERROR] measure_spread: %s\n" % msg)
    sys.exit(1)


def read_columns(path, min_cols):
    """Numeric rows of a `//`- or `#`-commented whitespace table; returns the list of rows."""
    rows = []
    with open(path) as fh:
        for ln, line in enumerate(fh, 1):
            s = line.strip()
            if not s or s.startswith("//") or s.startswith("#"):
                continue
            toks = s.split()
            try:
                vals = [float(t) for t in toks]
            except ValueError:
                die("%s line %d: non-numeric row %r" % (path, ln, s))
            if len(vals) < min_cols:
                die("%s line %d: %d columns, need at least %d" % (path, ln, len(vals), min_cols))
            rows.append(vals)
    if not rows:
        die("%s: no data rows" % path)
    ts = [r[0] for r in rows]
    if any(b <= a for a, b in zip(ts, ts[1:])):
        die("%s: time column not strictly increasing" % path)
    return rows


def interp(ts, vs, t):
    """Linear interpolation of (ts ascending, vs) at t; refuses extrapolation."""
    if t < ts[0] or t > ts[-1]:
        die("t = %g outside the PATO output range [%g, %g]" % (t, ts[0], ts[-1]))
    lo, hi = 0, len(ts) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if ts[mid] <= t:
            lo = mid
        else:
            hi = mid
    if t == ts[lo]:
        return vs[lo]
    if t == ts[hi]:
        return vs[hi]
    return vs[lo] + (vs[hi] - vs[lo]) * (t - ts[lo]) / (ts[hi] - ts[lo])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ta_plot")
    ap.add_argument("ta_surfacepatch")
    ap.add_argument("fiat_t")
    ap.add_argument("--write-expected", metavar="JSON", default=None)
    args = ap.parse_args()

    plot = read_columns(args.ta_plot, 8)
    surf = read_columns(args.ta_surfacepatch, 3)
    fiat = read_columns(args.fiat_t, 8)
    src = {"plot": plot, "surf": surf}
    tF = [r[0] for r in fiat if r[0] <= T_END]
    fiat = [r for r in fiat if r[0] <= T_END]
    if tF[0] != 0.0 or tF[-1] != T_END:
        die("FIAT time grid is %g..%g, expected 0..%g" % (tF[0], tF[-1], T_END))

    print("# PATO AblationTestCase_1.0 vs FIAT (%s)" % args.fiat_t)
    print("# Ta_plot %s (%d rows, t = %g..%g), Ta_surfacePatch %s (%d rows)"
          % (args.ta_plot, len(plot), plot[0][0], plot[-1][0], args.ta_surfacepatch, len(surf)))
    print("# PATO interpolated linearly onto FIAT's grid (%d points, 0..%g s). IC: FIAT %g K, PATO %g K."
          % (len(tF), T_END, fiat[0][1], plot[0][1]))
    print("%-16s %14s %6s %14s %6s %8s %10s %10s" % ("probe", "max|dT| 0-60s", "@t[s]", "max|dT| 5-60s", "@t[s]", "dT@60s", "T_FIAT@60", "T_PATO@60"))
    diffs = {}
    per_probe = []
    for label, fcol, sname, scol in PAIRS:
        rows = src[sname]
        ts = [r[0] for r in rows]
        vs = [r[scol] for r in rows]
        Pi = [interp(ts, vs, t) for t in tF]
        F = [r[fcol] for r in fiat]
        d = [p - f for p, f in zip(Pi, F)]
        diffs[label] = d
        i_full = max(range(len(d)), key=lambda i: abs(d[i]))
        idx5 = [i for i in range(len(d)) if tF[i] >= 5.0]
        i5 = max(idx5, key=lambda i: abs(d[i]))
        print("%-16s %14.2f %6.0f %14.2f %6.0f %8.2f %10.2f %10.2f"
              % (label, abs(d[i_full]), tF[i_full], abs(d[i5]), tF[i5], d[-1], F[-1], Pi[-1]))
        per_probe.append({"probe": label, "fiat_col": fcol, "pato_file": sname, "pato_col": scol,
                          "max_abs_dT_0_60": abs(d[i_full]), "t_of_max": tF[i_full],
                          "max_abs_dT_5_60": abs(d[i5]), "t_of_max_5_60": tF[i5], "dT_at_60": d[-1]})
    print("#")
    print("# Window sweep: worst-probe max|dT| over t_cut <= t <= 60 s (surface column excluded)")
    print("%8s %13s  probe @t[s]" % ("t_cut[s]", "worst|dT|[K]"))
    sweep = []
    gated = [p for p in PAIRS if p[0] != "0 mm (surface)"]
    for tc in T_CUTS:
        idx = [i for i in range(len(tF)) if tF[i] >= tc]
        best = None
        for label, _, _, _ in gated:
            d = diffs[label]
            i = max(idx, key=lambda i: abs(d[i]))
            if best is None or abs(d[i]) > best[0]:
                best = (abs(d[i]), label, tF[i])
        print("%8.0f %13.2f  %s @%.0f" % (tc, best[0], best[1], best[2]))
        sweep.append({"t_cut": tc, "worst_abs_dT": best[0], "probe": best[1], "t": best[2]})

    measured = max(max(abs(x) for x in diffs[label]) for label, _, _, _ in gated)
    budget = int(math.ceil(1.5 * measured - 1e-12))
    print("#")
    print("# measured spread (non-surface probes, 0-60 s) = %.4f K; budget = ceil(1.5 x %.4f) = %d K" % (measured, measured, budget))

    if args.write_expected:
        mtime = datetime.datetime.fromtimestamp(os.path.getmtime(args.ta_plot))
        dts = [b - a for a, b in zip([r[0] for r in plot], [r[0] for r in plot][1:])]
        dts.sort()
        ic_vals = plot[0][1:] + surf[0][2:3]          # every probe at t = 0 (uniform IC in a real run)
        if max(ic_vals) - min(ic_vals) > 0:
            print("# note: the t = 0 probe values are not uniform (%g..%g K); ic_K is the deepest Ta_plot probe"
                  % (min(ic_vals), max(ic_vals)))
        out = {
            "budget_K": budget,
            "measured_spread_K": measured,
            "pato_run": {
                "date": mtime.date().isoformat(),
                "ic_K": plot[0][-1],
                "ic_K_range": [min(ic_vals), max(ic_vals)],
                "output_dt_s": round(dts[len(dts) // 2], 9),
                "t_end_s": plot[-1][0],
                "ta_plot": os.path.abspath(args.ta_plot),
                "ta_surfacepatch": os.path.abspath(args.ta_surfacepatch),
            },
            "fiat_t": os.path.abspath(args.fiat_t),
            "budget_rule": "ceil(1.5 x max |T_PATO - T_FIAT| over 0-60 s, probes 1/2/4/8/16/50 mm)",
            "per_probe": per_probe,
            "window_sweep": sweep,
            "written": datetime.date.today().isoformat(),
        }
        with open(args.write_expected, "w") as fh:
            json.dump(out, fh, indent=1, sort_keys=True)
            fh.write("\n")
        print("# wrote %s" % args.write_expected)


if __name__ == "__main__":
    main()
