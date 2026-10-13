#!/usr/bin/env python3
"""
verify_exact.py -- gates of the three exact-solution benchmarks of the type-2 closure (2026-10-10, made gates 2026-10-13).

Usage: verify_exact.py <kind> <case> [<case> ...] [--record]
  kind: pressure_transient | transpiration | sharp_front; the cases form one mesh ladder (coarsest first), each a run
  case of make_exact_case.py (case.json, OUTPUT/field*.tec); budgets in <kind>/expected.json beside this script.
Checks (exact solutions: darcy_exact.py; errors: exact_compare.py):
  pressure_transient  max|p - p_ex|/dp at the snapshot nearest t = 0.5 L^2/D: the finest error <= e_max and the order on
                      the finest pair >= order_min (the line solve is second order in space; the step is backward Euler,
                      first order in time, and dt ~ dx^2 at fixed vnn)
  transpiration       max|T - T_ex| and |mdot - mdot_ex|/mdot_ex at the steady state: finest <= eT_max, emdot_max, and
                      the T order on the finest pair >= order_min (red before 2026-10-13: the open face took the gas
                      enthalpy at the cell temperature, order 1.0)
  sharp_front         L2 T error and mass-based front-position error, averaged over t >= 5 s: finest <= eL2_max,
                      exf_max, and the front order on the finest pair inside [order_lo, order_hi] (a cell-by-cell front
                      is first order by design)
--record writes the measured values into expected.json (never the budgets). FUSS_GATE_SELFTEST=1 sets every budget to 0
and the order windows to empty (the proven-red hook). Exit 0 on PASS, 1 on FAIL. Python 3.6, no numpy.
"""
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tools"))
import darcy_tools as dt      # noqa: E402
import exact_compare as ec    # noqa: E402


def order(e1, e2, n1, n2):
    return math.log(e1 / e2) / math.log(n2 / n1) if e1 > 0 and e2 > 0 else float("nan")


def measure(kind, cases):
    out = []
    for c in cases:
        meta, rows = ec.compare(c)
        if meta["kind"] != kind:
            sys.exit("verify_exact.py: %s is a %s case, not %s" % (c, meta["kind"], kind))
        n = meta["ncells"]
        if kind == "pressure_transient":
            snaps = [(t, e) for t, e in rows if isinstance(t, float)]
            t, e = min(snaps, key=lambda r: abs(r[0] / meta["tau_p"] - 0.5))
            out.append({"n": n, "t": t, "e": e["e_p"]})
        elif kind == "transpiration":
            snaps = [(t, e) for t, e in rows if isinstance(t, float)]
            md = [e for t, e in rows if t == "mdot"][0]
            out.append({"n": n, "t": snaps[-1][0], "e": snaps[-1][1]["e_T"], "emdot": md["e_mdot_rel"]})
        else:
            late = [e for t, e in rows if isinstance(t, float) and t >= 5.0]
            out.append({"n": n, "e": sum(r["e_L2"] for r in late) / len(late),
                        "exf": sum(r["e_xf_rel"] for r in late) / len(late), "nsnap": len(late)})
    return out


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    record = "--record" in sys.argv
    if len(args) < 3:
        sys.exit(__doc__)
    kind, cases = args[0], [os.path.abspath(c) for c in args[1:]]
    path = os.path.join(HERE, kind, "expected.json")
    exp = json.load(open(path))
    res = measure(kind, cases)
    sel = dt.selftest()
    ok = True
    print("exact benchmark %s: %s" % (kind, "  ".join("N%d" % r["n"] for r in res)))
    for r in res:
        if kind == "pressure_transient":
            print("  %4d cells  t %.4g s  max|p - p_ex|/dp %.4e" % (r["n"], r["t"], r["e"]))
        elif kind == "transpiration":
            print("  %4d cells  t %.4g s  max|T - T_ex| %.4e K  mdot error %.3e" % (r["n"], r["t"], r["e"], r["emdot"]))
        else:
            print("  %4d cells  (%d snapshots, t >= 5 s)  L2 T error %.4f K  front error %.3e" % (r["n"], r["nsnap"], r["e"], r["exf"]))
    a, b = res[-2], res[-1]
    p_e = order(a["e"], b["e"], a["n"], b["n"])
    checks = []
    if kind in ("pressure_transient", "transpiration"):
        e_max = 0.0 if sel else exp["e_max"]
        omin = 1.0e9 if sel else exp["order_min"]
        checks.append(("error on the finest grid %.4e <= %.4e" % (b["e"], e_max), b["e"] <= e_max))
        checks.append(("order on the finest pair %.3f >= %g" % (p_e, omin), p_e >= omin))
        if kind == "transpiration":
            m_max = 0.0 if sel else exp["emdot_max"]
            checks.append(("mdot error on the finest grid %.3e <= %.3e" % (b["emdot"], m_max), b["emdot"] <= m_max))
    else:
        p_f = order(a["exf"], b["exf"], a["n"], b["n"])
        lo, hi = (1.0e9, -1.0e9) if sel else (exp["order_lo"], exp["order_hi"])
        e_max = 0.0 if sel else exp["eL2_max"]
        f_max = 0.0 if sel else exp["exf_max"]
        checks.append(("L2 T error on the finest grid %.4f K <= %.4f K" % (b["e"], e_max), b["e"] <= e_max))
        checks.append(("front error on the finest grid %.3e <= %.3e" % (b["exf"], f_max), b["exf"] <= f_max))
        checks.append(("front order on the finest pair %.3f in [%g, %g]" % (p_f, lo, hi), lo <= p_f <= hi))
        print("  T L2 order on the finest pair %.3f (printed, not gated)" % p_e)
    for text, good in checks:
        ok &= good
        print("  %s  %s" % (text, "PASS" if good else "FAIL"))
    if record:
        exp["measured"] = res
        with open(path, "w") as fh:
            json.dump(exp, fh, indent=1, sort_keys=True)
            fh.write("\n")
        print("  recorded the measured values in %s" % path)
    print("%s: %s" % (kind, "PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
