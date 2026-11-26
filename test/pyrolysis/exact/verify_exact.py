#!/usr/bin/env python3
"""
verify_exact.py -- gates of the three exact-solution benchmarks of the type-2 closure (2026-10-10, made gates 2026-10-13).

Usage: verify_exact.py <kind> <case> [<case> ...]
  kind: pressure_transient | transpiration | sharp_front; the cases form one mesh ladder (coarsest first), each a run
  case of make_exact_case.py (case.json, OUTPUT/field*.tec); the expected orders in <kind>/expected.json beside this
  script.
Checks (exact solutions: darcy_exact.py; errors: exact_compare.py), each by the convergence rule of tools/conv_rule.py
(plan 17 P17-D1): the error falls on every rung and its order on the finest pair lies in [0.8 p, 1.2 p]:
  pressure_transient  max|p - p_ex|/dp at the snapshot nearest t = 0.5 L^2/D, p = 2 (the line solve is second order in
                      space; the step is backward Euler, first order in time, and dt ~ dx^2 at fixed vnn)
  transpiration       max|T - T_ex| (p = 2; red before 2026-10-13: the open face took the gas enthalpy at the cell
                      temperature, order 1.0) and |mdot - mdot_ex|/mdot_ex (p = 2: the Darcy flux of that line solve)
  sharp_front         L2 T error and mass-based front-position error, averaged over t >= 5 s, p = 1 (a cell-by-cell
                      front is first order by design; the L2 T error follows the front)
FUSS_GATE_SELFTEST=1 empties the order windows (the proven-red hook). Exit 0 on PASS, 1 on FAIL. Python 3.6, no numpy.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tools"))
import conv_rule as cr        # noqa: E402
import exact_compare as ec    # noqa: E402


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
    if len(args) < 3:
        sys.exit(__doc__)
    kind, cases = args[0], [os.path.abspath(c) for c in args[1:]]
    path = os.path.join(HERE, kind, "expected.json")
    exp = json.load(open(path))
    res = measure(kind, cases)
    ok = True
    print("exact benchmark %s: %s" % (kind, "  ".join("N%d" % r["n"] for r in res)))
    for r in res:
        if kind == "pressure_transient":
            print("  %4d cells  t %.4g s  max|p - p_ex|/dp %.4e" % (r["n"], r["t"], r["e"]))
        elif kind == "transpiration":
            print("  %4d cells  t %.4g s  max|T - T_ex| %.4e K  mdot error %.3e" % (r["n"], r["t"], r["e"], r["emdot"]))
        else:
            print("  %4d cells  (%d snapshots, t >= 5 s)  L2 T error %.4f K  front error %.3e" % (r["n"], r["nsnap"], r["e"], r["exf"]))
    ns = [r["n"] for r in res]
    qty = {"pressure_transient": [("max|p - p_ex|/dp", "e")],
           "transpiration": [("max|T - T_ex| [K]", "e"), ("mdot error", "emdot")],
           "sharp_front": [("L2 T error [K]", "e"), ("front error", "exf")]}[kind]
    for label, key in qty:
        good, text = cr.check(label, [r[key] for r in res], ns, exp["order"][key])
        ok &= good
        print("  %s  %s" % (text, "PASS" if good else "FAIL"))
    print("%s: %s" % (kind, "PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
