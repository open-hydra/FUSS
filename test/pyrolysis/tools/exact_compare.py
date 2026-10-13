#!/usr/bin/env python3
"""exact_compare.py -- compare FUSS runs of make_exact_case.py cases with darcy_exact.py.

Usage: exact_compare.py <case-dir> [<case-dir> ...]   (each a run case: case.json, OUTPUT/field*.tec)
Prints, per case, the error against the exact solution at every snapshot, and, when several cases
of the same kind are given in refinement order, the observed orders between successive cases.
Errors: pressure_transient max_i |p - p_ex| / (p_w - p0); transpiration max_i |T - T_ex| [K] and
|mdot - mdot_ex|/mdot_ex; sharp_front max_i |T - T_ex| [K], front position error and the
surface gas flux error. Python 3.6, no numpy.
"""
import glob
import json
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import darcy_tools as dt    # noqa: E402
import darcy_exact as dx    # noqa: E402


def snapshots(case):
    files = glob.glob(os.path.join(case, "OUTPUT", "field*.tec"))
    out = []
    for f in files:
        dims, t, d = dt.read_tec(f)
        out.append((t, f, dims, d))
    out.sort(key=lambda r: r[0])
    return out


def row_i(dims, d, name):
    """Values of cell row j = 1, k = 1 along i."""
    ni = dims[0] - 1
    return d[name][:ni]


def centres(dims, d):
    ni = dims[0] - 1
    x = d["x"]
    return [0.5 * (x[i] + x[i + 1]) for i in range(ni)]


def err_pressure_transient(meta, t, xc, p):
    L, D, p0, pw = meta["L"], meta["D"], meta["p0"], meta["pw"]
    return max(abs(p[i] - dx.pressure_transient_linear(xc[i], t, L, D, p0, pw)) for i in range(len(p))) / abs(pw - p0)


def compare(case):
    meta = json.load(open(os.path.join(case, "case.json")))
    kind = meta["kind"]
    snaps = snapshots(case)
    rows = []
    if kind == "pressure_transient":
        for t, f, dims, d in snaps:
            xc = centres(dims, d)
            p = row_i(dims, d, "p_pg")
            rows.append((t, {"e_p": err_pressure_transient(meta, t, xc, p)}))
    elif kind == "transpiration":
        md, A, B, Pe, lin = dx.transpiration_steady(meta["L"], meta["k"], meta["cpg"], meta["W"], meta["mu"],
                                                    meta["K"], meta["pin"], meta["pout"], meta["h1"], meta["T1"],
                                                    meta["h2"], meta["T2"])
        for t, f, dims, d in snaps:
            xc = centres(dims, d)
            T = row_i(dims, d, "T")
            p = row_i(dims, d, "p_pg")
            eT = max(abs(T[i] - dx.transpiration_T(xc[i], meta["L"], A, B, Pe, lin)) for i in range(len(T)))
            ep = max(abs(p[i] - dx.transpiration_p(xc[i], meta["L"], md, meta["W"], meta["mu"], meta["K"], meta["pin"],
                                                   A, B, Pe, lin)) for i in range(len(p)))
            # mdot from the budget's outflow rate at the last two rows (steady): M_gas_out per area
            rows.append((t, {"e_T": eT, "e_p_Pa": ep}))
        # mdot from the last snapshot: FUSS's cell mass-flux magnitude |Gamma grad p| (uniform at steady state)
        t, f, dims, d = snaps[-1]
        mf = row_i(dims, d, "mdot_pg")
        rate = sum(mf[1:-1]) / (len(mf) - 2)
        rows.append(("mdot", {"mdot_fuss": rate, "mdot_exact": md, "e_mdot_rel": abs(rate - md) / md}))
    elif kind == "sharp_front":
        ac, av = meta["alpha_c"], meta["alpha_v"]
        lam = meta["lam"]
        for t, f, dims, d in snaps:
            if t <= 0.0:
                continue
            xc = centres(dims, d)
            T = row_i(dims, d, "T")
            Tex = [dx.sharp_front_T(x, t, lam, ac, av, meta["rho_v"], meta["rho_c"], meta["cp"], meta["cpg"],
                                    meta["Ts"], meta["Tp"], meta["T0"]) for x in xc]
            e = [T[i] - Tex[i] for i in range(len(T))]
            h = xc[1] - xc[0]
            eT = max(abs(v) for v in e)
            eL2 = math.sqrt(sum(v * v for v in e) * h / meta["L"])
            # front position from the decomposed mass (sub-cell): sum (rho_v - rho_s) dx / (rho_v - rho_c)
            rs = row_i(dims, d, "rho_s")
            xf = sum((meta["rho_v"] - r) for r in rs) * h / (meta["rho_v"] - meta["rho_c"])
            xf_ex = dx.sharp_front_position(t, lam, ac)
            rows.append((t, {"e_T": eT, "e_L2": eL2, "xf": xf, "xf_ex": xf_ex, "e_xf_rel": abs(xf - xf_ex) / xf_ex}))
    return meta, rows


def meta_kind_is_sf(m):
    return m.get("kind") == "sharp_front"


def main():
    cases = sys.argv[1:]
    if not cases:
        sys.exit(__doc__)
    results = []
    for c in cases:
        meta, rows = compare(c)
        results.append((c, meta, rows))
        print("== %s  (%s, %d cells, vnn %g)" % (c, meta["kind"], meta["ncells"], meta["vnn"]))
        if meta["kind"] == "sharp_front":
            late = [e for t, e in rows if t >= 5.0]
            if late:
                print("   time averages over %d snapshots, t >= 5 s:  max-norm %.4f K  L2 %.4f K  front %.3e" % (
                      len(late), sum(e["e_T"] for e in late) / len(late), sum(e["e_L2"] for e in late) / len(late),
                      sum(e["e_xf_rel"] for e in late) / len(late)))
                results[-1] = (c, meta, rows + [("avg", {"e_Tavg": sum(e["e_T"] for e in late) / len(late),
                                                       "e_L2avg": sum(e["e_L2"] for e in late) / len(late),
                                                       "e_xfavg": sum(e["e_xf_rel"] for e in late) / len(late)})])
            continue
        for t, e in rows:
            if t == "mdot":
                print("   mdot FUSS %.10e  exact %.10e  rel err %.3e" % (e["mdot_fuss"], e["mdot_exact"], e["e_mdot_rel"]))
                continue
            print("   t %12.6g  " % t + "  ".join("%s %s" % (k, ("%.4e" % v) if isinstance(v, float) else v) for k, v in e.items()))
    # orders between successive cases at matching snapshot index
    if len(results) > 1:
        print("== observed orders between successive cases (same snapshot index)")
        for (c1, m1, r1), (c2, m2, r2) in zip(results[:-1], results[1:]):
            ratio_h = m2["ncells"] / m1["ncells"]
            if meta_kind_is_sf(m1):
                a1 = [e for t, e in r1 if t == "avg"]
                a2 = [e for t, e in r2 if t == "avg"]
                if a1 and a2:
                    e1, e2 = a1[0], a2[0]
                    print("   %s -> %s  time-averaged  max-norm order %.3f  L2 order %.3f  front order %.3f" % (
                          os.path.basename(c1), os.path.basename(c2),
                          math.log(e1["e_Tavg"] / e2["e_Tavg"]) / math.log(ratio_h),
                          math.log(e1["e_L2avg"] / e2["e_L2avg"]) / math.log(ratio_h),
                          math.log(e1["e_xfavg"] / e2["e_xfavg"]) / math.log(ratio_h)))
                continue
            for (t1, e1), (t2, e2) in zip(r1, r2):
                if t1 == "avg":
                    print("   %s -> %s  time-averaged  max-norm order %.3f  L2 order %.3f  front order %.3f" % (
                          os.path.basename(c1), os.path.basename(c2),
                          math.log(e1["e_Tavg"] / e2["e_Tavg"]) / math.log(ratio_h),
                          math.log(e1["e_L2avg"] / e2["e_L2avg"]) / math.log(ratio_h),
                          math.log(e1["e_xfavg"] / e2["e_xfavg"]) / math.log(ratio_h)))
                    continue
                if meta_kind_is_sf(m1):
                    continue
                if t1 == "mdot":
                    if e1["e_mdot_rel"] > 0 and e2["e_mdot_rel"] > 0:
                        print("   %s -> %s  mdot  order %.3f" % (os.path.basename(c1), os.path.basename(c2),
                              math.log(e1["e_mdot_rel"] / e2["e_mdot_rel"]) / math.log(ratio_h)))
                    continue
                parts = []
                for k in e1:
                    if k.startswith("e_") and e1[k] and e2[k] and e1[k] > 0 and e2[k] > 0 and not math.isnan(e1[k]):
                        parts.append("%s %.3f" % (k, math.log(e1[k] / e2[k]) / math.log(ratio_h)))
                print("   %s -> %s  t %.6g/%.6g  %s" % (os.path.basename(c1), os.path.basename(c2), t1, t2, "  ".join(parts)))


if __name__ == "__main__":
    main()
