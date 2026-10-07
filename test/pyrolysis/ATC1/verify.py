#!/usr/bin/env python3
"""Gate 4.4 of plan 06: Ablation Test Case 1 against the FIAT reference.

A 0.05 m TACOT slab, 100 cells graded 0.1 toward the heated face, whose
temperature follows Tw_time.dat (298 -> 1644 K in 0.1 s, held to 60 s); the
back face is adiabatic; T(0) = 298 K (the ATC1 specification, [ATC1] v1.1;
PATO's tutorial starts at 300 K). The oracle is the FIAT solution bundled with
PATO 3.1 (tutorials/1D/AblationTestCase_1.0/data/ref/FIAT/T: t = 0..60 s every
1 s; columns t, T(0 mm), T(1), T(2), T(4), T(8), T(16), T(50 mm)); a copy is
in ref/FIAT_T.dat. The surface column is the Dirichlet input and is not gated.

FUSS samples each FIAT depth with TWO index-position probes whose cell centres
bracket x = L - d (probes.json, written by make_case.py); the value at the exact
depth is the linear interpolation in x between them, then linear in time onto
FIAT's 1 s grid. The 50 mm back face is the first cell (adiabatic face, the
cell value is the face value to well below 0.1 K).

Pass: max_{t = 1..60 s, depths 1, 2, 4, 8, 16, 50 mm} |T_FUSS - T_FIAT| <= budget_K
(expected.json; 21 K = 1.5 x the 13.76 K PATO-vs-FIAT spread measured on
2026-09-15 at PATO's 300 K IC, plan 06 section 4.4; the 298 K re-measurement of
2026-10-07 gave 13.18 K, ref/pato-298K/spread_298.txt).

Also reported, not gated: the 98 % / 2 % rho_s front depths at the last field
write and mdot_pg of the surface cell.

Usage: verify.py [case_dir]   Exit 0 on PASS, 1 on FAIL. Python 3.6, no numpy.
"""
import json, os, re, sys

DEPTHS_MM = [1, 2, 4, 8, 16]            # bracketed by probe pairs
FIAT_COL = {1: 2, 2: 3, 4: 4, 8: 5, 16: 6, 50: 7}
BUDGET_DEFAULT = 21.0


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


def read_fiat(path):
    rows = []
    for line in open(path):
        p = line.split()
        if len(p) >= 8:
            try:
                rows.append([float(x) for x in p])
            except ValueError:
                pass
    return rows


def read_tec(path):
    with open(path) as fh:
        lines = fh.read().split("\n")
    names = [t.strip('"') for t in lines[0].split("=", 1)[1].split() if t.strip('"')]
    m = re.search(r"I=(\d+),\s*J=(\d+),\s*K=(\d+)", lines[1])
    I, J, K = (int(x) for x in m.groups())
    vals = []
    for line in lines[2:]:
        s = line.strip()
        if not s:
            continue
        try:
            vals.append(float(s))
        except ValueError:
            pass
    n_nodal = I * J * K
    n_cell = (I - 1) * (J - 1) * max(K - 1, 1)
    data, off = {}, 0
    for name in names:
        n = n_nodal if name in ("x", "y", "z") else n_cell
        data[name] = vals[off:off + n]
        off += n
    return (I - 1, J - 1, max(K - 1, 1)), data


def probe_x(pj):
    """name -> cell-centre x from probes.json (list of dicts or dict)."""
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
    case = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
    here = os.path.dirname(os.path.abspath(__file__))
    pj = json.load(open(os.path.join(case, "probes.json")))
    L = pj["L"]
    px = probe_x(pj)
    fiat_path = os.path.join(here, "ref", "FIAT_T.dat")
    if not os.path.isfile(fiat_path) and os.environ.get("PATO_DIR"):
        fiat_path = os.path.join(os.environ["PATO_DIR"], "tutorials", "1D", "AblationTestCase_1.0", "data", "ref", "FIAT", "T")
    fiat = read_fiat(fiat_path)
    expected_path = os.path.join(here, "expected.json")
    budget = json.load(open(expected_path)).get("budget_K", BUDGET_DEFAULT) if os.path.isfile(expected_path) else BUDGET_DEFAULT

    # FUSS series at the exact depths
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
    series[50] = read_series(os.path.join(case, "OUTPUT", "back.txt"))

    print("ATC1 gate: FUSS (bracketing probes, exact depths) vs FIAT, t = 1..60 s, budget %.1f K" % budget)
    worst_all = 0.0
    ok = True
    for d in DEPTHS_MM + [50]:
        worst, t_worst, n = 0.0, None, 0
        for row in fiat:
            t = row[0]
            if t < 1.0 or t > 60.0:
                continue
            v = interp_t(series[d], t)
            if v is None:
                continue
            n += 1
            dT = abs(v - row[FIAT_COL[d]])
            if dT > worst:
                worst, t_worst = dT, t
        if n == 0:
            print("FAIL: no FUSS sample inside FIAT's window for depth %d mm" % d); return 1
        worst_all = max(worst_all, worst)
        print("  %2d mm: max |dT| = %6.2f K at t = %4.0f s  (%d samples)" % (d, worst, t_worst, n))
    print("  max over the six probes: %.2f K  (budget %.1f K)" % (worst_all, budget))
    ok = worst_all <= budget

    # not gated: pyrolysis front and surface blowing at the last field write
    try:
        files = sorted((f for f in os.listdir(os.path.join(case, "OUTPUT")) if re.match(r"field\d*\.tec$", f)),
                       key=lambda f: int(re.sub(r"\D", "", f) or 10 ** 9))
        (ni, nj, nk), fd = read_tec(os.path.join(case, "OUTPUT", files[-1]))
        rho_v, rho_c = 280.0, 220.0
        xs = [px[n] for n in sorted(px, key=lambda n: px[n])]
        rho = [fd["rho_s"][i] for i in range(ni)]          # j = 1, k = 1 line
        xc = [(fd["x"][i] + fd["x"][i + 1]) * 0.5 for i in range(ni)]
        def front(frac):
            thr = rho_c + frac * (rho_v - rho_c)
            for i in range(ni - 1, -1, -1):
                if rho[i] > thr:
                    return L - xc[i]
            return L
        print("  (not gated) %s: rho_s 98 %% front %.2f mm, 2 %% front %.2f mm below the surface; surface mdot_pg = %.3e kg/m2/s"
              % (files[-1], front(0.98) * 1e3, front(0.02) * 1e3, fd["mdot_pg"][ni - 1]))
    except Exception as e:  # reporting only
        print("  (not gated) front report unavailable: %s" % e)

    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
