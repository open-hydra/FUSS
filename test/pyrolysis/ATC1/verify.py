#!/usr/bin/env python3
"""Gate 4.4 of plan 06, re-anchored by plan 17 (2026-11-26) on PATO: Ablation Test Case 1 against PATO 3.1.

GATE (plan 17, the lead's reading (a) of P17-D2): max over the six probes (1, 2, 4, 8, 16 mm, the 50 mm back face) and
FIAT's 1 s grid t = 1..60 s of |T_FUSS - T_PATO| <= the PATO-FIAT spread, max over the same probes and grid (0..60 s)
of |T_PATO - T_FIAT|, computed here with tools/measure_spread.py's pairing from ref/pato-298K (PATO 3.1, the unmodified
tutorial at 298 K; $FUSS_PATO_ATC1_OUT replaces that directory) and ref/FIAT_T.dat (13.18 K,
ref/pato-298K/spread_298.txt). FUSS_GATE_SELFTEST=1 sets the tolerance to 0. The text below, from "Gate" to the end
of the Richardson paragraph, is the history of the anchor this replaced (P17-D8: no Richardson check is gated).

A 0.05 m TACOT slab, 100 cells graded 0.1 toward the heated face, whose
temperature follows Tw_time.dat (298 -> 1644 K in 0.1 s, held to 60 s); the
back face is adiabatic; T(0) = 298 K (the ATC1 specification, [ATC1] v1.1;
PATO's tutorial starts at 300 K). FUSS samples each FIAT depth with TWO
index-position probes whose cell centres bracket x = L - d (probes.json, written
by make_case.py); the value at the exact depth is the linear interpolation in x
between them, then linear in time. The 50 mm back face is the first cell.

Gate (expected.json, verification record plan-bucket/records/2026-10-08-verification):
the anchor is f_star, the Richardson extrapolation of FUSS itself on 100/200/400
graded cells (observed order 0.7-1.0 behind the pyrolysis front, 1.6-2.3 ahead of
it), at FIAT's instants t = 1..60 s and the depths 1, 2, 4, 8, 16, 50 mm where the
three meshes are asymptotic. For every such point
    |(T_FUSS - f_star) - f100_minus_fstar| <= band_K
i.e. the 100-cell discretisation error measured on 2026-10-08 (7.4 K at 1 mm /
1 s, 2.5-4 K deep and late) may move by at most band_K (1 K). A change of the
model (plan 07 type 2, plan 08) or of the mesh/probe layout needs a new anchor.

Reported, not gated: max |T_FUSS - T_FIAT| over the same window (the former
budget), f_star - FIAT, and the 98 % / 2 % rho_s front depths at the last field
write with mdot_pg of the surface cell. FIAT (bundled with PATO 3.1, ref/FIAT_T.dat,
columns t, T(0 mm), T(1), T(2), T(4), T(8), T(16), T(50 mm)) is 4-17 K cooler than
the converged FUSS during each front passage; the exact oracles of the record put
the type-2 gas-storage shift at 2-4 K and the table-hbar data difference at 0.15 K,
so the FIAT deficit cannot be attributed without running FIAT.

Usage: verify.py [case_dir]   Exit 0 on PASS, 1 on FAIL. Python 3.6, no numpy.
"""
import json, os, re, sys

DEPTHS_MM = [1, 2, 4, 8, 16]            # bracketed by probe pairs
FIAT_COL = {1: 2, 2: 3, 4: 4, 8: 5, 16: 6, 50: 7}


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


def read_cols(path):
    """a PATO probe file (scalar/Ta_plot, Ta_surfacePatch): a //-header line, then rows t v1 .. vn"""
    rows = []
    for line in open(path):
        s = line.strip()
        if not s or s.startswith("//"):
            continue
        try:
            rows.append([float(x) for x in s.split()])
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

    # the gate: FUSS against PATO on FIAT's grid, the PATO-FIAT spread as tolerance (measure_spread.py's pairing)
    pato_dir = os.environ.get("FUSS_PATO_ATC1_OUT") or os.path.join(here, "ref", "pato-298K")
    Ppl = read_cols(os.path.join(pato_dir, "Ta_plot"))
    Psu = read_cols(os.path.join(pato_dir, "Ta_surfacePatch"))
    pato = {d: [(r[0], r[c]) for r in rows] for d, (rows, c) in
            {1: (Ppl, 1), 2: (Ppl, 2), 4: (Ppl, 3), 8: (Ppl, 4), 16: (Ppl, 6), 50: (Psu, 2)}.items()}
    spread = max((abs(interp_t(pato[d], row[0]) - row[FIAT_COL[d]]), d, row[0])
                 for d in pato for row in fiat if 0.0 <= row[0] <= 60.0)
    tol = 0.0 if os.environ.get("FUSS_GATE_SELFTEST", "0") == "1" else spread[0]
    worst = (0.0, None, None)
    for d in DEPTHS_MM + [50]:
        for row in fiat:
            t = row[0]
            if t < 1.0 or t > 60.0:
                continue
            v = interp_t(series[d], t)
            if v is None and series[d] and 0.0 <= t - series[d][-1][0] <= 0.2:
                v = series[d][-1][1]        # the run ends within one probe interval of t (60 s): last sample
            if v is None:
                print("FAIL: no FUSS sample at t = %g s for depth %d mm" % (t, d)); return 1
            e = abs(v - interp_t(pato[d], t))
            if e > worst[0]:
                worst = (e, d, t)
    ok = worst[0] <= tol
    print("ATC1 gate: max over the six probes and FIAT's 1 s grid of |T_FUSS - T_PATO| = %.2f K at %s mm, %g s; "
          "tolerance = the PATO-FIAT spread %.2f K (at %s mm, %g s; PATO %s)  %s"
          % (worst[0], worst[1], worst[2], tol, spread[1], spread[2], pato_dir, "PASS" if ok else "FAIL"))

    # informational: the former metric against FIAT
    worst_all = 0.0
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
        worst_all = max(worst_all, worst)
        print("  (not gated) %2d mm: max |T_FUSS - T_FIAT| = %6.2f K at t = %4.0f s  (%d samples)" % (d, worst, t_worst, n))
    print("  (not gated) max over the six probes vs FIAT: %.2f K" % worst_all)

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

    print("ATC1: PASS" if ok else "ATC1: FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
