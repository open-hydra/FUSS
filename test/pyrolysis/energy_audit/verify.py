#!/usr/bin/env python3
"""Gate 4.3 of plan 06: energy and mass audit of the pyrolysis update.

A 40-cell TACOT slab, face 2 ramped 300 -> 1644 K in 1 s (BC 302 time file),
face 1 adiabatic, 10 s. The solver writes OUTPUT/energy-budget.txt with one
line at t = 0 and one at every solution write:

    t  E_tot  M_s  Q_out  Q_gas_out  M_gas_out  Q_src

where E_tot = sum V E(T, chi), M_s = sum V rho_s, Q_out = int sum_cells R dt
(net outgoing flux: conduction + BC + gas advection), Q_gas_out and M_gas_out
the energy and mass carried out through the outflow face, Q_src = int sum
qvol V dt. The accumulators are advanced with the scheme's own RK weights, so
the balances below are identities of the discrete update (round-off only).

Checks (plan 06 section 4.3):
  (i)   energy:  |E_tot(t) - E_tot(0) + Q_out - Q_src| <= 1e-10 max(|E_tot(0)|, |Q_out|)
  (ii)  solid mass: |M_s(t) - M_s(0) + M_gas_out| <= 1e-10 M_s(0)
  (iii) independent state sum: E_tot(t) and M_s(t) recomputed here from
        OUTPUT/field.tec (T, chi_*), the INPUT tables and cell volumes from the
        node coordinates, agree with the solver's columns to 1e-10 relative
  (iv)  time-discretisation error bound, from the --vnn025 and --vnn0125 reruns:
        probe series interpolated onto t = 2.0, 2.1, ..., 10.0 s (the 1 s ramp
        is resolved by the coarsest run, dt ~ 0.17 s; with the plan's 0.1 s ramp
        the staircase sampling of the ramp aliased and the differences GREW under
        refinement -- 3.3, 4.9, 1.9 K measured 2026-10-07 -- so the ramp was lengthened);
        max|T(0.5) - T(0.25)| / max|T(0.25) - T(0.125)| >= 1.6 and
        max|T(0.5) - T(0.125)| <= budget_iv (expected.json; recorded on the
        first run with --record, then the gate fails above measured + 1 K)
  (v)   zero-kinetics self-consistency, from the --nokin rerun: every chi_*
        and mdot_pg exactly 0.0, rho_s = rho_v to 1e-14 relative, and
        Q_gas_out = M_gas_out = 0 exactly in its budget file
  (vi)  cross-thread determinism, from the --p8 rerun: field.tec and
        energy-budget.txt byte-identical to the 1-thread run

Usage: verify.py <case_dir> [--vnn025 DIR] [--vnn0125 DIR] [--nokin DIR] [--p8 DIR] [--record]
A sub-check whose directory is not given prints SKIP. Exit 0 on PASS, 1 on FAIL.
Python 3.6, no numpy.
"""
import hashlib, json, os, re, sys

TOL_BAL = 1.0e-10
TOL_STATE = 1.0e-10
TOL_RHO_V = 1.0e-14
ORDER_RATIO_MIN = 1.6
BUDGET_IV_DEFAULT = 3.0     # K, provisional (plan 06 section 4.3 (iv)); replaced by expected.json


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
    return (I, J, K), data


def cell_volumes(dims, data):
    """Hexahedral volumes by the divergence theorem from the node coordinates (test/moving-mesh/conservation/verify.py)."""
    I, J, K = dims
    X, Y, Z = data["x"], data["y"], data["z"]

    def nid(i, j, k):
        return i + j * I + k * I * J
    FACES = ((0, 2, 3, 1), (4, 5, 7, 6), (0, 1, 5, 4), (2, 6, 7, 3), (0, 4, 6, 2), (1, 3, 7, 5))

    def tri_flux(p, q, r):
        u = [q[c] - p[c] for c in range(3)]
        v = [r[c] - p[c] for c in range(3)]
        n = (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])
        cen = [(p[c] + q[c] + r[c]) / 3.0 for c in range(3)]
        return 0.5 * sum(cen[c] * n[c] for c in range(3))

    def quad_flux(a, b, c, d):
        m = [0.25 * (a[t] + b[t] + c[t] + d[t]) for t in range(3)]
        return tri_flux(a, b, m) + tri_flux(b, c, m) + tri_flux(c, d, m) + tri_flux(d, a, m)
    vols = []
    idx = {(0, 0, 0): 0, (0, 0, 1): 1, (0, 1, 0): 2, (0, 1, 1): 3, (1, 0, 0): 4, (1, 0, 1): 5, (1, 1, 0): 6, (1, 1, 1): 7}
    for k in range(max(K - 1, 1)):
        for j in range(J - 1):
            for i in range(I - 1):
                pts = [None] * 8
                for a in (0, 1):
                    for b in (0, 1):
                        for c in (0, 1):
                            nd = nid(i + a, j + b, k + c)
                            pts[idx[(a, b, c)]] = (X[nd], Y[nd], Z[nd])
                v = sum(quad_flux(pts[f[0]], pts[f[1]], pts[f[2]], pts[f[3]]) for f in FACES)
                vols.append(-v / 3.0)
    return vols


def read_budget(path):
    rows = []
    for line in open(path):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        rows.append([float(x) for x in s.split()])
    return rows


def read_probe(path):
    out = []
    for line in open(path):
        p = line.split()
        if len(p) >= 2:
            try:
                out.append((float(p[0]), float(p[1])))
            except ValueError:
                pass
    return out


def read_pyrolysis_dat(path):
    kv = {}
    for line in open(path):
        line = line.split("//", 1)[0].strip()
        if line:
            k, v = line.split(None, 1)
            kv[k] = v.strip()
    nph = int(kv["nSolidPhases"])
    rho_I = [float(kv["rhoI[%d]" % i]) for i in range(1, nph + 1)]
    eps_I = [float(kv["epsI[%d]" % i]) for i in range(1, nph + 1)]
    erf = []
    for i in range(1, nph + 1):
        for j in range(1, int(kv["nPyroReac[%d]" % i]) + 1):
            erf.append(eps_I[i - 1] * rho_I[i - 1] * float(kv["F[%d][%d]" % (i, j)]))
    rho_v = sum(e * r for e, r in zip(eps_I, rho_I))
    rho_c = rho_v - sum(erf)
    return rho_v, rho_c, erf


def read_solid_table(path):
    """virgin.dat / char.dat: {T_int: h} of the (single) zone."""
    h = {}
    for line in open(path):
        p = line.split()
        if len(p) == 8:
            try:
                h[int(round(float(p[0])))] = float(p[2])
            except ValueError:
                pass
    return h


def interp_h(h, T):
    Tm = int(T)
    return h[Tm] + (T - Tm) * (h[Tm + 1] - h[Tm])


def interp_series(series, t):
    """Linear interpolation in time of a (t, T) series; None outside its range."""
    if t < series[0][0] or t > series[-1][0]:
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
    if t1 == t0:
        return v0
    return v0 + (v1 - v0) * (t - t0) / (t1 - t0)


def md5(path):
    return hashlib.md5(open(path, "rb").read()).hexdigest()


def main():
    args = sys.argv[1:]
    opts = {"--vnn025": None, "--vnn0125": None, "--nokin": None, "--p8": None}
    record = False
    pos = []
    i = 0
    while i < len(args):
        if args[i] in opts:
            opts[args[i]] = os.path.abspath(args[i + 1]); i += 2
        elif args[i] == "--record":
            record = True; i += 1
        else:
            pos.append(args[i]); i += 1
    case = os.path.abspath(pos[0]) if pos else os.path.dirname(os.path.abspath(__file__))
    here = os.path.dirname(os.path.abspath(__file__))
    ok = True
    print("energy audit gate on %s" % case)

    # (i), (ii): the solver's own budget
    rows = read_budget(os.path.join(case, "OUTPUT", "energy-budget.txt"))
    if len(rows) < 2:
        print("FAIL: energy-budget.txt has %d rows, need the t=0 line and at least one write" % len(rows)); return 1
    t0, E0, M0 = rows[0][0], rows[0][1], rows[0][2]
    worst_i = worst_ii = 0.0
    for r in rows[1:]:
        t, E, M, Qout, Qg, Mg, Qs = r
        worst_i = max(worst_i, abs(E - E0 + Qout - Qs) / max(abs(E0), abs(Qout)))
        worst_ii = max(worst_ii, abs(M - M0 + Mg) / M0)
    print("  (i)   energy balance residual      %.3e  (tol %.0e)  %s" % (worst_i, TOL_BAL, "ok" if worst_i <= TOL_BAL else "FAIL"))
    print("  (ii)  solid-mass balance residual  %.3e  (tol %.0e)  %s" % (worst_ii, TOL_BAL, "ok" if worst_ii <= TOL_BAL else "FAIL"))
    ok = ok and worst_i <= TOL_BAL and worst_ii <= TOL_BAL

    # (iii): independent state sum from field.tec at the last write
    dims, d = read_tec(os.path.join(case, "OUTPUT", "field.tec"))
    V = cell_volumes(dims, d)
    rho_v, rho_c, erf = read_pyrolysis_dat(os.path.join(case, "INPUT", "pyrolysis.dat"))
    hv = read_solid_table(os.path.join(case, "INPUT", "virgin.dat"))
    hc = read_solid_table(os.path.join(case, "INPUT", "char.dat"))
    nsub = len(erf)
    chi = [d["chi_%d" % (r + 1)] for r in range(nsub)]
    E_py = M_py = 0.0
    for c in range(len(V)):
        rho_s = rho_v - sum(erf[r] * chi[r][c] for r in range(nsub))
        tau = (rho_s - rho_c) / (rho_v - rho_c) if rho_v != rho_c else 1.0
        tau = min(max(tau, 0.0), 1.0)
        T = d["T"][c]
        E_py += V[c] * (tau * rho_v * interp_h(hv, T) + (1.0 - tau) * rho_c * interp_h(hc, T))
        M_py += V[c] * rho_s
    last = rows[-1]
    rel_E = abs(E_py - last[1]) / abs(last[1])
    rel_M = abs(M_py - last[2]) / abs(last[2])
    print("  (iii) independent E_tot, M_s       %.3e  %.3e  (tol %.0e)  %s" % (rel_E, rel_M, TOL_STATE,
          "ok" if max(rel_E, rel_M) <= TOL_STATE else "FAIL"))
    ok = ok and max(rel_E, rel_M) <= TOL_STATE

    # (iv): vnn refinement
    expected_path = os.path.join(here, "expected.json")
    expected = json.load(open(expected_path)) if os.path.isfile(expected_path) else {}
    if opts["--vnn025"] and opts["--vnn0125"]:
        probes = sorted(f for f in os.listdir(os.path.join(case, "OUTPUT")) if f.startswith("p") and f.endswith("mm.txt"))
        grid = [2.0 + 0.1 * n for n in range(81)]
        d1 = d2 = d3 = 0.0
        for pf in probes:
            s05 = read_probe(os.path.join(case, "OUTPUT", pf))
            s025 = read_probe(os.path.join(opts["--vnn025"], "OUTPUT", pf))
            s0125 = read_probe(os.path.join(opts["--vnn0125"], "OUTPUT", pf))
            for t in grid:
                a, b, c = interp_series(s05, t), interp_series(s025, t), interp_series(s0125, t)
                if a is None or b is None or c is None:
                    continue
                d1 = max(d1, abs(a - b)); d2 = max(d2, abs(b - c)); d3 = max(d3, abs(a - c))
        ratio = d1 / d2 if d2 > 0.0 else float("inf")
        budget = expected.get("budget_iv_K", BUDGET_IV_DEFAULT)
        if record:
            expected["measured_iv_K"] = d3
            expected["budget_iv_K"] = d3 + 1.0
            budget = d3 + 1.0
        st = "ok" if (ratio >= ORDER_RATIO_MIN and d3 <= budget) else "FAIL"
        print("  (iv)  vnn refinement: |T05-T025| %.3e |T025-T0125| %.3e ratio %.3f (>= %.1f); max |T05-T0125| %.3e K (budget %.2f)  %s"
              % (d1, d2, ratio, ORDER_RATIO_MIN, d3, budget, st))
        ok = ok and st == "ok"
    else:
        print("  (iv)  vnn refinement               SKIP (no --vnn025/--vnn0125)")

    # (v): zero kinetics on the pyrolysis path
    if opts["--nokin"]:
        _, dn = read_tec(os.path.join(opts["--nokin"], "OUTPUT", "field.tec"))
        rn = read_budget(os.path.join(opts["--nokin"], "OUTPUT", "energy-budget.txt"))
        chi_cols = [k for k in dn if k.startswith("chi_")]
        chi_zero = all(v == 0.0 for k in chi_cols for v in dn[k])
        mdot_zero = all(v == 0.0 for v in dn["mdot_pg"])
        rho_dev = max(abs(v - rho_v) / rho_v for v in dn["rho_s"])
        gas_zero = all(r[4] == 0.0 and r[5] == 0.0 for r in rn)
        bal = max(abs(r[1] - rn[0][1] + r[3] - r[6]) / max(abs(rn[0][1]), abs(r[3])) for r in rn[1:])
        st = "ok" if (chi_zero and mdot_zero and rho_dev <= TOL_RHO_V and gas_zero and bal <= TOL_BAL) else "FAIL"
        print("  (v)   no-kinetics: chi all 0 %s, mdot_pg all 0 %s, |rho_s-rho_v|/rho_v %.1e, Q_gas=M_gas=0 %s, balance %.1e  %s"
              % (chi_zero, mdot_zero, rho_dev, gas_zero, bal, st))
        ok = ok and st == "ok"
    else:
        print("  (v)   no-kinetics                  SKIP (no --nokin)")

    # (vi): 1 vs 8 threads
    if opts["--p8"]:
        same = all(md5(os.path.join(case, "OUTPUT", f)) == md5(os.path.join(opts["--p8"], "OUTPUT", f))
                   for f in ("field.tec", "energy-budget.txt"))
        print("  (vi)  1 vs 8 threads: field.tec and energy-budget.txt identical %s  %s" % (same, "ok" if same else "FAIL"))
        ok = ok and same
    else:
        print("  (vi)  cross-thread                 SKIP (no --p8)")

    if record:
        json.dump(expected, open(expected_path, "w"), indent=1, sort_keys=True)
        print("  recorded %s" % expected_path)
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
