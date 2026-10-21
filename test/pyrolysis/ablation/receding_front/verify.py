#!/usr/bin/env python3
"""
P-D3's exact case (plan 12 section D): a receding ablating wall (BC 306) sweeps out the WALL state.

    python3 verify.py FRONT TWIN_PYRO TWIN_LEGACY     # the type-1 runs (generate.py's front, twin_pyro, twin_legacy)
    python3 verify.py --t2 T2                         # the type-2 run (generate.py's t2)

front: a constant-property line (rho_v 280, cp 1000, k 0.4, kinetics frozen) recedes at sdot = 1/700 m/s under an
absorbed flux q = 4e5 W/m2; the continuum steady wall T_w* = T_inf + q/(rho c sdot) = 1300 K, whatever k.

  (G1)  |T_w(last) - 1300| / 1000 <= 2e-3: the closed form (the transient left at t* = 14.9 is +-1e-3; a cell-state
        sweep gives +0.20)
  (G2)  every step n: E_abl(n) - E_abl(n-1) = sum_faces A sdot (t_n - t_(n-1)) rho_v h_v(T_w^(n-1)), T_w^(n-1) the face's
        seb.txt row of step n-1 (step 1: the surface cell of ic.tec, since a new run's Tw_prev is 0); h_v = 1000 (T -
        298.15): relative 1e-12 (a cell-state sweep misses by 16 %)
  (G3)  E_tot + E_abl + Q_out - Q_src = E_tot(0) on every row: 1e-12 of max(|E_tot(0)|, max |Q_out|)
  (G4)  M_abl = rho_v A_tot sdot t on every row (1e-13 relative); s_total of bc306-state.dat = sdot t_end (1e-13); the
        surface nodes of field.tec at x = L - s_total (1e-15 m)
  (G5)  the discrete oracle: a numpy replica of the scheme (start-of-step dt, the s/L node move of nodeshift, the
        donor-cell remap sweeping Tw_prev at the 306 face, SSP-RK3 of the two-point conduction with T_w = T_c + q/G,
        G = 2k/dx): FUSS's T_w at every step and T in every cell of the last field, max |diff| / 1000 <= 1e-12
        (measured 4.3e-14 and 1.3e-13 on two builds, before and after upstream's diffusive-metric default; a bound
        of 1.2 times the first measurement would break on a code-generation move)
  (G6)  printed: the last profile against the continuum curve, the cell-to-cell decay against 1/(1 + Pe)
  (twin) twin_pyro (Robin hconv 400 W/m2/K to 2300 K: T_w* 1300 K again, the pyrolysis remap) |err| <= 2e-3, and
        twin_legacy (the same without [FUSS-Pyrolysis]: the legacy remap sweeps the cell's own h) err >= 5e-2: the
        gate shows on every run that it tells the two sweeps apart

t2 (type 2, K = 1e-12 m2, p_w = 101325 Pa, 20 steps with a field file each): for every step n >= 2 (step 1 sweeps the
set-up state, which no file holds), from the surface cell's eta_st of field file n-1 and the face's T_w of step n-1:
  (T2-a) E_abl(n) - E_abl(n-1) = sum_faces A delta_n [rho_v h_v(T_w) + eta_st p_w h_g(T_w) - eps_g p_w], 1e-12
  (T2-b) G_abl(n) - G_abl(n-1) = sum_faces A delta_n eta_st p_w, 1e-12: the stored gas leaves at p_w, not at p_c
         (a sweep at p_c misses by (p_c - p_w)/p_w, 2e-4 here)
  (bal-m) M_s + M_g + M_gas_out + M_abl + G_abl constant, 1e-12 of itself; (bal-e) as (G3); p_c - p_w printed.
  h_g = 1000 (T - 298.15) and eps_g = 0.8: the --isothermal gas table and pyrolysis.dat of ../../exact/tables.

The per-step identities (G2, T2-a, T2-b) stop at about 3e-13: FUSS's swept volume is a difference of node positions
near L = 0.1 m for a step of about 6e-5 m (2e-16 x 0.1 / 6e-5); hence 1e-12, against 0.16 (G2) and 2e-4 (T2-b) for
the wrong state.

The gate's selftest raises the flux by 1 % in front's INPUT (env_front.dat): (G1) and (G5) must fail, the rest pass.
Python 3 + numpy.
"""
import glob
import math
import os
import re
import sys

import numpy as np

NX, L, NY, DY, DZ = 200, 0.1, 2, 0.1, 0.1
A_FACE = DY * DZ
A_TOT = NY * A_FACE
K_COND, RHO_V, CP = 0.4, 280.0, 1000.0
SDOT = 1.0 / 700.0
Q, TINF, TWSTAR, DT_REF = 4.0e5, 300.0, 1300.0, 1000.0
PW, EPS_G, H0 = 101325.0, 0.8, 298.15
TOL_G1, TOL_G2, TOL_G3, TOL_G4, TOL_S, TOL_X = 2.0e-3, 1.0e-12, 1.0e-12, 1.0e-13, 1.0e-13, 1.0e-15
TOL_G5 = 1.0e-12
TOL_TWIN_PYRO, MIN_TWIN_LEGACY = 2.0e-3, 5.0e-2
TOL_T2A, TOL_T2B, TOL_BALM, TOL_BALE = 1.0e-12, 1.0e-12, 1.0e-12, 1.0e-12

bad = []


def check(tag, ok, text):
    print("%-8s %s  %s" % ("(" + tag + ")", "ok  " if ok else "FAIL", text))
    if not ok:
        bad.append(tag)


def h_v(T):
    return CP * (T - H0)


def read_cols(path):
    """A FUSS text output with a `# name name ...` header: {name: array}."""
    names, rows = None, []
    for line in open(path):
        if line.startswith("#"):
            if names is None:
                names = line[1:].split("(")[0].split()
            continue
        if line.strip():
            rows.append([float(v) for v in line.split()])
    a = np.array(rows)
    return {n: a[:, k] for k, n in enumerate(names)}


def read_seb(path):
    """seb.txt: {face j: (times, T_w)} for the 306 faces (two rows a step, j = 1, 2)."""
    out = {}
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        v = line.split()
        out.setdefault(int(v[3]), []).append((float(v[0]), float(v[6])))
    return {j: (np.array([r[0] for r in rows]), np.array([r[1] for r in rows])) for j, rows in out.items()}


def read_field(path):
    """Every zone of a Tecplot BLOCK file: [{name: array}] (cells and nodes in file order, i fastest)."""
    txt = open(path).read()
    names = txt.split("\n", 1)[0].split("=", 1)[1].replace('"', " ").split()
    zones = []
    parts = re.split(r"\n\s*ZONE", txt)
    for part in parts[1:]:
        head, body = part.split("\n", 1)
        I, J, K = (int(v) for v in re.search(r"I=(\d+), J=(\d+), K=(\d+)", head).groups())
        vals = np.array(body.split(), float)
        nn, nc = I * J * K, (I - 1) * (J - 1) * max(K - 1, 1)
        z, off = {"_dims": (I, J, K)}, 0
        for n, nm in enumerate(names):
            size = nn if n < 3 else nc
            z[nm] = vals[off:off + size]
            off += size
        zones.append(z)
    return zones


def ic_surface_T(run):
    """The surface cell (i = NX) of each line j in INPUT/ic.tec: step 1 sweeps it (Tw_prev = 0 on a new run)."""
    lines = open(os.path.join(run, "INPUT", "ic.tec")).read().split("\n")
    I, J, K = (int(v) for v in re.search(r"I=(\d+), J=(\d+), K=(\d+)", lines[1]).groups())
    first = 2 + 3 * I * J * K
    T = [float(lines[first + n]) for n in range((I - 1) * (J - 1) * max(K - 1, 1))]
    return {j: T[(j - 1) * NX + NX - 1] for j in range(1, NY + 1)}


def replica(nsteps, T_init, q=Q):
    """FUSS's 1-D receding step in numpy (plan 12 research R3, replica_np.py), the wall sweep, constant flux.
    Cell 0 is the surface cell; returns the per-step T_w (stage 3, Tw_prev) and the last T."""
    k, rc = K_COND, RHO_V * CP
    a = k / rc
    x = L * (1.0 - np.arange(NX + 1) / NX)
    T = np.array(T_init, float)
    twp, tws = -1.0, []
    for _ in range(nsteps):
        dx = x[:-1] - x[1:]
        dt = np.min(0.5 * 0.5 / (a * (1.0 / dx ** 2 + 1.0 / DY ** 2)))
        delta = SDOT * dt
        xo = x.copy()
        d = delta * (xo / xo[0])
        x = xo - d
        vold, vnew = xo[:-1] - xo[1:], x[:-1] - x[1:]
        Eo = rc * T
        acc = Eo * vold
        eup = Eo.copy()
        eup[0] = rc * (twp if twp > 0 else T[0])
        acc -= eup * d[:-1]
        acc[:-1] += Eo[1:] * d[1:-1]
        T = (acc / vnew) / rc
        dx = vnew
        G = 2.0 * k / dx[0]
        dist = 0.5 * (dx[:-1] + dx[1:])

        def resid(TT):
            F = k * (TT[1:] - TT[:-1]) / dist
            R = np.zeros(NX)
            R[:-1] += F
            R[1:] -= F
            R[0] += q
            return R, TT[0] + q / G
        E0 = rc * T
        R0, _ = resid(T)
        u1 = E0 + dt * R0 / dx
        R1, _ = resid(u1 / rc)
        u2 = 0.75 * E0 + 0.25 * (u1 + dt * R1 / dx)
        R2, tw3 = resid(u2 / rc)
        T = (E0 / 3.0 + (2.0 / 3.0) * (u2 + dt * R2 / dx)) / rc
        twp = tw3
        tws.append(twp)
    return np.array(tws), T


def last_Tw(run):
    st = [l.split() for l in open(os.path.join(run, "OUTPUT", "bc306-state.dat")) if not l.startswith("#")]
    return np.array([float(v[6]) for v in st])


def type1(front, twin_pyro, twin_legacy):
    eb = read_cols(os.path.join(front, "OUTPUT", "energy-budget.txt"))
    seb = read_seb(os.path.join(front, "OUTPUT", "seb.txt"))
    t = eb["t"]
    nsteps = len(t) - 1
    tw_last = np.array([seb[j][1][-1] for j in sorted(seb)])
    err = (tw_last - TWSTAR) / DT_REF
    check("G1", np.all(np.abs(err) <= TOL_G1),
          "T_w(t = %.6f s, step %d) = %s K: (T_w - 1300)/1000 = %s (|.| <= %.0e)"
          % (t[-1], nsteps, ", ".join("%.9f" % v for v in tw_last), ", ".join("%+.2e" % v for v in err), TOL_G1))

    tw0 = ic_surface_T(front)
    dE = np.diff(eb["E_abl"])
    worst2 = 0.0
    for n in range(1, nsteps + 1):
        delta = SDOT * (t[n] - t[n - 1])
        pred = 0.0
        for j in sorted(seb):
            tw = tw0[j] if n == 1 else seb[j][1][n - 2]
            pred += A_FACE * delta * RHO_V * h_v(tw)
        worst2 = max(worst2, abs(dE[n - 1] - pred) / abs(pred))
    check("G2", worst2 <= TOL_G2, "swept energy against rho_v h_v(T_w^(n-1)) A sdot dt, %d steps: worst %.2e (<= %.0e)"
          % (nsteps, worst2, TOL_G2))

    bal = eb["E_tot"] + eb["E_abl"] + eb["Q_out"] - eb["Q_src"]
    scale = max(abs(eb["E_tot"][0]), np.max(np.abs(eb["Q_out"])))
    g3 = np.max(np.abs(bal - bal[0])) / scale
    check("G3", g3 <= TOL_G3, "E_tot + E_abl + Q_out - Q_src drift %.2e of %.4e J (<= %.0e)" % (g3, scale, TOL_G3))

    mab = RHO_V * A_TOT * SDOT * t
    g4 = np.max(np.abs(eb["M_abl"][1:] - mab[1:]) / mab[1:])
    check("G4", g4 <= TOL_G4, "M_abl = rho_v A_tot sdot t, %d rows: worst %.2e (<= %.0e)" % (nsteps, g4, TOL_G4))
    st = [l.split() for l in open(os.path.join(front, "OUTPUT", "bc306-state.dat")) if not l.startswith("#")]
    s_tot = np.array([float(v[5]) for v in st])
    gs = np.max(np.abs(s_tot - SDOT * t[-1])) / (SDOT * t[-1])
    check("G4", gs <= TOL_S, "s_total = sdot t_end = %.15e m: worst %.2e (<= %.0e)" % (SDOT * t[-1], gs, TOL_S))
    z = read_field(os.path.join(front, "OUTPUT", "field.tec"))[0]
    I, J, K = z["_dims"]
    xs = z["x"].reshape(K, J, I)[:, :, -1]
    gx = np.max(np.abs(xs - (L - s_tot[0])))
    check("G4", gx <= TOL_X, "surface nodes at x = L - s_total: worst %.2e m (<= %.0e)" % (gx, TOL_X))

    T_f = z["T"].reshape(1, NY, NX)[0]
    ic = read_field(os.path.join(front, "INPUT", "ic.tec"))[0]["T"].reshape(1, NY, NX)[0]
    tws, T_r = replica(nsteps, ic[0][::-1])
    g5w = max(np.max(np.abs(seb[j][1] - tws)) for j in sorted(seb)) / DT_REF
    g5t = max(np.max(np.abs(T_f[j][::-1] - T_r)) for j in range(NY)) / DT_REF
    check("G5", max(g5w, g5t) <= TOL_G5, "the numpy replica: T_w every step %.2e, T every cell of the last field %.2e "
          "(of 1000 K; <= %.0e)" % (g5w, g5t, TOL_G5))

    xc = 0.5 * (z["x"].reshape(K, J, I)[0, 0, :-1] + z["x"].reshape(K, J, I)[0, 0, 1:])
    xw = L - s_tot[0]
    land = TINF + (TWSTAR - TINF) * np.exp(-(xw - xc) * SDOT * RHO_V * CP / K_COND)
    near = (xw - xc) <= 5.0e-3
    pe = SDOT * (xc[1] - xc[0]) * RHO_V * CP / K_COND
    r = (T_f[0][-2:-12:-1] - TINF) / (T_f[0][-1:-11:-1] - TINF)
    print("(G6)     printed: last profile vs the continuum curve within 5 mm of the wall: max |dT| %.3f K; decay ratio of "
          "the 10 cells next to the wall %.4f .. %.4f (semi-discrete upwind 1/(1 + Pe) = %.4f, Pe %.3f)"
          % (np.max(np.abs(T_f[0][near] - land[near])), r.min(), r.max(), 1.0 / (1.0 + pe), pe))

    ep = (last_Tw(twin_pyro) - TWSTAR) / DT_REF
    el = (last_Tw(twin_legacy) - TWSTAR) / DT_REF
    check("twin", np.all(np.abs(ep) <= TOL_TWIN_PYRO),
          "twin_pyro (Robin wall, pyrolysis remap): (T_w - 1300)/1000 = %s (|.| <= %.0e)"
          % (", ".join("%+.2e" % v for v in ep), TOL_TWIN_PYRO))
    check("twin", np.all(el >= MIN_TWIN_LEGACY),
          "twin_legacy (Robin wall, legacy remap: the cell state): (T_w - 1300)/1000 = %s (>= %.0e)"
          % (", ".join("%+.2e" % v for v in el), MIN_TWIN_LEGACY))


def type2(run):
    eb = read_cols(os.path.join(run, "OUTPUT", "energy-budget.txt"))
    seb = read_seb(os.path.join(run, "OUTPUT", "seb.txt"))
    t = eb["t"]
    nsteps = len(t) - 1
    files = {}
    for f in glob.glob(os.path.join(run, "OUTPUT", "field*.tec")):
        m = re.search(r"field(\d*)\.tec$", f)
        files[int(m.group(1)) if m.group(1) else nsteps] = f
    if sorted(files) != list(range(1, nsteps + 1)):
        check("T2", False, "field files for steps 1..%d expected, found %s" % (nsteps, sorted(files)))
        return
    eta = {}
    pc = {}
    for n, f in files.items():
        z = read_field(f)[0]
        eta[n] = z["eta_st"].reshape(NY, NX)[:, -1]
        pc[n] = z["p_pg"].reshape(NY, NX)[:, -1]
    dE, dG = np.diff(eb["E_abl"]), np.diff(eb["G_abl"])
    wa = wb = 0.0
    for n in range(2, nsteps + 1):
        delta = SDOT * (t[n] - t[n - 1])
        pe_ = pg_ = 0.0
        for jj, j in enumerate(sorted(seb)):
            tw = seb[j][1][n - 2]
            g = eta[n - 1][jj] * PW
            pe_ += A_FACE * delta * (RHO_V * h_v(tw) + g * h_v(tw) - EPS_G * PW)
            pg_ += A_FACE * delta * g
        wa = max(wa, abs(dE[n - 1] - pe_) / abs(pe_))
        wb = max(wb, abs(dG[n - 1] - pg_) / abs(pg_))
    check("T2-a", wa <= TOL_T2A, "swept energy at (T_w, p_w), steps 2..%d: worst %.2e (<= %.0e)" % (nsteps, wa, TOL_T2A))
    check("T2-b", wb <= TOL_T2B, "swept gas eta_st p_w, steps 2..%d: worst %.2e (<= %.0e)" % (nsteps, wb, TOL_T2B))
    m = eb["M_s"] + eb["M_g"] + eb["M_gas_out"] + eb["M_abl"] + eb["G_abl"]
    bm = np.max(np.abs(m - m[0])) / abs(m[0])
    check("bal-m", bm <= TOL_BALM, "M_s + M_g + M_gas_out + M_abl + G_abl drift %.2e (<= %.0e)" % (bm, TOL_BALM))
    bal = eb["E_tot"] + eb["E_abl"] + eb["Q_out"] - eb["Q_src"]
    scale = max(abs(eb["E_tot"][0]), np.max(np.abs(eb["Q_out"])))
    be = np.max(np.abs(bal - bal[0])) / scale
    check("bal-e", be <= TOL_BALE, "E_tot + E_abl + Q_out - Q_src drift %.2e of %.4e J (<= %.0e)" % (be, scale, TOL_BALE))
    print("(p)      printed: p_c - p_w at step %d = %s Pa (%s of p_w): what a sweep at p_c would miss"
          % (nsteps, ", ".join("%.3f" % (v - PW) for v in pc[nsteps]),
             ", ".join("%.2e" % ((v - PW) / PW) for v in pc[nsteps])))


def main():
    args = sys.argv[1:]
    if args[:1] == ["--t2"] and len(args) == 2:
        type2(args[1])
    elif len(args) == 3:
        type1(*args)
    else:
        sys.exit(__doc__)
    print("FAIL: %s" % " ".join(sorted(set(bad))) if bad else "PASS")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
