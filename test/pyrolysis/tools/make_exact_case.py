#!/usr/bin/env python3
"""
make_exact_case.py -- write a FUSS case whose exact solution is in darcy_exact.py.

Exact-solution benchmarks of the type-2 (Darcy) closure (2026-10-10), the three problems of
darcy_exact.py, on a 1-D slab of L = 0.05 m (nx cells along i, 2 along j, 1 along k, as every
1-D pyrolysis case), with the --isothermal table set of make_tacot_tables.py (solid cp 1000,
k 0.4, rho_v 280; gas M 0.022, mu 3e-5, h = 1000 (T - 298.15), eps_g 0.8):

  pressure_transient  sealed at x = 0 (301 q 0), open at x = L (303 hconv 0, p_w = p0 + dp),
                      kinetics off, T0 = 298.15 K (h_g = 0), p0 = 101325 Pa; K chosen so that the
                      pressure diffusivity D = K p0/(eps mu) is --diff-ratio times alpha = k/(rho c);
                      run to --tend-frac of L^2/D with --nout snapshots. Small dp: the linear limit.
  transpiration       both faces open (303, per-face p_w): p_in at x = 0, p_out at x = L, Robin heat
                      exchange (h1, Tref1) and (h2, Tref2); kinetics off; T(0) = the exact steady
                      profile, run for --tend (default 20 slow e-folding times) to the discrete steady state.
  sharp_front         semi-infinite virgin slab at T0, face 1 held at Ts (302, open to the gas), face
                      2 sealed (301 q 0); ONE reaction, phase 2, F = 0.5 (rho_c = 220), E = 0, A = --A,
                      order 1, threshold T = Tp: the solid decomposes at once where T >= Tp. --closure
                      type2 (K_c, K_v, open-bc-types 302) or type1 (outflow-direction -i).

Usage: make_exact_case.py <case-dir> --kind KIND --tables <iso-dir> [options, see --help]
Python 3.6, no third-party modules; imports make_case.py, make_darcy_case.py, darcy_exact.py.
"""
import argparse
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_case as mc          # noqa: E402
import make_darcy_case as mdc   # noqa: E402
import darcy_exact as dx        # noqa: E402

L = 0.05
P0 = 101325.0
NY, DY, DZ = 2, 0.1, 0.1


def die(msg):
    sys.stderr.write("[ERROR] make_exact_case: %s\n" % msg)
    sys.exit(1)


def read_key(path, key):
    for s in open(path):
        tok = s.split("//", 1)[0].split()
        if len(tok) >= 2 and tok[0] == key:
            return tok[1]
    die("%s has no key %s" % (path, key))


def gas_consts(path):
    """M, mu and cp_g (from two h rows) of the first zone of an isothermal gas.dat."""
    rows = []
    for s in open(path):
        tok = s.split()
        if len(tok) == 5:
            try:
                rows.append([float(v) for v in tok])
            except ValueError:
                pass
        if len(rows) == 2:
            break
    (T1, _, M, h1, mu), (T2, _, _, h2, _) = rows
    return M, mu, (h2 - h1) / (T2 - T1)


def write_ic_T(path, xs, ys, zs, Tcells):
    """ic.tec with a per-cell T along i (the same in every j, k)."""
    ni, nj, nk = len(xs), len(ys), len(zs)
    out = [' VARIABLES ="x" "y" "z" T matID',
           " ZONE  T = B1-SP, I=%d, J=%d, K=%d, DATAPACKING=BLOCK, VARLOCATION=([1-3]=NODAL,[4-5]=CELLCENTERED),"
           " SOLUTIONTIME=0.000000000000000E+000, STRANDID = 0" % (ni, nj, nk)]
    for axis in range(3):
        for k in range(nk):
            for j in range(nj):
                for i in range(ni):
                    out.append(" %.15E" % (xs[i], ys[j], zs[k])[axis])
    for k in range(nk - 1):
        for j in range(nj - 1):
            for i in range(ni - 1):
                out.append(" %.15E" % Tcells[i])
    for _ in range((ni - 1) * (nj - 1) * (nk - 1)):
        out.append(" %.15E" % 1.0)
    with open(path, "w") as fh:
        fh.write("\n".join(out) + "\n")


def freeze_solid(path, factor):
    """Scale cp, h and k of every row of a solid table by `factor`: alpha = k/(rho cp) is unchanged (so the VNN time
    step is), while rho cp grows by `factor` and the temperature no longer responds to the pressure work eps dp."""
    out = []
    for s in open(path):
        tok = s.split()
        if len(tok) == 8:
            try:
                v = [float(t) for t in tok]
            except ValueError:
                out.append(s)
                continue
            v[1] *= factor; v[2] *= factor; v[3] *= factor; v[4] *= factor; v[5] *= factor
            out.append(" ".join("%.9E" % x if i else "%.1f" % x for i, x in enumerate(v)) + "\n")
        else:
            out.append(s)
    with open(path, "w") as fh:
        fh.write("".join(out))


def single_reaction(path, Tp, A, F):
    """Replace every reaction of pyrolysis.dat by ONE reaction of phase 2: F, A, E = 0, order 1, threshold Tp."""
    out = []
    for s in open(path):
        tok = s.split("//", 1)[0].split()
        k = tok[0] if tok else ""
        if k.startswith("nPyroReac["):
            out.append("%s %d\n" % (k, 1 if k == "nPyroReac[2]" else 0))
        elif k[:2] in ("F[", "A[", "E[", "m[", "n[", "T[", "h[") and "][" in k:
            continue
        else:
            out.append(s)
    out += ["F[2][1] %.15g\n" % F, "A[2][1] %.15g\n" % A, "E[2][1] 0\n", "m[2][1] 1\n", "n[2][1] 0\n",
            "T[2][1] %.15g\n" % Tp, "h[2][1] 0\n"]
    with open(path, "w") as fh:
        fh.write("".join(out))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case_dir")
    ap.add_argument("--kind", required=True, choices=["pressure_transient", "transpiration", "sharp_front"])
    ap.add_argument("--tables", required=True)
    ap.add_argument("--ncells", type=int, default=40)
    ap.add_argument("--vnn", type=float, default=0.5)
    ap.add_argument("--fuss-sh", default=mc.DONOR_FUSS_SH)
    ap.add_argument("--masterdir", default=None, help="absolute MASTERDIR for FUSS.sh (a case outside the tree)")
    # pressure_transient
    ap.add_argument("--dp", type=float, default=1.0, help="p_w - p0 [Pa]")
    ap.add_argument("--diff-ratio", type=float, default=1.0, help="D / alpha")
    ap.add_argument("--tend-frac", type=float, default=0.5, help="t_end / (L^2/D)")
    ap.add_argument("--nout", type=int, default=25)
    ap.add_argument("--freeze-factor", type=float, default=1.0e12,
                    help="pressure_transient: solid cp, h and k times this (alpha kept): T frozen against the pressure "
                         "work eps dp (the isothermal series ignores it; 1 = the coupled model)")
    # transpiration
    ap.add_argument("--K", type=float, default=4.7e-13)
    ap.add_argument("--pin", type=float, default=2.0 * P0)
    ap.add_argument("--pout", type=float, default=P0)
    ap.add_argument("--h1", type=float, default=50.0)
    ap.add_argument("--T1", type=float, default=300.0)
    ap.add_argument("--h2", type=float, default=50.0)
    ap.add_argument("--T2", type=float, default=1000.0)
    ap.add_argument("--tend", type=float, default=None)
    # sharp_front
    ap.add_argument("--closure", choices=["type1", "type2"], default="type2")
    ap.add_argument("--T0", type=float, default=300.0)
    ap.add_argument("--Ts", type=float, default=1200.0)
    ap.add_argument("--Tp", type=float, default=600.0)
    ap.add_argument("--A", type=float, default=1.0e9)
    ap.add_argument("--Kc", type=float, default=1.0e-11)
    ap.add_argument("--Kv", type=float, default=1.0e-18)
    args = ap.parse_args()

    for f in mc.TABLE_FILES:
        if not os.path.isfile(os.path.join(args.tables, f)):
            die("%s not in --tables %s" % (f, args.tables))
    case = os.path.abspath(args.case_dir)
    for d in ("INPUT", "MESH"):
        os.makedirs(os.path.join(case, d), exist_ok=True)
    for f in mc.TABLE_FILES:
        mc.shutil.copy(os.path.join(args.tables, f), os.path.join(case, "INPUT", f))
    pyro = os.path.join(case, "INPUT", "pyrolysis.dat")
    name, rho_v, _ = mc.read_pyrolysis_dat(pyro)
    virgin = mc.read_table_zone1(os.path.join(case, "INPUT", "virgin.dat"), 7)
    eps = float(read_key(pyro, "eps_g_v"))
    if abs(eps - float(read_key(pyro, "eps_g_c"))) > 0:
        die("the exact problems need eps_g_v = eps_g_c (the --isothermal set)")
    k_s, cp_s = 0.4, 1000.0
    W, mu, cpg = gas_consts(os.path.join(case, "INPUT", "gas.dat"))
    alpha = k_s / (rho_v * cp_s)

    nx = args.ncells
    xs = [L * i / nx for i in range(nx + 1)]
    xs[-1] = L
    xc = [0.5 * (xs[i] + xs[i + 1]) for i in range(nx)]
    ys = [DY * j for j in range(NY + 1)]
    zs = [0.0, DZ]
    kind = args.kind
    meta = {"kind": kind, "L": L, "ncells": nx, "vnn": args.vnn, "p0": P0, "eps": eps, "W": W, "mu": mu,
            "cpg": cpg, "k": k_s, "cp": cp_s, "rho_v": rho_v, "alpha": alpha}
    darcy = None
    pyro_sec = ["[FUSS-Pyrolysis]", "gas-closure = type2"]
    sealed = (301, mdc.VAL_FMT % 0.0)

    if kind == "pressure_transient":
        mc.strip_kinetics(pyro)
        if args.freeze_factor != 1.0:
            for f in ("virgin.dat", "char.dat"):
                freeze_solid(os.path.join(case, "INPUT", f), args.freeze_factor)
        meta["freeze_factor"] = args.freeze_factor
        D = args.diff_ratio * alpha
        K = D * eps * mu / P0
        mdc.write_permeability(pyro, [K, 0, 0, 0, K, 0, 0, 0, K], [K, 0, 0, 0, K, 0, 0, 0, K])
        T0 = 298.15
        tau = L * L / D
        tend = args.tend_frac * tau
        pw = P0 + args.dp
        face2 = (303, mdc.VAL4_FMT % (0.0, 0.0, T0, pw))
        spec = {1: lambda i, j, k: sealed, 2: lambda i, j, k: face2}
        Tcells = [T0] * nx
        darcy = ["[FUSS-Darcy]", "solver = tridiag", "tol = 1e-12", "max-iter = 100000",
                 "p-init = %.1f" % P0, "p-wall = %.1f" % P0, "quasi-steady = false"]
        timing = ["res-threshold = 1e-30", "time-threshold = %.15g" % tend, "iter-threshold = 100000000"]
        sol_dt = tend / args.nout
        meta.update({"T0": T0, "pw": pw, "dp": args.dp, "D": D, "K": K, "tau_p": tau, "t_end": tend})
        head = ["Transient gas pressure in a sealed/open porous slab, small amplitude (exact: darcy_exact.py 1)",
                "%d cells, K %.6g m2, D = K p0/(eps mu) = %.6g m2/s = %g alpha, p_w - p0 = %g Pa," % (nx, K, D, args.diff_ratio, args.dp),
                "T0 %.2f K (h_g = 0), to %.6g s = %g L^2/D, vnn %g" % (T0, tend, args.tend_frac, args.vnn)]
    elif kind == "transpiration":
        mc.strip_kinetics(pyro)
        K = args.K
        mdc.write_permeability(pyro, [K, 0, 0, 0, K, 0, 0, 0, K], [K, 0, 0, 0, K, 0, 0, 0, K])
        md, A, B, Pe, lin = dx.transpiration_steady(L, k_s, cpg, W, mu, K, args.pin, args.pout,
                                                    args.h1, args.T1, args.h2, args.T2)
        f1 = (303, mdc.VAL4_FMT % (args.h1, 0.0, args.T1, args.pin))
        f2 = (303, mdc.VAL4_FMT % (args.h2, 0.0, args.T2, args.pout))
        spec = {1: lambda i, j, k: f1, 2: lambda i, j, k: f2}
        Tcells = [dx.transpiration_T(x, L, A, B, Pe, lin) for x in xc]
        tau_slow = L * L / (math.pi ** 2 * alpha)
        tend = args.tend if args.tend is not None else 20.0 * tau_slow
        darcy = ["[FUSS-Darcy]", "solver = tridiag", "tol = 1e-12", "max-iter = 100000",
                 "p-init = %.15g" % (0.5 * (args.pin + args.pout)), "p-wall = %.1f" % P0, "quasi-steady = false"]
        timing = ["res-threshold = 1e-30", "time-threshold = %.15g" % tend, "iter-threshold = 100000000"]
        sol_dt = tend / 4.0
        meta.update({"K": K, "pin": args.pin, "pout": args.pout, "h1": args.h1, "T1": args.T1, "h2": args.h2,
                     "T2": args.T2, "mdot": md, "Pe": Pe, "t_end": tend})
        head = ["Steady transpiration through a porous slab (exact: darcy_exact.py 2)",
                "%d cells, K %.6g m2, p_in %.6g Pa at x = 0, p_out %.6g Pa at x = L, mdot %.6g kg/m2/s, Pe %.4f," % (nx, K, args.pin, args.pout, md, Pe),
                "Robin (h %g, Tref %g) at x = 0 and (h %g, Tref %g) at x = L; T(0) = the exact profile, to %.6g s" % (args.h1, args.T1, args.h2, args.T2, tend)]
    else:
        F = 0.5
        single_reaction(pyro, args.Tp, args.A, F)
        rho_c = rho_v - F * 0.1 * 1200.0
        Kc, Kv = args.Kc, args.Kv
        mdc.write_permeability(pyro, [Kv, 0, 0, 0, Kv, 0, 0, 0, Kv], [Kc, 0, 0, 0, Kc, 0, 0, 0, Kc])
        face1 = (302, mdc.VAL_FMT % args.Ts)
        spec = {1: lambda i, j, k: face1, 2: lambda i, j, k: sealed}
        Tcells = [args.T0] * nx
        tend = args.tend if args.tend is not None else 30.0
        if args.closure == "type2":
            darcy = ["[FUSS-Darcy]", "solver = tridiag", "tol = 1e-10", "max-iter = 100000",
                     "p-init = %.1f" % P0, "p-wall = %.1f" % P0, "open-bc-types = 302", "quasi-steady = false"]
        else:
            pyro_sec = ["[FUSS-Pyrolysis]", "gas-closure = type1", "outflow-direction = -i", "p-gas = %.1f" % P0]
        timing = ["res-threshold = 1e-30", "time-threshold = %.15g" % tend, "iter-threshold = 100000000"]
        sol_dt = tend / 60.0
        ac, av = k_s / (rho_c * cp_s), alpha
        lam = dx.sharp_front_lambda(ac, av, k_s, k_s, rho_v, rho_c, cp_s, cpg, args.Ts, args.Tp, args.T0, 0.0)
        meta.update({"closure": args.closure, "T0": args.T0, "Ts": args.Ts, "Tp": args.Tp, "A": args.A, "F": F,
                     "rho_c": rho_c, "Kc": Kc, "Kv": Kv, "alpha_c": ac, "alpha_v": av, "lam": lam, "Q": 0.0,
                     "t_end": tend})
        head = ["Decomposition front with gas transpiration through the char (exact: darcy_exact.py 3)",
                "%d cells, T0 %g K, surface %g K (302) from t = 0, instant decomposition at Tp %g K (A %g, E 0)," % (nx, args.T0, args.Ts, args.Tp, args.A),
                "rho_v %g -> rho_c %g, closure %s, K_c %g, K_v %g, lam %.8f, to %g s" % (rho_v, rho_c, args.closure, Kc, Kv, lam, tend)]

    write_ic_T(os.path.join(case, "INPUT", "ic.tec"), xs, ys, zs, Tcells)
    mdc.write_ic(os.path.join(case, "MESH", "mesh.tec"), xs, ys, zs, "Block1", False, 0.0)
    spec.update({3: lambda i, j, k: (300, None), 4: lambda i, j, k: (300, None),
                 5: lambda i, j, k: (0, None), 6: lambda i, j, k: (0, None)})
    mdc.write_bc(os.path.join(case, "INPUT", "bc.txt"), nx, NY, spec)
    with open(os.path.join(case, "INPUT", "phase.txt"), "w") as fh:
        fh.write("solid phase\n%s 1\n" % name)
    mc.write_properties_dat(os.path.join(case, "INPUT", "properties.dat"), name, virgin, rho_v)

    lines = ["; " + "=" * 77] + [";  " + s for s in head] + [";", ";  Generated by test/pyrolysis/tools/make_exact_case.py -- edit that, not this.", "; " + "=" * 77]
    lines += ["", "[FUSS-Parameters]", "newrun = true"] + timing
    lines += ["", "[FUSS-Numerics]", "time-scheme = RK3", "vnn = %g" % args.vnn, "time-accurate = true",
              "integration-variables = cons", "irs = false", ""] + pyro_sec
    if darcy:
        lines += [""] + darcy
    lines += ["", "[FUSS-IO]", "sol-dtime     = %.15g" % sol_dt, "sol-overwrite = false", "sol-format    = tecplot ascii",
              "ini-diter     = 100000000"]
    with open(os.path.join(case, "input.ini"), "w") as fh:
        fh.write("\n".join(lines) + "\n")

    if args.masterdir:
        with open(args.fuss_sh) as fh:
            src = fh.read().split("\n")
        out = []
        for s in src:
            out.append("MASTERDIR=%s" % args.masterdir if s.startswith("MASTERDIR=") else s)
        with open(os.path.join(case, "FUSS.sh"), "w") as fh:
            fh.write("\n".join(out))
        os.chmod(os.path.join(case, "FUSS.sh"), 0o755)
    else:
        depth, _ = mc.case_depth(case, "adiabatic0d", 1000.0, None)
        mc.write_fuss_sh(os.path.join(case, "FUSS.sh"), args.fuss_sh, depth)
    with open(os.path.join(case, "case.json"), "w") as fh:
        json.dump(meta, fh, indent=1, sort_keys=True)
        fh.write("\n")
    print("make_exact_case.py: %s -> %s  (%s)" % (kind, case, "; ".join(h.strip() for h in head[1:])))


if __name__ == "__main__":
    main()
