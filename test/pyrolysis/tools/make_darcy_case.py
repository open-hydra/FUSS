#!/usr/bin/env python3
"""
make_darcy_case.py -- write a complete FUSS Darcy (type-2 gas closure) verification case.

Implements plan 07 section 4: the exact-solution cases of the pressure solver, in the standard
test-case layout of make_case.py (FUSS.sh, input.ini, MESH/mesh.tec, INPUT/{ic.tec, bc.txt,
phase.txt, properties.dat, pyrolysis.dat, virgin.dat, char.dat, gas.dat}). The kinetics are ALWAYS
stripped (every nPyroReac 0): the gas source is [FUSS-Darcy] pi-source, a constant pi_tot.

Usage:
    make_darcy_case.py <case-dir> --kind {sealed_nosource,mass_balance_1d,sealed_profile,tensor_linear_2d}
                       --tables <dir> [--ncells N] [--quasi-steady] [--K K] [--tensor45 K1 K2] [--identity]
                       [--ny N] [--depth N] [--fuss-sh <path>]

    sealed_nosource   gate 4.1(a): 1-D slab, N cells (default 100) x L/N, L = 0.05 m, both i-faces
                      301 q = 0 (impermeable), pi-source = 0, T(0) = 1000 K, up to 200 material
                      steps -- the run ends at its FIRST step, because nothing moves and the
                      convergence residual (T, chi and p increments) is exactly zero; p must stay
                      at p-init to 1e-14 and every face flux at 0.
    mass_balance_1d   gate 4.1(b): face 1 sealed (301 q = 0), face 2 open (303 hconv = 0, qrad = 0,
                      Tref = 1000, p_w = 101325), pi-source = 0.1 kg/m3/s, K = 1e-11 m2 isotropic,
                      T(0) = 1000 K, 2 s time-accurate; --quasi-steady for (b1), without for (b2).
                      Steady: mdot(x_f) = pi x_f at every face.
    sealed_profile    gate 4.2: as mass_balance_1d with K = 1e-13 m2, quasi-steady, tol 1e-13, the
                      solid cp of the tables 1e12 (T frozen), run as material steps = Picard
                      iterations on the lagged Gamma until p is stationary (400 steps or the
                      residual stop); exact p^2 = p_w^2 + (pi R T mu/(W K)) (L^2 - x^2).
    tensor_linear_2d  gate 4.3: 2-D square 0.05 x 0.05 m, N x N cells (default 20; the plan's 40
                      takes ~10 min, see KINDS), one block, four 303 faces (hconv 0, qrad 0,
                      Tref 300) each with its own p_w on the field p^2 = p0^2 + a x (p0 = 101325,
                      a = 2e11 Pa^2/m), K = R diag(K1, K2) R^T with R the 45-degree rotation in
                      x-y (--tensor45 1e-13 1e-15, the default), or the diagonal K = diag(K1, K2)
                      with --identity (the red control: u_y = 0 exactly), pi-source = 0,
                      quasi-steady, solver pseudo with pseudo-vnn 0.25, tol 1e-13, T(0) = 300 K.

    --tables <dir>    the isothermal-solid + constant-gas set of make_tacot_tables.py
                      (--solid-isothermal --constant-gas 0.022 4e-5 0 [--cp-override 1e12]); the
                      permeability lines of the copy are rewritten from --K / --tensor45.
    --ny N            the 1-D kinds: N identical lines along j (default 2), each 0.1 m tall; four lines
                      split into two blocks side by side make a correction group whose blocks share
                      their line positions (plan 11 stage 3b: test/pyrolysis/darcy/mass_balance_1d_2x).
    --depth N         MASTERDIR depth below the repo root (default from the case path under test/).

Python 3.6, no third-party modules. Imports make_case.py (same directory) for the shared writers.
"""

import argparse
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_case as mc  # noqa: E402

L_SLAB = 0.05            # m, every kind
DZ = 0.1                 # m, one cell along k (ndir = 2 like every existing case)
DY_1D = 0.1              # m, two cells along j in the 1-D kinds (make_case.py's geometry)
NY_1D = 2
P0 = 101325.0
R_UNIV = 8.314471469     # FUSS_Parameters_m
A_LIN = 2.0e11           # Pa^2/m, the p^2 slope of tensor_linear_2d

KINDS = {
    "sealed_nosource": dict(ncells=100, T0=1000.0, quasi=False, K=1.0e-11, pi=0.0, tol=1.0e-8, iters=200, solver="tridiag"),
    "mass_balance_1d": dict(ncells=100, T0=1000.0, quasi=False, K=1.0e-11, pi=0.1, tol=1.0e-8, iters=None, solver="tridiag"),
    "sealed_profile": dict(ncells=100, T0=1000.0, quasi=True, K=1.0e-13, pi=0.1, tol=1.0e-13, iters=400, solver="tridiag"),
    # tensor_linear_2d: 20 x 20 cells and pseudo-vnn 0.25 (measured 2026-10-08: the damped Jacobi sweep
    # at the plan's 0.5 diverges on the 45-degree 100:1 tensor at 40 x 40 -- the lagged cross terms
    # make 0.5 marginal -- while 0.25 converges in 4.6e4 sweeps per step at 40 x 40 (40 s per step,
    # ~10 min for the Picard steps) and in ~1.1e4 at 20 x 20; the exactness test does not depend on N)
    "tensor_linear_2d": dict(ncells=20, T0=300.0, quasi=True, K=None, pi=0.0, tol=1.0e-13, iters=200, solver="pseudo",
                             pseudo_vnn=0.25),
}

VAL_FMT = "    %13.6E,"
VAL3_FMT = "    %13.6E, %13.6E, %13.6E,"
# p_w with 16 significant digits: the 4.3 field is exact only for exact face values (a 7-digit p_w is off by
# up to 0.05 Pa, which showed as a 4e-7 relative error and |u_y|/|u_x| off by 3e-7, measured 2026-10-08)
VAL4_FMT = "    %13.6E, %13.6E, %13.6E, %.15E,"


def die(msg):
    sys.stderr.write("[ERROR] make_darcy_case: %s\n" % msg)
    sys.exit(1)


def fmt(v):
    return " %.15E" % v


# ---------------------------------------------------------------- mesh / ic (generate.py formats)
def write_ic(path, xs, ys, zs, zone_name, with_vars, T0):
    ni, nj, nk = len(xs), len(ys), len(zs)
    lines = []
    if with_vars:
        lines.append(' VARIABLES ="x" "y" "z" T matID')
        varloc = "VARLOCATION=([1-3]=NODAL,[4-5]=CELLCENTERED)"
    else:
        lines.append(' VARIABLES ="x" "y" "z"')
        varloc = "VARLOCATION=([1-3]=NODAL)"
    lines.append(" ZONE  T = %s, I=%d, J=%d, K=%d, DATAPACKING=BLOCK, %s,"
                 " SOLUTIONTIME=0.000000000000000E+000, STRANDID = 0" % (zone_name, ni, nj, nk, varloc))
    for axis in range(3):
        for k in range(nk):
            for j in range(nj):
                for i in range(ni):
                    lines.append(fmt((xs[i], ys[j], zs[k])[axis]))
    ncell = (ni - 1) * (nj - 1) * (nk - 1)
    if with_vars:
        for _ in range(ncell):
            lines.append(fmt(T0))
        for _ in range(ncell):
            lines.append(fmt(1.0))
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


# ---------------------------------------------------------------- bc.txt
def write_bc(path, nx, ny, face_spec):
    """face_spec[f] = function(i, j, k) -> (type, value-line or None); records in face order 1..6."""
    out = []

    def rec(i, j, k, f):
        typ, line = face_spec[f](i, j, k)
        out.append("%8d%8d%8d%8d%8d%8d" % (1, i, j, k, f, typ))
        if line is not None:
            out.append(line)

    for j in range(1, ny + 1):
        rec(1, j, 1, 1)
    for j in range(1, ny + 1):
        rec(nx, j, 1, 2)
    for i in range(1, nx + 1):
        rec(i, 1, 1, 3)
    for i in range(1, nx + 1):
        rec(i, ny, 1, 4)
    for j in range(1, ny + 1):
        for i in range(1, nx + 1):
            rec(i, j, 1, 5)
    for j in range(1, ny + 1):
        for i in range(1, nx + 1):
            rec(i, j, 1, 6)
    with open(path, "w") as fh:
        fh.write("\n".join(out) + "\n")
    return len([l for l in out if len(l.split()) == 6 and l.split()[0] == "1"])


# ---------------------------------------------------------------- permeability
def write_permeability(path, Kv, Kc):
    """Rewrite the K_v / K_c lines of a pyrolysis.dat with the 9-component tensors (row-major)."""
    out, seen = [], set()
    lines = {"K_v": "K_v %s\n" % " ".join("%.15g" % v for v in Kv), "K_c": "K_c %s\n" % " ".join("%.15g" % v for v in Kc)}
    with open(path) as fh:
        for s in fh:
            key = s.split("//", 1)[0].split()
            k = key[0] if key else ""
            if k in lines:
                out.append(lines[k])
                seen.add(k)
            else:
                out.append(s)
    for k in ("K_v", "K_c"):
        if k not in seen:
            out.append(lines[k])
    with open(path, "w") as fh:
        fh.write("".join(out))


def rotated_tensor(k1, k2, identity):
    """K = R diag(k1, k2, k2) R^T with R the 45-degree rotation about z (or diag(k1, k2, k2) with identity)."""
    if identity:
        return [k1, 0.0, 0.0, 0.0, k2, 0.0, 0.0, 0.0, k2]
    s = 0.5 * (k1 + k2)
    d = 0.5 * (k1 - k2)
    return [s, d, 0.0, d, s, 0.0, 0.0, 0.0, k2]


# ---------------------------------------------------------------- input.ini
def write_ini(path, kind, cfg, subst, probes_T=True):
    head = {
        "sealed_nosource": [
            "Sealed, source-free slab under the type-2 gas closure  (plan 07, gate 4.1(a))",
            "",
            "%(ncells)d cells x %(dx).6g m, both i-faces impermeable (BC 301, q 0), no gas source, T(0) at",
            "%(T0)g K, up to %(iters)d material steps (the run ends at its first step: nothing moves and the",
            "convergence residual is exactly zero). The pressure must stay at p-init to round-off and every",
            "face flux at zero (verify.py reads OUTPUT/field.tec and OUTPUT/darcy.txt).",
        ],
        "mass_balance_1d": [
            "Steady gas mass balance on a 1-D line  (plan 07, gate 4.1(b%(variant)s))",
            "",
            "%(ncells)d cells x %(dx).6g m, face 1 sealed (BC 301, q 0), face 2 open (BC 303, hconv 0,",
            "qrad 0, Tref %(T0)g, p_w %(pw)g Pa), constant gas production pi-source %(pi)g kg/m3/s,",
            "K %(K)g m2, constant synthetic gas, isothermal solid tables, T(0) at %(T0)g K, 2 s",
            "time-accurate%(qs)s. At steady state the mass flux through the face at x is pi x",
            "(verify.py rebuilds the face fluxes from OUTPUT/field.tec).",
        ],
        "sealed_profile": [
            "Analytic steady pressure profile  (plan 07, gate 4.2, grid N%(ncells)03d)",
            "",
            "%(ncells)d cells x %(dx).6g m, face 1 sealed, face 2 open at p_w %(pw)g Pa, pi-source %(pi)g,",
            "K %(K)g m2, constant gas (M 0.022, mu 4e-5, h 0), solid cp 1e12 (T frozen at %(T0)g K),",
            "quasi-steady: each material step is one Picard iteration on the lagged Gamma(p), run",
            "until p is stationary (%(iters)d steps or the residual stop). Exact:",
            "p^2 = p_w^2 + (pi R T mu/(W K)) (L^2 - x^2); verify.py fits the order over the grids.",
        ],
        "tensor_linear_2d": [
            "Orthotropic permeability on a 2-D square  (plan 07, gate 4.3%(ctrl)s)",
            "",
            "%(ncells)d x %(ncells)d cells, four open faces (BC 303, hconv 0, qrad 0, Tref %(T0)g) each with",
            "its own p_w on p^2 = p0^2 + a x (p0 %(pw)g Pa, a %(a)g Pa^2/m), K %(Kdesc)s, no gas source,",
            "quasi-steady, solver pseudo (open faces in two index directions), tol %(tol)g; the solid cp of",
            "1e12 freezes T. Exact: p^2 linear in x, u_y/u_x = K_yx/K_xx in every cell (verify.py).",
        ],
    }[kind]
    lines = ["; " + "=" * 77]
    for s in head:
        s = s % subst if "%(" in s else s
        lines.append(";  " + s if s else ";")
    lines.append(";")
    lines.append(";  Generated by test/pyrolysis/tools/make_darcy_case.py -- edit that, not this.")
    lines.append("; " + "=" * 77)
    if cfg["iters"] is not None:
        tparams = ["res-threshold = %s" % cfg["res"], "time-threshold = 1e30", "iter-threshold = %d" % cfg["iters"]]
    else:
        tparams = ["res-threshold = 1e-30", "time-threshold = %s" % cfg["tend"], "iter-threshold = 100000000"]
    darcy = ["[FUSS-Darcy]",
             "solver = %s" % cfg["solver"],
             "tol = %g" % cfg["tol"],
             "max-iter = %d" % cfg["max_iter"],
             "p-init = %.1f" % P0,
             "p-wall = %.1f" % P0,
             "quasi-steady = %s" % ("true" if cfg["quasi"] else "false"),
             "pi-source = %g" % cfg["pi"]]
    if cfg.get("pseudo_vnn") is not None:
        darcy.append("pseudo-vnn = %g" % cfg["pseudo_vnn"])
    lines += [
        "",
        "[FUSS-Parameters]",
        "newrun = true",
    ] + tparams + [
        "",
        "[FUSS-Numerics]",
        "time-scheme = RK3",
        "vnn = 0.5",
        "time-accurate = true",
        "integration-variables = cons",
        "irs = false",
        "",
        "[FUSS-Pyrolysis]",
        "gas-closure = type2",
        "",
    ] + darcy + [
        "",
        "[FUSS-IO]",
        "sol-dtime     = %s" % cfg["sol_dtime"],
        "sol-overwrite = true",
        "sol-format    = tecplot ascii",
        "ini-diter     = 100000000",
    ]
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case_dir")
    ap.add_argument("--kind", required=True, choices=sorted(KINDS))
    ap.add_argument("--tables", required=True, help="directory holding pyrolysis.dat virgin.dat char.dat gas.dat")
    ap.add_argument("--ncells", type=int, default=None)
    ap.add_argument("--quasi-steady", action="store_true", help="mass_balance_1d: variant (b1)")
    ap.add_argument("--K", type=float, default=None, help="isotropic permeability [m2] (1-D kinds)")
    ap.add_argument("--tensor45", type=float, nargs=2, default=[1.0e-13, 1.0e-15], metavar=("K1", "K2"),
                    help="tensor_linear_2d: principal permeabilities of the 45-degree tensor")
    ap.add_argument("--identity", action="store_true", help="tensor_linear_2d: K = diag(K1, K2) (red control)")
    ap.add_argument("--ny", type=int, default=None, help="1-D kinds: identical lines along j (default 2)")
    ap.add_argument("--pseudo-vnn", type=float, default=None, help="damping of the Jacobi sweep (tensor_linear_2d default 0.25)")
    ap.add_argument("--fuss-sh", default=mc.DONOR_FUSS_SH)
    ap.add_argument("--depth", type=int, default=None)
    args = ap.parse_args()

    kind = args.kind
    cfg = dict(KINDS[kind])
    if args.ncells is not None:
        if args.ncells < 2:
            die("--ncells must be at least 2")
        cfg["ncells"] = args.ncells
    if args.K is not None:
        if args.K <= 0:
            die("--K must be > 0")
        cfg["K"] = args.K
    if kind == "mass_balance_1d":
        cfg["quasi"] = args.quasi_steady
    if args.ny is not None:
        if kind == "tensor_linear_2d" or args.ny < 1:
            die("--ny is for the 1-D kinds, at least 1")
        cfg["ny"] = args.ny
    if args.pseudo_vnn is not None:
        if not 0.0 < args.pseudo_vnn <= 1.0:
            die("--pseudo-vnn must be in (0, 1]")
        cfg["pseudo_vnn"] = args.pseudo_vnn
    cfg["max_iter"] = 1000000 if kind == "tensor_linear_2d" else 100000
    # the convergence residual of a quasi-steady run is sum((p - p^n)/p-ref)^2: the Picard
    # increments bottom out at the line solve's round-off (~5e-13 relative at 160 cells,
    # measured); the stop compares the ROOT of that sum with res-threshold, so 1e-10 (a last
    # increment of at most ~1e-10 relative per cell) is reachable and the cap is not needed
    cfg["res"] = "1e-10"
    cfg["tend"] = "2.0"
    cfg["sol_dtime"] = "0.1" if kind == "mass_balance_1d" else "1e30"
    for f in mc.TABLE_FILES:
        if not os.path.isfile(os.path.join(args.tables, f)):
            die("%s not found in --tables %s" % (f, args.tables))

    case = os.path.abspath(args.case_dir)
    os.makedirs(os.path.join(case, "INPUT"), exist_ok=True)
    os.makedirs(os.path.join(case, "MESH"), exist_ok=True)
    for f in mc.TABLE_FILES:
        mc.shutil.copy(os.path.join(args.tables, f), os.path.join(case, "INPUT", f))
    pyro = os.path.join(case, "INPUT", "pyrolysis.dat")
    mc.strip_kinetics(pyro)
    if kind == "tensor_linear_2d":
        Kt = rotated_tensor(args.tensor45[0], args.tensor45[1], args.identity)
        write_permeability(pyro, Kt, Kt)
        Kdesc = ("diag(%g, %g) (red control)" % tuple(args.tensor45) if args.identity
                 else "R diag(%g, %g) R^T, R the 45-degree rotation" % tuple(args.tensor45))
    else:
        K = cfg["K"]
        write_permeability(pyro, [K, 0, 0, 0, K, 0, 0, 0, K], [K, 0, 0, 0, K, 0, 0, 0, K])
        Kdesc = "%g m2 isotropic" % K
    name, rho_v, nreac = mc.read_pyrolysis_dat(pyro)
    virgin = mc.read_table_zone1(os.path.join(case, "INPUT", "virgin.dat"), 7)

    # mesh
    nx = cfg["ncells"]
    xs = [L_SLAB * i / nx for i in range(nx + 1)]
    xs[-1] = L_SLAB
    if kind == "tensor_linear_2d":
        ny = nx
        ys = [L_SLAB * j / ny for j in range(ny + 1)]
        ys[-1] = L_SLAB
    else:
        ny = cfg.get("ny", NY_1D)
        ys = [DY_1D * j for j in range(ny + 1)]
    zs = [0.0, DZ]
    T0 = cfg["T0"]
    write_ic(os.path.join(case, "INPUT", "ic.tec"), xs, ys, zs, "B1-SP", True, T0)
    write_ic(os.path.join(case, "MESH", "mesh.tec"), xs, ys, zs, "Block1", False, T0)

    # boundary conditions
    if kind == "tensor_linear_2d":
        xc = [0.5 * (xs[i - 1] + xs[i]) for i in range(1, nx + 1)]     # xc[i-1] = centre of cell i

        def pw_at(x):
            return math.sqrt(P0 * P0 + A_LIN * x)

        def open_303(x):
            return (303, VAL4_FMT % (0.0, 0.0, T0, pw_at(x)))

        spec = {1: lambda i, j, k: open_303(0.0),
                2: lambda i, j, k: open_303(L_SLAB),
                3: lambda i, j, k: open_303(xc[i - 1]),
                4: lambda i, j, k: open_303(xc[i - 1]),
                5: lambda i, j, k: (0, None),
                6: lambda i, j, k: (0, None)}
    else:
        sealed = (301, VAL_FMT % 0.0)
        if kind == "sealed_nosource":
            face2 = sealed
        else:
            face2 = (303, VAL4_FMT % (0.0, 0.0, T0, P0))
        spec = {1: lambda i, j, k: sealed,
                2: lambda i, j, k: face2,
                3: lambda i, j, k: (300, None),
                4: lambda i, j, k: (300, None),
                5: lambda i, j, k: (0, None),
                6: lambda i, j, k: (0, None)}
    nrec = write_bc(os.path.join(case, "INPUT", "bc.txt"), nx, ny, spec)

    with open(os.path.join(case, "INPUT", "phase.txt"), "w") as fh:
        fh.write("solid phase\n%s 1\n" % name)
    mc.write_properties_dat(os.path.join(case, "INPUT", "properties.dat"), name, virgin, rho_v)

    subst = {"ncells": nx, "dx": L_SLAB / nx, "T0": T0, "iters": cfg["iters"] or 0, "pw": P0, "pi": cfg["pi"],
             "K": cfg.get("K") or 0.0, "qs": ", quasi-steady" if cfg["quasi"] else "",
             "variant": "1" if cfg["quasi"] else "2", "a": A_LIN, "Kdesc": Kdesc, "tol": cfg["tol"],
             "ctrl": ", red control K diagonal" if args.identity else ""}
    write_ini(os.path.join(case, "input.ini"), kind, cfg, subst)
    depth, how = mc.case_depth(case, "adiabatic0d", 1000.0, args.depth)
    masterdir = mc.write_fuss_sh(os.path.join(case, "FUSS.sh"), args.fuss_sh, depth)

    meta = {"kind": kind, "L": L_SLAB, "ncells": nx, "ny": ny, "T0": T0, "p0": P0, "pi": cfg["pi"],
            "quasi_steady": cfg["quasi"], "tol": cfg["tol"], "solver": cfg["solver"]}
    if kind == "tensor_linear_2d":
        meta.update({"a": A_LIN, "K1": args.tensor45[0], "K2": args.tensor45[1], "identity": args.identity, "K": Kt})
    else:
        meta["K"] = cfg["K"]
    with open(os.path.join(case, "case.json"), "w") as fh:
        json.dump(meta, fh, indent=1, sort_keys=True)
        fh.write("\n")

    print("make_darcy_case.py: %s -> %s" % (kind, case))
    print("  material %s (rho_v = %g), tables from %s, kinetics stripped, K %s" % (name, rho_v, args.tables, Kdesc))
    print("  cells %d x %d x 1, dx = %.6g m, dy = %.6g m; bc.txt records = %d" % (nx, ny, L_SLAB / nx, ys[1] - ys[0], nrec))
    print("  T0 = %g K, pi-source = %g, quasi-steady = %s, solver = %s, tol = %g, steps = %s"
          % (T0, cfg["pi"], cfg["quasi"], cfg["solver"], cfg["tol"], cfg["iters"] if cfg["iters"] else "to 2 s"))
    print("  FUSS.sh: %s (depth %d from %s)" % (masterdir, depth, how))


if __name__ == "__main__":
    main()
