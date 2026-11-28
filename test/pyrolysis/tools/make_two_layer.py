#!/usr/bin/env python3
"""
make_two_layer.py -- write a FUSS pyrolysis case of a two-material slab in ONE block, with an exact solution.

The face between the two materials lies inside the block, so the case exercises the interior-face conductivity of
the pyrolysis path at a material jump (Diffusive_Flux with the two cells' k_s) and the time step of the cells beside
it. The kinetics are absent (no reaction in pyrolysis.dat) and every property is constant, so the energy equation is
the linear heat equation of a composite slab and tools/two_layer_exact.py gives its exact solution.

    x = 0 (face 1, BC 301 q = 0)        x = A                          x = L (face 2, BC 302 at TW from t = 0)
    |  material 1 HIGH: k 8.0, rho_s 2800, cp 4000  |  material 2 LOW: k 0.4, rho_s 280, cp 1000  |
                     k ratio 20, alpha_HIGH = alpha_LOW / 2

alpha_HIGH < alpha_LOW: the LOW cells set the global step of the time-accurate run, so the LOW cell beside the jump
runs at the step its own k_s gives. With the mean of the two k_s at the jump (10.5 x k_LOW) that cell's update leaves
the stable range of RK3 above vnn ~ 0.08; with the two half cells in series (1.9 x k_LOW) it stays inside it at the
nominal limit. Uniform cells (equal widths on both sides), N along i (even: the jump is the face N/2), the case
geometry of make_case.py otherwise (2 cells of 0.1 m along j, 1 along k).

Variants (the defaults write the case above):
    --k-ratio R          HIGH takes k = 0.4 R, rho_s = 140 R, cp 4000 (alpha_HIGH = alpha_LOW / 2 at every R)
    --cells-high NH --cells-low NL   NH cells in HIGH, NL in LOW (h_HIGH / h_LOW = NL / NH), instead of --ncells
    --ny NY --dy DY      NY cells of DY m along j (the field stays uniform along j: faces 3, 4 are 300)

Written, in the layout of make_case.py (whose writers it uses):
    <case>/FUSS.sh, input.ini, case.json (the parameters tools/two_layer_exact.py and verify.py read)
    <case>/MESH/mesh.tec, INPUT/ic.tec (T0, matID 1 or 2 per cell), INPUT/bc.txt, INPUT/phase.txt
    <case>/INPUT/{pyrolysis,virgin,char,gas,properties}.dat   two materials; gas.dat = the --tables gas zones, once
                                                              per material

Usage:
    make_two_layer.py <case-dir> --ncells N [--vnn V] [--time-scheme RK3|RK2|euler] [--tend T] [--sol-dtime D]
                      [--type2] [--tables DIR] [--depth N] [--k-ratio R] [--cells-high NH --cells-low NL]
                      [--ny NY --dy DY]
    --type2   gas-closure = type2, the heated face open (open-bc-types = 302, p-wall): the pores' gas expands as the
              slab heats and leaves there; the gas storage and its enthalpy enter E, so there is no exact solution: the
              type-2 variant serves the block-independence twin (the cut of test/tools/split_blocks.py). One open
              direction: the pressure is solved line by line (a sealed block takes the damped Jacobi sweeps, ~3e4 a
              step here).
Python 3.6, no numpy.
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_case as mc  # noqa: E402

L = 0.01                     # m, slab thickness
A_FRAC = 0.5                 # the jump at x = A_FRAC L (face N/2)
T0 = 300.0
TW = 1000.0
MATS = [  # (name, k, rho_s, cp): material 1 at the back, 2 at the heated face
    ("HIGH", 8.0, 2800.0, 4000.0),
    ("LOW", 0.4, 280.0, 1000.0),
]
T_REF = 298.15               # h = cp (T - T_REF), the datum of the --isothermal tables
T_TOP = 1500                 # the tables stop here (the run stays in [T0, TW])
K_PERM = 1.0e-11             # m2, isotropic permeability of both materials (type 2 only)
DEFAULT_TABLES = os.path.join(mc.REPO, "test", "pyrolysis", "conduction", "graded100", "INPUT")


def write_mesh_ic(case, xs, matid, ny=mc.NY):
    """MESH/mesh.tec (nodes) and INPUT/ic.tec (T0 and the matID of every cell): make_case.write_ic's format."""
    mc.write_ic(os.path.join(case, "MESH", "mesh.tec"), xs, "Block1", False, T0, ny=ny)
    n = mc.write_ic(os.path.join(case, "INPUT", "ic.tec"), xs, "B1-SP", True, T0, ny=ny)
    lines = open(os.path.join(case, "INPUT", "ic.tec")).read().split("\n")
    nx = len(xs) - 1
    ncell = mc.NK * ny * nx
    # the last ncell value lines are the matIDs, cell order i fastest
    first = len(lines) - 1 - ncell
    for c in range(ncell):
        lines[first + c] = mc.fmt(float(matid[c % nx]))
    with open(os.path.join(case, "INPUT", "ic.tec"), "w") as fh:
        fh.write("\n".join(lines))
    return n


def solid_table(path):
    """virgin.dat / char.dat: one zone per material, the T rows of the --tables virgin.dat, constant cp and k,
    h = cp (T - T_REF), emissivity and absorptivity 0.8 (no wall reads them: 301 and 302)"""
    Ts = [int(round(T)) for T in sorted(mc.read_table_zone1(path, 7)) if T <= T_TOP]
    out = ['TITLE = "FUSS pyrolysis table"',
           'VARIABLES = "T" "cp" "h" "ki" "kj" "kk" "emissivity" "absorptivity"']
    for name, k, _rho, cp in MATS:
        out += ['ZONE T="%s"' % name, "I=%d, F=POINT" % len(Ts)]
        for T in Ts:
            out.append("%.1f %.9E %.9E %.9E %.9E %.9E %.9E %.9E" % (T, cp, cp * (T - T_REF), k, k, k, 0.8, 0.8))
    return "\n".join(out) + "\n"


def gas_table(path):
    """gas.dat: the zones of the --tables gas.dat (one per p level), once per material, titles renamed, rows up to
    T_TOP"""
    head, zones, cur = [], [], None
    for line in open(path):
        if line.startswith("ZONE"):
            cur = [line]
            zones.append(cur)
        elif cur is None:
            head.append(line)
        else:
            cur.append(line)
    out = list(head)
    for name, _k, _rho, _cp in MATS:
        for z in zones:
            title = z[0]
            p = title.split("p=", 1)[1]
            out.append('ZONE T="%s p=%s' % (name, p))
            rows = [r for r in z[2:] if r.strip() and float(r.split()[0]) <= T_TOP]
            out.append("I=%d, F=POINT\n" % len(rows))
            out.extend(rows)
    return "".join(out)


def pyrolysis_dat(type2):
    out = ["// make_two_layer.py: two materials, no reaction (nPyroReac 0), rho_s = sum eps_I rho_I",
           "// (two phases of eps_I 0.1 each, eps_g 0.8, as the --isothermal TACOT tables)"]
    for name, _k, rho, _cp in MATS:
        out += ["material %s" % name, "nSolidPhases 2",
                "rhoI[1] %r" % (rho / 0.2), "epsI[1] 0.1", "nPyroReac[1] 0",
                "rhoI[2] %r" % (rho / 0.2), "epsI[2] 0.1", "nPyroReac[2] 0",
                "eps_g_v 0.8", "eps_g_c 0.8"]
        if type2:   # required by type 2
            out += ["K_v %g 0 0 0 %g 0 0 0 %g" % (K_PERM, K_PERM, K_PERM), "K_c %g 0 0 0 %g 0 0 0 %g" % (K_PERM, K_PERM, K_PERM)]
    return "\n".join(out) + "\n"


def properties_dat():
    lines = ['TITLE = "Mass Thermodynamic Properties"',
             'VARIABLES = "Temperature", "Cp", "Density", "Conductivity", "Energy"']
    for name, k, rho, cp in MATS:
        lines += ['ZONE T="%s"' % name, "I=4000, F=POINT"]
        for T in range(1, 4001):
            lines.append("%.1f %f %f %f %f" % (T, cp, rho, k, rho * cp * T))
    return "\n".join(lines) + "\n"


def write_ini(path, n, vnn, scheme, tend, sol_dtime, type2, cells_text=None, kratio_text="20"):
    pyro = (["[FUSS-Pyrolysis]", "gas-closure = type2", "", "[FUSS-Darcy]", "solver = tridiag",
             "p-init = 101325.0", "p-wall = 101325.0", "open-bc-types = 302"] if type2 else
            ["[FUSS-Pyrolysis]", "gas-closure = type1", "outflow-direction = +i", "p-gas = 101325.0"])
    body = [
        "; " + "=" * 77,
        ";  Two-material slab in one block (k ratio %s), no kinetics, constant properties: the composite-slab"
        % kratio_text,
        (";  heat equation, exact solution tools/two_layer_exact.py. %d equal cells, the jump at face %d." % (n, n // 2))
        if cells_text is None else (";  heat equation, exact solution tools/two_layer_exact.py. %s." % cells_text),
        ";  Face 1 adiabatic (301), face 2 at %g K from t = 0 (302), T(0) = %g K." % (TW, T0),
        ";  Gas closure %s%s." % ("type2" if type2 else "type1",
                                  ", the heated face open to the gas (302)" if type2 else ""),
        ";",
        ";  Generated by test/pyrolysis/tools/make_two_layer.py -- edit that, not this.",
        "; " + "=" * 77,
        "",
        "[FUSS-Parameters]",
        "newrun = true",
        "res-threshold = 1e-30",
        "time-threshold = %s" % mc.fmt_time(tend),
        "iter-threshold = 100000000",
        "",
        "[FUSS-Numerics]",
        "time-scheme = %s" % scheme,
        "vnn = %r" % vnn,
        "time-accurate = true",
        "integration-variables = cons",
        "",
    ] + pyro + [
        "",
        "[FUSS-IO]",
        "sol-dtime     = %s" % mc.fmt_time(sol_dtime),
        "sol-overwrite = false",
        "sol-format    = tecplot ascii",
        "ini-diter     = 100000000",
    ]
    with open(path, "w") as fh:
        fh.write("\n".join(body) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case_dir")
    ap.add_argument("--ncells", type=int, default=None)
    ap.add_argument("--vnn", type=float, default=0.05)
    ap.add_argument("--time-scheme", default="RK3", choices=("RK3", "RK2", "euler"))
    ap.add_argument("--tend", type=float, default=8.0)
    ap.add_argument("--sol-dtime", type=float, default=2.0)
    ap.add_argument("--type2", action="store_true")
    ap.add_argument("--tables", default=DEFAULT_TABLES, help="the --isothermal tables (virgin.dat, gas.dat)")
    ap.add_argument("--fuss-sh", default=mc.DONOR_FUSS_SH)
    ap.add_argument("--depth", type=int, default=None)
    ap.add_argument("--k-ratio", type=float, default=None, help="k_HIGH / k_LOW (default 20)")
    ap.add_argument("--cells-high", type=int, default=None)
    ap.add_argument("--cells-low", type=int, default=None)
    ap.add_argument("--ny", type=int, default=mc.NY)
    ap.add_argument("--dy", type=float, default=mc.DY)
    args = ap.parse_args()
    kratio_text = "20"
    if args.k_ratio is not None and args.k_ratio != 20.0:
        R = args.k_ratio
        MATS[0] = ("HIGH", 0.4 * R, 140.0 * R, 4000.0)
        kratio_text = "%g" % R
    mc.DY = args.dy
    cells_text = None
    if args.cells_high is not None or args.cells_low is not None:
        if args.cells_high is None or args.cells_low is None or args.ncells is not None:
            mc.die("--cells-high and --cells-low go together, without --ncells")
        nh, nl = args.cells_high, args.cells_low
        n = nh + nl
        A = A_FRAC * L
        xs = [A * i / nh for i in range(nh)] + [A + (L - A) * i / nl for i in range(nl)] + [L]
        ia = nh
        cells_text = "%d + %d cells (h_HIGH / h_LOW = %g), the jump at face %d" % (nh, nl, nl / nh, nh)
    else:
        if args.ncells is None:
            mc.die("--ncells N, or --cells-high NH --cells-low NL")
        n = args.ncells
        if n < 2 or n % 2:
            mc.die("--ncells must be even (the jump is the face N/2)")
        xs = [L * i / n for i in range(n + 1)]
        xs[-1] = L
        ia = int(round(A_FRAC * n))
    if args.ny != mc.NY or args.dy != 0.1:
        cells_text = (cells_text if cells_text is not None else
                      "%d equal cells, the jump at face %d" % (n, n // 2)) + "; %d cells of %g m along j" % (args.ny, args.dy)
    case = args.case_dir
    for d in ("INPUT", "MESH"):
        os.makedirs(os.path.join(case, d), exist_ok=True)

    matid = [1 if i < ia else 2 for i in range(n)]
    write_mesh_ic(case, xs, matid, ny=args.ny)
    nrec = mc.write_bc(os.path.join(case, "INPUT", "bc.txt"), n, (301, mc.VAL_FMT % 0.0), (302, mc.VAL_FMT % TW),
                       ny=args.ny)
    with open(os.path.join(case, "INPUT", "phase.txt"), "w") as fh:
        fh.write("solid phase\n" + "".join("%s 1\n" % m[0] for m in MATS))
    for name, text in (("virgin.dat", solid_table(os.path.join(args.tables, "virgin.dat"))),
                       ("char.dat", solid_table(os.path.join(args.tables, "virgin.dat"))),
                       ("gas.dat", gas_table(os.path.join(args.tables, "gas.dat"))),
                       ("pyrolysis.dat", pyrolysis_dat(args.type2)),
                       ("properties.dat", properties_dat())):
        with open(os.path.join(case, "INPUT", name), "w") as fh:
            fh.write(text)
    write_ini(os.path.join(case, "input.ini"), n, args.vnn, args.time_scheme, args.tend, args.sol_dtime, args.type2,
              cells_text, kratio_text)
    cfg = {"L": L, "A": xs[ia], "T0": T0, "TW": TW, "ncells": n,
           "materials": [{"name": m[0], "k": m[1], "rho": m[2], "cp": m[3]} for m in MATS],
           "vnn": args.vnn, "time_scheme": args.time_scheme, "type2": args.type2}
    if cells_text is not None:
        cfg.update({"cells": [ia, n - ia], "ny": args.ny, "dy": args.dy})
    with open(os.path.join(case, "case.json"), "w") as fh:
        json.dump(cfg, fh, indent=1, sort_keys=True)
        fh.write("\n")
    depth, how = mc.case_depth(case, "conduction", T0, args.depth)
    master = mc.write_fuss_sh(os.path.join(case, "FUSS.sh"), args.fuss_sh, depth)
    aL = MATS[1][1] / (MATS[1][2] * MATS[1][3])
    hmin = min(xs[i + 1] - xs[i] for i in range(n))
    dt = args.vnn / (aL * (1.0 / hmin ** 2 + 1.0 / mc.DY ** 2))
    print("make_two_layer.py: %d cells (jump at face %d, x = %g m), %s %s vnn %g, dt ~ %.3e s -> ~%d steps; "
          "bc.txt records %d; FUSS.sh %s (%s)" % (n, ia, xs[ia], "type2" if args.type2 else "type1",
                                                  args.time_scheme, args.vnn, dt, int(args.tend / dt), nrec, master, how))


if __name__ == "__main__":
    main()
