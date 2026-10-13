#!/usr/bin/env python3
"""make_abl_demo_base.py -- the two bases of the ablating-wall guard demos (plan 08 section 2.10).

    make_abl_demo_base.py <case-dir>            BASE_ABL      = test/pyrolysis/ablation/demo_base
    make_abl_demo_base.py <case-dir> --pyro     BASE_ABL_PYRO = test/pyrolysis/ablation/demo_base_pyro

Mesh: one block, L = 0.05 m along i in 20 equal cells, one cell in j (0.1 m) and k (0.1 m), so the line is 1-D
(ndir = 1). T(0) = 300 K, vnn 0.3, 5 steps. Faces: 1 = BC 301 (adiabatic), 2 = BC 306, 3/4 = 300, 5/6 = 0.

BASE_ABL: the legacy material `mat` (k 1 W/m/K, rho 1000 kg/m3, cp 1000 J/kg/K, so alpha 1e-6 m2/s; the first step is
dt = 0.5 vnn dx^2/alpha = 0.9375 s), no [FUSS-Pyrolysis]: the 306 face runs in kinematic mode. Its line 2 is
`0.0, 300.0, none, env_demo.dat`; INPUT/env_demo.dat = two chemistry-off rows (hconv 10 W/m2/K, T_edge 400 K), and
with no pyrolysis material the face radiates nothing, so it is the BC 303 kernel.

BASE_ABL_PYRO: the same mesh under [FUSS-Pyrolysis] gas-closure = type2 with the tables of
test/pyrolysis/darcy/mass_balance_1d/INPUT/ (TACOT without kinetics, gas table 1.01325 .. 1.01325e6 Pa), the 306 line
`0.5, 300.0, bprime_demo.dat, env_demo.dat`, env rows with the chemistry on (rhoeUeCH0 0.1, h_r 2e6, p_w 101325) and
INPUT/bprime_demo.dat: `// kind: equilibrium`, two p-blocks (2e5, 5e4 Pa) of two B'g sub-blocks (10, 0.001) of two
T rows (300, 4000 K), B'c 0.1, h_w = -1e6 + 1000 T.

Each base runs to completion as committed (the demos must fail for their own reason, run_guard_demos.sh checks that).
Python 3.6, no third-party modules.
"""
import argparse
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_case as mc  # noqa: E402

N = 20
L = 0.05
T0 = 300.0
PYRO_TABLES = os.path.join(mc.REPO, "test", "pyrolysis", "darcy", "mass_balance_1d", "INPUT")
ENV_KIN = [(0.0, 101325.0, 0.0, 0.0, 0.0, 0.0, 10.0, 400.0), (1.0e6, 101325.0, 0.0, 0.0, 0.0, 0.0, 10.0, 400.0)]
ENV_PYRO = [(0.0, 101325.0, 0.1, 2.0e6, 0.0, 1.0, 0.0, 300.0), (1.0e6, 101325.0, 0.1, 2.0e6, 0.0, 1.0, 0.0, 300.0)]


def bprime_demo(path):
    with open(path, "w") as fh:
        fh.write("// kind: equilibrium\n// axes: p Bprime-g T\n")
        for p in (2.0e5, 5.0e4):
            for bg in (10.0, 0.001):
                for T in (300.0, 4000.0):
                    fh.write("%r %r %r %r %r\n" % (p, bg, 0.1, T, -1.0e6 + 1000.0 * T))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case_dir")
    ap.add_argument("--pyro", action="store_true")
    ap.add_argument("--depth", type=int, default=4)
    args = ap.parse_args()
    case = os.path.abspath(args.case_dir)
    os.makedirs(os.path.join(case, "INPUT"), exist_ok=True)
    os.makedirs(os.path.join(case, "MESH"), exist_ok=True)
    xs = [L * i / N for i in range(N + 1)]
    xs[-1] = L
    mc.write_ic(os.path.join(case, "INPUT", "ic.tec"), xs, "B1-SP", True, T0, ny=1)
    mc.write_ic(os.path.join(case, "MESH", "mesh.tec"), xs, "Block1", False, T0, ny=1)
    if args.pyro:
        for f in mc.TABLE_FILES:
            shutil.copy(os.path.join(PYRO_TABLES, f), os.path.join(case, "INPUT", f))
        name, rho_v, _ = mc.read_pyrolysis_dat(os.path.join(case, "INPUT", "pyrolysis.dat"))
        virgin = mc.read_table_zone1(os.path.join(case, "INPUT", "virgin.dat"), 7)
        mc.write_properties_dat(os.path.join(case, "INPUT", "properties.dat"), name, virgin, rho_v)
        bprime_demo(os.path.join(case, "INPUT", "bprime_demo.dat"))
        mc.write_env_file(os.path.join(case, "INPUT", "env_demo.dat"), ENV_PYRO)
        line2 = "     0.5, 300.0, bprime_demo.dat, env_demo.dat"
    else:
        name = "mat"
        with open(os.path.join(case, "INPUT", "properties.dat"), "w") as fh:
            fh.write('TITLE = "Mass Thermodynamic Properties"\n'
                     'VARIABLES = "Temperature", "Cp", "Density", "Conductivity", "Energy"\n'
                     'ZONE T="mat"\nI=4000, F=POINT\n')
            for T in range(1, 4001):
                fh.write("%.1f %f %f %f %f\n" % (T, 1000.0, 1000.0, 1.0, 1.0e6 * T))
        mc.write_env_file(os.path.join(case, "INPUT", "env_demo.dat"), ENV_KIN)
        line2 = "     0.0, 300.0, none, env_demo.dat"
    with open(os.path.join(case, "INPUT", "phase.txt"), "w") as fh:
        fh.write("solid phase\n%s 1\n" % name)
    mc.write_bc(os.path.join(case, "INPUT", "bc.txt"), N, (301, mc.VAL_FMT % 0.0), (306, line2), ny=1)
    head = ["; " + "=" * 77,
            ";  Ablating-wall guard demo base%s  (plan 08 section 2.10; test/run_guard_demos.sh)" % (" with pyrolysis" if args.pyro else ""),
            ";",
            ";  L 0.05 m in 20 equal cells, one cell across: 1-D. Face 2 is the ablating wall BC 306,",
            (";  type-2 TACOT without kinetics, the toy B' table INPUT/bprime_demo.dat, chemistry on."
             if args.pyro else
             ";  kinematic mode (no pyrolysis): legacy material k 1, rho 1000, cp 1000; chemistry off."),
            ";  5 steps. The demos edit a copy; this base itself runs to completion.",
            ";",
            ";  Generated by test/pyrolysis/tools/make_abl_demo_base.py -- edit that, not this.",
            "; " + "=" * 77]
    body = head + ["", "[FUSS-Parameters]", "newrun = true", "res-threshold = 1e-30", "time-threshold = 1.0e6",
                   "iter-threshold = 5", "", "[FUSS-Numerics]", "time-scheme = RK3", "vnn = 0.3", "time-accurate = true",
                   "integration-variables = cons", "irs = false", ""]
    if args.pyro:
        body += ["[FUSS-Pyrolysis]", "gas-closure = type2", "", "[FUSS-Darcy]", "solver = tridiag",
                 "p-init = 101325.0", "p-wall = 101325.0", ""]
    body += ["[FUSS-Ablation]", "seb-tol = 1.0e-8", "", "[FUSS-IO]", "sol-diter     = 100000000", "sol-overwrite = true",
             "sol-format    = tecplot ascii", "ini-diter     = 100000000"]
    with open(os.path.join(case, "input.ini"), "w") as fh:
        fh.write("\n".join(body) + "\n")
    mc.write_fuss_sh(os.path.join(case, "FUSS.sh"), mc.DONOR_FUSS_SH, args.depth)
    print("make_abl_demo_base.py -> %s (%s)" % (case, "pyrolysis, type 2" if args.pyro else "kinematic, legacy material"))


if __name__ == "__main__":
    main()
