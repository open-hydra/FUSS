#!/usr/bin/env python3
"""make_abl_demo_base.py -- the two bases of the ablating-wall guard demos (plan 08 section 2.10).

    make_abl_demo_base.py <case-dir>              BASE_ABL      = test/pyrolysis/ablation/demo_base
    make_abl_demo_base.py <case-dir> --pyro       BASE_ABL_PYRO = test/pyrolysis/ablation/demo_base_pyro
    make_abl_demo_base.py <case-dir> --two-block  BASE_ABL_2BLK = test/pyrolysis/ablation/demo_base_2blk

BASE_ABL_2BLK (plan 08 stage B; the node-map rules need two blocks): block 1 is 4 x 4 cells over 0.02 x 0.02 m,
block 2 beside it in x is 4 x 8 cells over 0.02 x 0.04 m (taller), one cell of 0.01 m in z, the legacy material of
BASE_ABL in both. Block 1's face 2 meets block 2's face 1 over j = 1..4 (101); block 2's face 1 above that and every
other side face are adiabatic (301); block 1's face 1 is the receding 306 face (law nodeshift, sdot-override 1e-5),
whose nodes no other block shares. The demos move the 306 faces to the bottom (face 3), where the corner node at
x = 0.02 is shared: with only block 1's bottom ablating, the patch is open (rule (a)); with both, the shared corner
ends lines of 0.02 and 0.04 m (rule (b)).

Mesh: one block, L = 0.05 m along i in 20 equal cells, one cell in j (0.1 m) and k (0.1 m), so the line is 1-D
(ndir = 1). T(0) = 300 K, vnn 0.3, 5 steps. Faces: 1 = BC 301 (adiabatic), 2 = BC 306, 3/4 = 300, 5/6 = 0.

BASE_ABL: the legacy material `mat` (k 1 W/m/K, rho 1000 kg/m3, cp 1000 J/kg/K, so alpha 1e-6 m2/s; the first step is
dt = 0.5 vnn dx^2/alpha = 0.9375 s), no [FUSS-Pyrolysis]: the 306 face runs in kinematic mode. Its line 2 is
`0.0, 300.0, none, env_demo.dat`; INPUT/env_demo.dat = two chemistry-off rows (hconv 10 W/m2/K, T_edge 400 K), and
with no pyrolysis material the face radiates nothing, so it is the BC 303 kernel. The face recedes:
[FUSS-MeshMotion] law = nodeshift with [FUSS-Ablation] sdot-override = 1e-5 m/s (plan 08 stage B; 5 steps move it
47 um of its 50 mm).

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


def fmt(v):
    return mc.fmt(v)


def two_block(case, depth):
    """BASE_ABL_2BLK (see the module docstring)."""
    blocks = [(0.00, 0.02, 4, 0.02, 4), (0.02, 0.02, 4, 0.04, 8)]   # x0, lx, nx, ly, ny
    ic, mesh = [' VARIABLES ="x" "y" "z" T matID'], [' VARIABLES ="x" "y" "z"']
    for b, (x0, lx, nx, ly, ny) in enumerate(blocks, 1):
        xs = [x0 + lx * i / nx for i in range(nx + 1)]
        ys = [ly * j / ny for j in range(ny + 1)]
        zs = [0.0, 0.01]
        for out, name, with_vars in ((ic, "B%d-SP" % b, True), (mesh, "Block%d" % b, False)):
            varloc = "VARLOCATION=([1-3]=NODAL,[4-5]=CELLCENTERED)" if with_vars else "VARLOCATION=([1-3]=NODAL)"
            out.append(" ZONE  T = %s, I=%d, J=%d, K=%d, DATAPACKING=BLOCK, %s, SOLUTIONTIME=0.000000000000000E+000,"
                       " STRANDID = 0" % (name, nx + 1, ny + 1, 2, varloc))
            for axis in range(3):
                for k in range(2):
                    for j in range(ny + 1):
                        for i in range(nx + 1):
                            out.append(fmt((xs[i], ys[j], zs[k])[axis]))
            if with_vars:
                out += [fmt(T0)] * (nx * ny) + [fmt(1.0)] * (nx * ny)
    open(os.path.join(case, "INPUT", "ic.tec"), "w").write("\n".join(ic) + "\n")
    open(os.path.join(case, "MESH", "mesh.tec"), "w").write("\n".join(mesh) + "\n")
    out = []

    def rec(b, i, j, k, f, typ, second=None):
        out.append("%8d%8d%8d%8d%8d%8d" % (b, i, j, k, f, typ))
        if second is not None:
            out.append(second)
    q0 = mc.VAL_FMT % 0.0
    conn = lambda bs, i_s, js, fs: "%8d%8d%8d%8d%8d%8d%8d%8d%8d" % (bs, i_s, js, 1, fs, 1, 0, 0, 1)
    for j in range(1, 5):
        rec(1, 1, j, 1, 1, 306, "     0.0, 300.0, none, env_demo.dat")
    for j in range(1, 5):
        rec(1, 4, j, 1, 2, 101, conn(2, 1, j, 1))
    for i in range(1, 5):
        rec(1, i, 1, 1, 3, 301, q0)
    for i in range(1, 5):
        rec(1, i, 4, 1, 4, 301, q0)
    for f in (5, 6):
        for j in range(1, 5):
            for i in range(1, 5):
                rec(1, i, j, 1, f, 0)
    for j in range(1, 9):
        if j <= 4:
            rec(2, 1, j, 1, 1, 101, conn(1, 4, j, 2))
        else:
            rec(2, 1, j, 1, 1, 301, q0)
    for j in range(1, 9):
        rec(2, 4, j, 1, 2, 301, q0)
    for i in range(1, 5):
        rec(2, i, 1, 1, 3, 301, q0)
    for i in range(1, 5):
        rec(2, i, 8, 1, 4, 301, q0)
    for f in (5, 6):
        for j in range(1, 9):
            for i in range(1, 5):
                rec(2, i, j, 1, f, 0)
    open(os.path.join(case, "INPUT", "bc.txt"), "w").write("\n".join(out) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case_dir")
    ap.add_argument("--pyro", action="store_true")
    ap.add_argument("--two-block", action="store_true")
    ap.add_argument("--depth", type=int, default=4)
    args = ap.parse_args()
    if args.pyro and args.two_block:
        sys.exit("--pyro and --two-block are separate bases")
    case = os.path.abspath(args.case_dir)
    os.makedirs(os.path.join(case, "INPUT"), exist_ok=True)
    os.makedirs(os.path.join(case, "MESH"), exist_ok=True)
    xs = [L * i / N for i in range(N + 1)]
    xs[-1] = L
    if not args.two_block:
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
    if args.two_block:
        two_block(case, args.depth)
    else:
        mc.write_bc(os.path.join(case, "INPUT", "bc.txt"), N, (301, mc.VAL_FMT % 0.0), (306, line2), ny=1)
    head = ["; " + "=" * 77,
            ";  Ablating-wall guard demo base%s  (plan 08 section 2.10; test/run_guard_demos.sh)"
            % (" with pyrolysis" if args.pyro else (", two blocks" if args.two_block else "")),
            ";",
            (";  Two blocks: 4 x 4 cells over 0.02 m and, beside it, 4 x 8 over 0.02 x 0.04 m; block 1's"
             "\n;  outer face is the ablating wall BC 306." if args.two_block else
             ";  L 0.05 m in 20 equal cells, one cell across: 1-D. Face 2 is the ablating wall BC 306,"),
            (";  type-2 TACOT without kinetics, the toy B' table INPUT/bprime_demo.dat, chemistry on."
             if args.pyro else
             ";  kinematic mode (no pyrolysis): legacy material k 1, rho 1000, cp 1000; chemistry off;\n"
             ";  the face recedes at the TEST ONLY sdot-override of 1e-5 m/s (law nodeshift)."),
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
    if not args.pyro:
        body += ["[FUSS-MeshMotion]", "law = nodeshift", ""]
    body += ["[FUSS-Ablation]", "seb-tol = 1.0e-8"] + ([] if args.pyro else ["sdot-override = 1.0e-5"]) + \
            ["", "[FUSS-IO]", "sol-diter     = 100000000", "sol-overwrite = true",
             "sol-format    = tecplot ascii", "ini-diter     = 100000000"]
    with open(os.path.join(case, "input.ini"), "w") as fh:
        fh.write("\n".join(body) + "\n")
    mc.write_fuss_sh(os.path.join(case, "FUSS.sh"), mc.DONOR_FUSS_SH, args.depth)
    print("make_abl_demo_base.py -> %s (%s)" % (case, "pyrolysis, type 2" if args.pyro else "kinematic, legacy material"))


if __name__ == "__main__":
    main()
