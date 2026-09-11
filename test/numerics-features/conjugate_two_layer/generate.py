#!/usr/bin/env python3
"""
Generate the conjugate two-layer case (plan 10, D4).

WHAT THIS CASE IS FOR
---------------------
It answers one question: does FUSS's block-to-block interface flux actually use
the NEIGHBOUR's geometry?

`BC_Connection` sets the interface temperature as a conductivity-and-metric
weighted blend of the two sides,

    T_int = (k1 T1 D1 + k2 T2 D2) / (k1 D1 + k2 D2)

where D1 and D2 come from each side's metric tensor dotted with the face
normal. D2 is read from bc%Mg -- the only place in the solver that field is
used. For a two-layer slab the same expression is exactly the flux-continuity
interface temperature between two cell centres at dx1/2 and dx2/2, so IF the
neighbour's metric is right, the scheme reproduces the analytic profile to
round-off, and if it is wrong the profile kinks at the interface.

No existing test could tell: every multi-block case in the suite has the SAME
cell spacing on both sides of every interface, which is the one configuration
where getting the neighbour's metric wrong is invisible. Real TPS stacks do not
look like that.

THE DESIGN, AND WHY EACH CHOICE MATTERS
---------------------------------------
  * dx2 / dx1 = 4. Unequal spacing across the interface is the whole point.
  * k1 / k2 = 10. Different materials, so the weights differ by conductivity
    as well as by spacing: k1/dx1 = 10000 against k2/dx2 = 250, a 40:1 split.
  * Both together make the case sensitive. Substituting the LOCAL metric for
    the neighbour's moves the weights from 40:1 to 10:1 and the interface
    temperature by about 4.3 K -- four orders of magnitude above the round-off
    the correct scheme achieves.
  * The exact solution is piecewise LINEAR, so it is representable on the grid
    exactly. There is no discretisation error to hide behind: the expected
    result is machine precision, not a tolerance.

Analytic, with R = L/k the thermal resistance of each layer:

    q     = (T_hot - T_cold) / (R1 + R2)
    T_int = T_hot - q R1
    T(x)  linear within each layer

Run this from the directory it lives in; it rewrites INPUT/ and MESH/.
"""

import os
import shutil
import stat

HERE = os.path.dirname(os.path.abspath(__file__))
DONOR = os.path.abspath(os.path.join(HERE, "..", "..", "moving-mesh", "freestream"))

# ------------------------------------------------------------------ geometry
L1, NX1 = 0.10, 20          # dx1 = 5.0e-3
L2, NX2 = 0.10, 5           # dx2 = 2.0e-2   -> 4x coarser across the interface
HY, NY = 0.02, 4
WZ = 0.01

# ----------------------------------------------------------------- materials
CP, RHO = 1000.0, 1000.0
K1, K2 = 50.0, 5.0
TMAX, NTAB = 2000, 2000

# ---------------------------------------------------------------------- state
T_HOT, T_COLD = 1000.0, 300.0
T_INIT = 300.0

DX1, DX2 = L1 / NX1, L2 / NX2
R1, R2 = L1 / K1, L2 / K2
Q = (T_HOT - T_COLD) / (R1 + R2)
T_INT = T_HOT - Q * R1


def t_exact(x):
    """Analytic steady profile. x measured from the hot face."""
    if x <= L1:
        return T_HOT - Q * x / K1
    return T_INT - Q * (x - L1) / K2


def fmt(v):
    return " %.15E" % v


def block_nodes(x0, lx, nx):
    xs = [x0 + lx * i / nx for i in range(nx + 1)]
    ys = [HY * j / NY for j in range(NY + 1)]
    zs = [0.0, WZ]
    return xs, ys, zs


def zone_text(name, xs, ys, zs, matid, with_vars):
    ni, nj, nk = len(xs), len(ys), len(zs)
    varloc = ("VARLOCATION=([1-3]=NODAL,[4-5]=CELLCENTERED)" if with_vars
              else "VARLOCATION=([1-3]=NODAL)")
    out = [" ZONE  T = %s, I=%d, J=%d, K=%d, DATAPACKING=BLOCK, %s,"
           " SOLUTIONTIME=0.000000000000000E+000, STRANDID = 0"
           % (name, ni, nj, nk, varloc)]

    for axis in range(3):
        for k in range(nk):
            for j in range(nj):
                for i in range(ni):
                    out.append(fmt((xs[i], ys[j], zs[k])[axis]))

    if with_vars:
        # Uniform initial condition -- deliberately NOT the analytic profile,
        # so the run has to converge onto it rather than start there.
        for _k in range(nk - 1):
            for _j in range(nj - 1):
                for _i in range(ni - 1):
                    out.append(fmt(T_INIT))
        for _ in range((nk - 1) * (nj - 1) * (ni - 1)):
            out.append(fmt(float(matid)))
    return out


def write_tec(path, with_vars):
    head = [' VARIABLES ="x" "y" "z" T matID' if with_vars
            else ' VARIABLES ="x" "y" "z"']
    b1 = block_nodes(0.0, L1, NX1)
    b2 = block_nodes(L1, L2, NX2)
    tag = "-SP" if with_vars else ""
    lines = head
    lines += zone_text("B1" + tag if with_vars else "Block1", *b1,
                       matid=1, with_vars=with_vars)
    lines += zone_text("B2" + tag if with_vars else "Block2", *b2,
                       matid=2, with_vars=with_vars)
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def write_properties(path):
    lines = ['TITLE = "Mass Thermodynamic Properties"',
             'VARIABLES = "Temperature", "Cp", "Density", "Conductivity", "Energy"']
    for name, k in (("layerA", K1), ("layerB", K2)):
        lines.append('ZONE T="%s"' % name)
        lines.append("I=%d, F=POINT" % NTAB)
        for i in range(1, NTAB + 1):
            t = float(i)
            lines.append("%.1f %f %f %f %f" % (t, CP, RHO, k, RHO * CP * t))
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def write_bc(path):
    """Records ordered block-major, face-minor -- the order Mod_BC_Fluxes groups by.

    Block 1 face 2 and block 2 face 1 are the conjugate interface, written as
    type 101 with an identity orientation map (d11 d12 d21 d22 = 1 0 0 1).
    """
    out = []

    def rec(b, i, j, k, f, typ, second=None):
        out.append("%8d%8d%8d%8d%8d%8d" % (b, i, j, k, f, typ))
        if second is not None:
            out.append(second)

    def val(v):
        return "    %13.6E," % v

    def conn(bs, i_s, js, ks, fs):
        return "%8d%8d%8d%8d%8d%8d%8d%8d%8d" % (bs, i_s, js, ks, fs, 1, 0, 0, 1)

    # ---- block 1
    for j in range(1, NY + 1):                       # face 1: hot wall
        rec(1, 1, j, 1, 1, 302, val(T_HOT))
    for j in range(1, NY + 1):                       # face 2: -> block 2 face 1
        rec(1, NX1, j, 1, 2, 101, conn(2, 1, j, 1, 1))
    for i in range(1, NX1 + 1):                      # face 3
        rec(1, i, 1, 1, 3, 301, val(0.0))
    for i in range(1, NX1 + 1):                      # face 4
        rec(1, i, NY, 1, 4, 301, val(0.0))
    for j in range(1, NY + 1):                       # face 5
        for i in range(1, NX1 + 1):
            rec(1, i, j, 1, 5, 0)
    for j in range(1, NY + 1):                       # face 6
        for i in range(1, NX1 + 1):
            rec(1, i, j, 1, 6, 0)

    # ---- block 2
    for j in range(1, NY + 1):                       # face 1: -> block 1 face 2
        rec(2, 1, j, 1, 1, 101, conn(1, NX1, j, 1, 2))
    for j in range(1, NY + 1):                       # face 2: cold wall
        rec(2, NX2, j, 1, 2, 302, val(T_COLD))
    for i in range(1, NX2 + 1):
        rec(2, i, 1, 1, 3, 301, val(0.0))
    for i in range(1, NX2 + 1):
        rec(2, i, NY, 1, 4, 301, val(0.0))
    for j in range(1, NY + 1):
        for i in range(1, NX2 + 1):
            rec(2, i, j, 1, 5, 0)
    for j in range(1, NY + 1):
        for i in range(1, NX2 + 1):
            rec(2, i, j, 1, 6, 0)

    with open(path, "w") as fh:
        fh.write("\n".join(out) + "\n")


INI = """\
; =============================================================================
;  Conjugate two-layer slab with UNEQUAL cell spacing across the interface
;  (plan 10, D4)
;
;  Generated by generate.py -- edit that, not this.
;
;  Two blocks, two materials, and dx2/dx1 = {ratio:g}. The steady solution is
;  piecewise linear and therefore exactly representable on this grid, so the
;  expected error is round-off and any interface mistreatment shows up as a
;  kink rather than as a tolerance argument.
;
;      q     = (Thot - Tcold) / (L1/k1 + L2/k2) = {q:.6f} W/m2
;      T_int = {tint:.6f} K
;
;  Every other multi-block case in the suite has matched spacing across its
;  interfaces, which is exactly the configuration in which using the wrong
;  side's metric is invisible. See ../../../test/run_regression.sh history and
;  plan 10 D4.
; =============================================================================

[FUSS-Parameters]
newrun = true
res-threshold = 1e-10

[FUSS-Numerics]
time-scheme = RK3
vnn = 2.0
time-accurate = false
integration-variables = cons
irs = true
irs-beta = 0.5

[FUSS-IO]
sol-diter     = 100000
sol-overwrite = true
sol-format    = tecplot ascii
ini-diter     = 100000

;

[GRIB-meshgen]
outpath = MESH/
method = fmsh
type = 2Dplane
width = {wz}

[GRIB-Block1]
pivot = 0.0 0.0 0.0
length = {l1}
height = {hy}
nx = {nx1}
ny = {ny}
strx = U 1.0
stry = U 1.0

[GRIB-Block2]
pivot = {l1} 0.0 0.0
length = {l2}
height = {hy}
nx = {nx2}
ny = {ny}
strx = U 1.0
stry = U 1.0

;

[GPB-Phase1]
type = solid
material = layerA layerB
cp  = {cp} {cp}
rho = {rho} {rho}
k   = {k1} {k2}
Tmax = {tmax}

;

[ICB-Block1]
type = homogeneous
material = layerA
T = {tinit}

[ICB-Block2]
type = homogeneous
material = layerB
T = {tinit}

;

[BCB-Block1]
face1 = Twall-hot
face2 = connection
face3 = adiabatic
face4 = adiabatic
face5 = null
face6 = null

[BCB-Block2]
face1 = connection
face2 = Twall-cold
face3 = adiabatic
face4 = adiabatic
face5 = null
face6 = null

[Twall-hot]
type = wall
T = {thot}

[Twall-cold]
type = wall
T = {tcold}

[adiabatic]
type = wall
q = 0.0
"""


def main():
    os.makedirs(os.path.join(HERE, "INPUT"), exist_ok=True)
    os.makedirs(os.path.join(HERE, "MESH"), exist_ok=True)

    write_tec(os.path.join(HERE, "INPUT", "ic.tec"), True)
    write_tec(os.path.join(HERE, "MESH", "mesh.tec"), False)
    write_bc(os.path.join(HERE, "INPUT", "bc.txt"))
    write_properties(os.path.join(HERE, "INPUT", "properties.dat"))

    with open(os.path.join(HERE, "INPUT", "phase.txt"), "w") as fh:
        fh.write("solid phase\nlayerA 1\nlayerB 1\n")

    with open(os.path.join(HERE, "input.ini"), "w") as fh:
        fh.write(INI.format(l1=L1, l2=L2, nx1=NX1, nx2=NX2, ny=NY, hy=HY,
                            wz=WZ, cp=CP, rho=RHO, k1=K1, k2=K2, tmax=TMAX,
                            thot=T_HOT, tcold=T_COLD, tinit=T_INIT,
                            q=Q, tint=T_INT, ratio=DX2 / DX1))

    with open(os.path.join(DONOR, "FUSS.sh")) as fh:
        sh = fh.read().replace("MASTERDIR=../../..", "MASTERDIR=../../..")
    dst = os.path.join(HERE, "FUSS.sh")
    with open(dst, "w") as fh:
        fh.write(sh)
    os.chmod(dst, os.stat(dst).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    print("dx1 = %.4e   dx2 = %.4e   ratio = %g" % (DX1, DX2, DX2 / DX1))
    print("k1/dx1 = %.1f   k2/dx2 = %.1f   weight ratio = %.1f:1"
          % (K1 / DX1, K2 / DX2, (K1 / DX1) / (K2 / DX2)))
    print("q = %.6f W/m2   T_int = %.6f K" % (Q, T_INT))


if __name__ == "__main__":
    main()
