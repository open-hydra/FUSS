#!/usr/bin/env python3
"""
Generate the grid sequence for the receding-surface conduction test
(plan 05, gate 5.3).

WHAT THE CASE IS
----------------
A slab of constant-property solid whose surface recedes into the material at a
constant rate. The whole mesh translates rigidly at that rate, so the surface
stays on the low-x boundary and the problem is steady in the mesh frame.

    lab frame   :  dT/dt = alpha d2T/dx2 ,  hot face at x = s(t) = v t
    mesh frame  :  xi = x - v t  ->  dT/dt = v dT/dxi + alpha d2T/dxi2

whose steady solution, with T(0) = Tw and T(L) = Tinf, is

    T(xi) = Tinf + (Tw - Tinf) * [ exp(-xi/d) - exp(-L/d) ] / [ 1 - exp(-L/d) ]

with the decay length d = alpha / v.  This is the oracle; it is exact for the
finite domain, so no far-field truncation is folded into the measured error.

WHY RIGID TRANSLATION AND NOT A COMPRESSING NODE-SHIFT
-----------------------------------------------------
Under rigid translation every cell keeps its volume, the swept volumes are
exactly +/- A*v*dt on the two faces normal to the motion, and the ALE remap
collapses to first-order upwind advection at Courant number C = v dt / dx. Its
truncation error is therefore known in closed form BEFORE the run -- which is
what "state the tolerance derived from the discretisation order, not chosen to
pass" requires. A node-shift law that compresses the mesh has a spatially
varying mesh velocity and a cell size that changes as the run proceeds, so the
same statement would only be approximate. Node-shifting is the plan-08 law and
is verified there, driven by a real recession rate.

WHY THE PARAMETERS ARE WHAT THEY ARE
------------------------------------
  * L = 10 d, so exp(-L/d) ~ 4.5e-5: the far boundary sits far outside the
    thermal layer and cannot contaminate the near-surface error.
  * The grids resolve the decay length with 4, 8, 16 and 32 cells, i.e. mesh
    Peclet numbers Pe = dx/d of 0.25 down to 0.03125. Coarse enough that the
    first-order upwind term dominates the second-order conduction term, which
    is what makes the measured ORDER meaningful.
  * The initial condition is the analytic profile itself, so the run only has
    to relax the difference between the continuum and discrete steady states.
  * end time = 200 s is 10 relaxation times of the slowest mode,
    tau = 1 / ( alpha (pi/L)^2 + v^2/(4 alpha) ) = 19.9 s.

Run this from the directory it lives in; it rewrites every N* case in place.
"""

import math
import os
import shutil
import stat

HERE = os.path.dirname(os.path.abspath(__file__))
DONOR = os.path.abspath(os.path.join(HERE, "..", "freestream"))

# ---------------------------------------------------------------- material
CP, RHO, K = 460.0, 7850.0, 52.0
ALPHA = K / (RHO * CP)                    # 1.4400443e-05 m^2/s

# ---------------------------------------------------------------- problem
DECAY = 0.01                              # d = alpha / v      [m]
VEL = ALPHA / DECAY                       # recession rate     [m/s]
LX = 10.0 * DECAY                         # slab thickness     [m]
HY = 0.02                                 # lateral extent     [m]
WZ = 0.01                                 # plane width        [m]
NY = 4                                    # lateral cells (solution is 1-D)
TW = 1273.15                              # receding surface   [K]
TINF = 273.15                             # far field          [K]
T_END = 200.0                             # end time           [s]
VNN = 0.5

GRIDS = [40, 80, 160, 320]


def t_exact(xi):
    """Steady mesh-frame profile on the FINITE domain [0, LX]."""
    e_l = math.exp(-LX / DECAY)
    return TINF + (TW - TINF) * (math.exp(-xi / DECAY) - e_l) / (1.0 - e_l)


def t_exact_cell_average(xa, xb):
    """Exact average of the profile over [xa, xb].

    Compared against the cell-centre point value this removes an O(dx^2)
    ambiguity about what a finite-volume unknown represents. It matters here
    because the quantity being measured is itself only first order, so a
    second-order bookkeeping difference would show up as curvature in the
    convergence plot.
    """
    e_l = math.exp(-LX / DECAY)
    integ = DECAY * (math.exp(-xa / DECAY) - math.exp(-xb / DECAY))
    return TINF + (TW - TINF) * (integ / (xb - xa) - e_l) / (1.0 - e_l)


def fmt(v):
    return " %.15E" % v


def write_ic(path, nx, zone_name, with_vars):
    ni, nj, nk = nx + 1, NY + 1, 2
    xs = [LX * i / nx for i in range(ni)]
    ys = [HY * j / NY for j in range(nj)]
    zs = [0.0, WZ]

    lines = []
    if with_vars:
        lines.append(' VARIABLES ="x" "y" "z" T matID')
        varloc = "VARLOCATION=([1-3]=NODAL,[4-5]=CELLCENTERED)"
    else:
        lines.append(' VARIABLES ="x" "y" "z"')
        varloc = "VARLOCATION=([1-3]=NODAL)"
    lines.append(
        " ZONE  T = %s, I=%d, J=%d, K=%d, DATAPACKING=BLOCK, %s,"
        " SOLUTIONTIME=0.000000000000000E+000, STRANDID = 0"
        % (zone_name, ni, nj, nk, varloc)
    )

    # Nodal coordinates, i fastest then j then k (Fortran order).
    for coord, axis in ((xs, 0), (ys, 1), (zs, 2)):
        for k in range(nk):
            for j in range(nj):
                for i in range(ni):
                    v = (xs[i], ys[j], zs[k])[axis]
                    lines.append(fmt(v))

    if with_vars:
        # Cell-centred T (the analytic profile) and matID.
        for _k in range(nk - 1):
            for _j in range(NY):
                for i in range(nx):
                    lines.append(fmt(t_exact_cell_average(xs[i], xs[i + 1])))
        for _ in range((nk - 1) * NY * nx):
            lines.append(fmt(1.0))

    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def write_bc(path, nx):
    """Boundary records: 'b i j k f type' then a type-dependent value line.

    Face 1 (low x)  = receding surface, fixed temperature   -> 302
    Face 2 (high x) = far field,        fixed temperature   -> 302
    Faces 3,4       = lateral,          adiabatic           -> 301
    Faces 5,6       = plane normal,     null                -> 0, no value line

    The count must match Allocate_BC's 2*nj*nk + 2*ni*nk + 2*ni*nj exactly, and
    the face order must be 1,2,3,4,5,6 -- Check_BC only counts records, so a
    wrong order would be read as valid and silently mislabel the boundaries.
    """
    nk = 1
    out = []

    def rec(i, j, k, f, typ, val=None):
        out.append("%8d%8d%8d%8d%8d%8d" % (1, i, j, k, f, typ))
        if val is not None:
            out.append("    %13.6E," % val)

    for j in range(1, NY + 1):
        for k in range(1, nk + 1):
            rec(1, j, k, 1, 302, TW)
    for j in range(1, NY + 1):
        for k in range(1, nk + 1):
            rec(nx, j, k, 2, 302, TINF)
    for i in range(1, nx + 1):
        for k in range(1, nk + 1):
            rec(i, 1, k, 3, 301, 0.0)
    for i in range(1, nx + 1):
        for k in range(1, nk + 1):
            rec(i, NY, k, 4, 301, 0.0)
    for j in range(1, NY + 1):
        for i in range(1, nx + 1):
            rec(i, j, 1, 5, 0)
    for j in range(1, NY + 1):
        for i in range(1, nx + 1):
            rec(i, j, 1, 6, 0)

    with open(path, "w") as fh:
        fh.write("\n".join(out) + "\n")


INI = """\
; =============================================================================
;  Receding-surface conduction, {nx} cells across the slab  (plan 05, gate 5.3)
;
;  Generated by ../generate.py -- edit that, not this.
;
;  The mesh translates rigidly at v = {vel:.9E} m/s, so the problem is steady in
;  the mesh frame and the oracle is the closed-form profile
;
;      T(xi) = Tinf + (Tw - Tinf) [exp(-xi/d) - exp(-L/d)] / [1 - exp(-L/d)]
;
;  with d = alpha/v = {decay} m. See ../verify.py for what is measured and why
;  the expected order is ONE, not two: the ALE remap upwinds on the sign of the
;  swept volume, which under rigid translation is exactly first-order upwind
;  advection.
; =============================================================================

[FUSS-Parameters]
newrun = true
res-threshold = 1e-30
time-threshold = {t_end}
iter-threshold = 100000000

[FUSS-Numerics]
time-scheme = RK3
vnn = {vnn}
time-accurate = true
integration-variables = cons
irs = false

[FUSS-MeshMotion]
law = translation
vel = {vel:.9E} 0.0 0.0
gcl-tolerance = 1.0e-10

[FUSS-IO]
sol-diter     = 100000000
sol-overwrite = true
sol-format    = tecplot ascii
ini-diter     = 100000000

;

[GRIB-meshgen]
outpath = MESH/
method = fmsh
type = 2Dplane
width = {wz}

[GRIB-Block1]
pivot = 0.0 0.0 0.0
length = {lx}
height = {hy}
nx = {nx}
ny = {ny}
strx = U 1.0
stry = U 1.0

;

[GPB-Phase1]
type = solid
material = mat
cp  = {cp}
rho = {rho}
k   = {k}
Tmax = 2000

;

[ICB-Block1]
type = homogeneous
material = mat
T = {tinf}

;

[BCB-Block1]
face1 = Tsurface
face2 = Tfarfield
face3 = Wadiabatic
face4 = Wadiabatic
face5 = null
face6 = null

[Tsurface]
type = wall
T = {tw}

[Tfarfield]
type = wall
T = {tinf}

[Wadiabatic]
type = wall
q = 0.0
"""


def main():
    for nx in GRIDS:
        case = os.path.join(HERE, "N%03d" % nx)
        os.makedirs(os.path.join(case, "INPUT"), exist_ok=True)
        os.makedirs(os.path.join(case, "MESH"), exist_ok=True)

        write_ic(os.path.join(case, "INPUT", "ic.tec"), nx, "B1-SP", True)
        write_ic(os.path.join(case, "MESH", "mesh.tec"), nx, "Block1", False)
        write_bc(os.path.join(case, "INPUT", "bc.txt"), nx)

        # Same material as every other test case, so the property table is
        # identical -- copied rather than regenerated because building it needs
        # ATLAS/GPB, which is not part of this repository.
        for name in ("phase.txt", "properties.dat"):
            shutil.copy(os.path.join(DONOR, "INPUT", name),
                        os.path.join(case, "INPUT", name))

        with open(os.path.join(case, "input.ini"), "w") as fh:
            fh.write(INI.format(nx=nx, ny=NY, lx=LX, hy=HY, wz=WZ, cp=CP,
                                rho=RHO, k=K, tw=TW, tinf=TINF, vel=VEL,
                                decay=DECAY, t_end=T_END, vnn=VNN))

        # FUSS.sh with MASTERDIR fixed up for the extra directory level.
        with open(os.path.join(DONOR, "FUSS.sh")) as fh:
            sh = fh.read().replace("MASTERDIR=../../..", "MASTERDIR=../../../..")
        dst = os.path.join(case, "FUSS.sh")
        with open(dst, "w") as fh:
            fh.write(sh)
        os.chmod(dst, os.stat(dst).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

        dx = LX / nx
        print("N%03d  dx=%.6e  Pe=dx/d=%.5f" % (nx, dx, dx / DECAY))

    print("alpha = %.9e m2/s   v = %.9e m/s   d = alpha/v = %g m"
          % (ALPHA, VEL, DECAY))


if __name__ == "__main__":
    main()
