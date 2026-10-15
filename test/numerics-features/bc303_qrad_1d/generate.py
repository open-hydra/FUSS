#!/usr/bin/env python3
"""
Generate the one-dimensional BC 303 case with a closed form (plan 10, D14).

WHAT THIS CASE IS FOR
---------------------
BC 303 reads three values, hconv, qrad and Tref, and the dispatch once handed the
wall routine the never-assigned bc%qw instead of qrad, at BOTH call sites: the flux
(Mod_BC_Fluxes) and the wall output (IO_Wall), which recomputes Tw and qw with a
second call. The companion case ../bc303_qrad checks the BC's own identity
qw - hconv (Tref - Tw) = qrad on the wall output, so it sees the wall-output site
only: with the flux site reverted, the interior moves by tens of K and that
identity still closes. This case has a closed form, so it sees both:

  * the interior temperature (the flux site),
  * the wall temperature and flux of the wall output (the output site).

THE PROBLEM
-----------
A slab of length L, constant conductivity k, steady, 1-D: face 1 (x = 0) held at
T0 by BC 302, face 2 (x = L) BC 303 with hconv h, radiative input qrad and
reference temperature Tref; the side faces adiabatic. The steady profile is
linear,

    T_w = (k T0 / L + h Tref + qrad) / (k / L + h)       = 700 K
    q_w = k (T_w - T0) / L = h (Tref - T_w) + qrad       = 4.0e4 W/m2
    T(x) = T0 + (T_w - T0) x / L                          (310 K .. 690 K at the centres)

and it is a discrete steady state of FUSS's scheme, not an approximation: BC 302
and BC 303 both use the half-cell distance (2 k (M.n) with M.n = 1/dx on a
Cartesian cell) and the interior faces the two-point difference, all exact on a
linear profile. The tolerance therefore comes from how far the run converges
(res-threshold), not from a discretisation argument; see verify.py.

Run this from the directory it lives in; it rewrites INPUT/, MESH/, input.ini and
FUSS.sh.
"""

import os
import stat

HERE = os.path.dirname(os.path.abspath(__file__))
DONOR = os.path.abspath(os.path.join(HERE, "..", "..", "moving-mesh", "freestream"))

# ------------------------------------------------------------------ geometry
L, NX = 0.10, 20            # dx = 5.0e-3
HY, NY = 0.02, 4
WZ = 0.01

# ------------------------------------------------------------------ material
CP, RHO, K = 1000.0, 1000.0, 10.0
TMAX, NTAB = 2000, 2000

# --------------------------------------------------------------------- state
T0 = 300.0                          # BC 302 on face 1
HCONV, QRAD, TREF = 100.0, 1.0e4, 1000.0    # BC 303 on face 2
T_INIT = 300.0

DX = L / NX
T_W = (K * T0 / L + HCONV * TREF + QRAD) / (K / L + HCONV)
Q_W = K * (T_W - T0) / L


def t_exact(x):
    """Steady profile, x from the 302 face."""
    return T0 + (T_W - T0) * x / L


def fmt(v):
    return " %.15E" % v


def write_tec(path, with_vars):
    xs = [L * i / NX for i in range(NX + 1)]
    ys = [HY * j / NY for j in range(NY + 1)]
    zs = [0.0, WZ]
    ni, nj, nk = len(xs), len(ys), len(zs)
    if with_vars:
        lines = [' VARIABLES ="x" "y" "z" T matID',
                 " ZONE  T = B1-SP, I=%d, J=%d, K=%d, DATAPACKING=BLOCK, "
                 "VARLOCATION=([1-3]=NODAL,[4-5]=CELLCENTERED), "
                 "SOLUTIONTIME=0.000000000000000E+000, STRANDID = 0" % (ni, nj, nk)]
    else:
        lines = [' VARIABLES ="x" "y" "z"',
                 " ZONE  T = Block1, I=%d, J=%d, K=%d, DATAPACKING=BLOCK, "
                 "VARLOCATION=([1-3]=NODAL), "
                 "SOLUTIONTIME=0.000000000000000E+000, STRANDID = 0" % (ni, nj, nk)]
    for axis in range(3):
        for k in range(nk):
            for j in range(nj):
                for i in range(ni):
                    lines.append(fmt((xs[i], ys[j], zs[k])[axis]))
    if with_vars:
        # uniform start, NOT the steady profile: the run has to converge onto it
        ncell = (ni - 1) * (nj - 1) * (nk - 1)
        lines += [fmt(T_INIT)] * ncell
        lines += [fmt(1.0)] * ncell
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def write_properties(path):
    lines = ['TITLE = "Mass Thermodynamic Properties"',
             'VARIABLES = "Temperature", "Cp", "Density", "Conductivity", "Energy"',
             'ZONE T="slab"', "I=%d, F=POINT" % NTAB]
    for i in range(1, NTAB + 1):
        t = float(i)
        lines.append("%.1f %f %f %f %f" % (t, CP, RHO, K, RHO * CP * t))
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def write_bc(path):
    """Block-major, face-minor records, the order Mod_BC_Fluxes groups by."""
    out = []

    def rec(i, j, k, f, typ, second=None):
        out.append("%8d%8d%8d%8d%8d%8d" % (1, i, j, k, f, typ))
        if second is not None:
            out.append(second)

    def val(*v):
        return "".join("    %13.6E," % x for x in v)

    for j in range(1, NY + 1):                  # face 1: Dirichlet T0
        rec(1, j, 1, 1, 302, val(T0))
    for j in range(1, NY + 1):                  # face 2: hconv, qrad, Tref
        rec(NX, j, 1, 2, 303, val(HCONV, QRAD, TREF))
    for i in range(1, NX + 1):                  # face 3
        rec(i, 1, 1, 3, 301, val(0.0))
    for i in range(1, NX + 1):                  # face 4
        rec(i, NY, 1, 4, 301, val(0.0))
    for f in (5, 6):
        for j in range(1, NY + 1):
            for i in range(1, NX + 1):
                rec(i, j, 1, f, 0)
    with open(path, "w") as fh:
        fh.write("\n".join(out) + "\n")


INI = """\
; =============================================================================
;  One-dimensional slab with a closed form for BC 303 (plan 10, D14)
;
;  Generated by generate.py -- edit that, not this.
;
;  Face 1 at T0 = {t0:g} K (302), face 2 BC 303 with hconv = {h:g}, qrad = {qrad:g},
;  Tref = {tref:g}; k = {k:g} W/m/K, L = {l:g} m. Steady and linear:
;
;      T_w = (k T0/L + h Tref + qrad) / (k/L + h) = {tw:.6f} K
;      q_w = k (T_w - T0)/L                       = {qw:.6f} W/m2
;
;  verify.py checks the interior against T(x) (the flux site of the 303 dispatch)
;  and the wall output's Tw, qw (its output site).
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
length = {l}
height = {hy}
nx = {nx}
ny = {ny}
strx = U 1.0
stry = U 1.0

;

[GPB-Phase1]
type = solid
material = slab
cp  = {cp}
rho = {rho}
k   = {k}
Tmax = {tmax}

;

[ICB-Block1]
type = homogeneous
material = slab
T = {tinit}

;

[BCB-Block1]
face1 = Twall
face2 = convection
face3 = adiabatic
face4 = adiabatic
face5 = null
face6 = null

[Twall]
type = wall
T = {t0}

[convection]
type = wall
hconv = {h}
Tref  = {tref}
qrad  = {qrad}

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
        fh.write("solid phase\nslab 1\n")
    with open(os.path.join(HERE, "input.ini"), "w") as fh:
        fh.write(INI.format(t0=T0, h=HCONV, qrad=QRAD, tref=TREF, k=K, l=L, tw=T_W, qw=Q_W,
                            wz=WZ, hy=HY, nx=NX, ny=NY, cp=CP, rho=RHO, tmax=TMAX, tinit=T_INIT))
    with open(os.path.join(DONOR, "FUSS.sh")) as fh:
        sh = fh.read()
    dst = os.path.join(HERE, "FUSS.sh")
    with open(dst, "w") as fh:
        fh.write(sh)
    os.chmod(dst, os.stat(dst).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    print("dx = %.4e   T_w = %.6f K   q_w = %.6f W/m2" % (DX, T_W, Q_W))


if __name__ == "__main__":
    main()
