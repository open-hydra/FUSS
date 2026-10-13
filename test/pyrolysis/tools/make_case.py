#!/usr/bin/env python3
"""
make_case.py -- write a complete FUSS pyrolysis case directory.

Implements plan 06 (phase 2, pyrolysis on a static mesh) section 4: the four case kinds of
gates 4.2 (isothermal kinetics), 4.3 (energy and mass audit), 4.4 (ATC1 against FIAT) and
4.5 (ATC2-chemistryOff against a local PATO run), plus two verification kinds with exact
oracles (conduction: the slab series; adiabatic0d: the zero-dimensional decomposition ODE),
and plan 08's ablating-wall cases (BC 306: atc2x306, and atc2off --bc306-chemoff),
in the standard test-case layout:

    <case>/FUSS.sh            copy of test/moving-mesh/freestream/FUSS.sh, MASTERDIR fixed for depth
    <case>/input.ini          FUSS sections only (no ATLAS sections), probes by index-position
    <case>/Tw_time.dat        BC 302 time file, at the case ROOT (resolved relative to the run dir)
    <case>/probes.json        probe name, i, cell-centre x (atc kinds and conduction; verify.py
                              interpolates in x), n_cells, the grading r and delta_1 (r = 1 when
                              --uniform)
    <case>/MESH/mesh.tec      nodes only, zone Block1
    <case>/INPUT/ic.tec       zone B1-SP, cell-centred T and matID (= 1.0)
    <case>/INPUT/bc.txt       records in face order 1..6, exactly 2 nj nk + 2 ni nk + 2 nj ni of them
    <case>/INPUT/phase.txt    `solid phase` then `<material> 1`
    <case>/INPUT/properties.dat   the LEGACY table Load_Table still reads: T cp rho k h, T = 1..4000
    <case>/INPUT/{pyrolysis,virgin,char,gas}.dat   copied from --tables (make_tacot_tables.py output)

Geometry common to every kind: one block, cells along i; j: 2 cells of 0.1 m; k: 1 cell of 0.1 m.
ic.tec / mesh.tec / bc.txt use the formats of test/moving-mesh/receding/generate.py.

Usage:
    make_case.py <case-dir> --kind {isothermal,energy_audit,atc1,atc2off,atc2x306,conduction,adiabatic0d}
                 --tables <dir> [--T0 K] [--vnn V] [--no-kinetics] [--ncells N] [--uniform]
                 [--fuss-sh <path>] [--depth N] [--case 2.1] [--lambda L] [--bprime FILE]
                 [--bc306-chemoff]

    isothermal    4 cells x 1 mm, both i-faces BC 302 at T0, T(0) = T0, 30 s, fields every 1 s
                  (--T0 one of 700 / 900 / 1200; tables: the --isothermal set)
    energy_audit  40 cells x 1 mm, face 2 BC 302 ramp 300 -> 1644 K in 1 s, face 1 adiabatic,
                  T(0) = 300 K, 10 s, probes at i = 40, 39, 37, 33 (real TACOT tables)
    atc1          L = 0.05 m, N cells graded 0.1 toward the heated face 2 (BC 302 time file),
                  T(0) = 298 K, 60 s, bracketing probes at 1, 2, 4, 8, 16 mm plus the back face
    atc2off       as atc1 with T(0) = 300 K and face 2 = BC 303 `300.0, 0.0, 2800.0`,
                  bracketing probes at 1, 2, 4, 8, 12, 16, 24 mm plus the back face
    atc2x306      plan 08 section 3.3: atc2off's mesh, tables and probes under the type-2 closure
                  (always; PATO's model), face 2 = the ablating wall BC 306 with line 2
                  `<lambda>, 300.0, bprime_tacot_earth.dat, env_atc2<x>.dat`: the B' table copied from
                  --bprime, the env table = PATO's AblationTestCase_2.x BoundaryConditions completed with
                  the Ta constants (qRad 0, hconv 0, Tedge 300), 120 s, a field file every 0.1 s (the
                  fronts are measured from them); --case 2.1 = static mesh with no-char = true (PATO's
                  movingMesh no); --case 2.2 and 2.3 = the receding mesh, law = nodeshift (PATO's
                  movingMesh yes; 2.3 has h_r = 2.5e7 on the heating rows), and no index probes: a
                  probe is a cell, and every cell moves, so verify.py samples the field files.
    conduction    pure conduction on the pyrolysis path: mesh as atc1 (or N equal cells with
                  --uniform), the kinetics ALWAYS stripped from the copied pyrolysis.dat, face 2
                  stepped to 1644 K at t = 0 (BC 302 time file `0 1644` / `60 1644`), face 1
                  adiabatic, T(0) = 298 K, 60 s, the probes of atc1; oracle = the exact slab
                  series (use the --isothermal tables: constant cp and k)
    adiabatic0d   4 cells x 1 mm, both i-faces adiabatic (BC 301, q = 0), T(0) = T0 (--T0
                  required), 20 s with a field file every 0.5 s, no probes, tables unchanged
                  (real TACOT); the uniform state follows the type-1 decomposition ODE. With
                  --qvol Q, INPUT/st.tec carries the uniform volumetric source Q [W/m3] and the
                  run is 3 s with a field file every 0.05 s (--tend / --sol-dtime override)

    --no-kinetics rewrites the copied pyrolysis.dat with nPyroReac[i] 0 for every phase and the
                  reaction lines removed (gate 4.3 (v)); implied by --kind conduction.
    --ncells N    cells along i for atc1, atc2off and conduction (default 100): graded with the
                  same last/first ratio 0.1, r = 0.1^(1/(N-1)), delta_1 = L (1-r)/(1-r^N), the
                  last node exactly L; the bracketing probes are recomputed from the cell centres.
    --uniform     conduction only: N equal cells instead of the grading.
    --qvol Q      adiabatic0d only, Q > 0: INPUT/st.tec = `VARIABLES ="x" "y" "z" qvol`, the ic.tec
                  ZONE line with [4]=CELLCENTERED, the three nodal blocks of ic.tec, one %.15E
                  value of Q per cell (3 nn + nc + 2 lines; Read_vtk_tec reads it beside ic.tec).
    --tend T      override the kind's time-threshold [s].
    --sol-dtime D override the kind's sol-dtime [s].
    --depth N     MASTERDIR depth of the case below the repo root (default: from the case path when
                  it lies under <repo>/test/, else the kind's canonical path).
    --type2       plan 07: gas-closure = type2 with a [FUSS-Darcy] block (solver tridiag, p-init = p-wall
                  = 101325); the gas leaves through the 301 faces of adiabatic0d, the 302 faces of the
                  other kinds, atc2off's 303 by default; probes record T and p_pg.
    --permeability K  isotropic K_v = K_c = K [m2] into the copied pyrolysis.dat (the 0-D type-2 gates use
                  1e-7 so that the fixed-p oracle applies: the pressure excursion is then < 1 Pa).
    --case C      atc2x306 only: the Ablation Test Case variant (2.1, 2.2 or 2.3; default 2.1).
    --lambda L    atc2x306 only: the blowing-correction lambda of the 306 line (default 0.5; 0 = no
                  blowing correction, plan 08 gate 3.1 (b)).
    --bprime F    atc2x306 only: the B' table copied to INPUT/bprime_tacot_earth.dat (default
                  test/pyrolysis/ablation/tables/bprime_tacot_earth.dat).
    --heat-until T  atc2x306 only: the switch-off row of the env table at T instead of 60 s (the
                  chemistry ramps off over T .. T + 0.1 s as PATO's does over 60 .. 60.1 s); the short
                  closure cases use 1 s so that 2 s cover chemistry on, the ramp and chemistry off.
    --sol-overwrite {true,false}  override the kind's sol-overwrite.
    --bc306-chemoff  atc2off only (plan 08 gate 3.1 (a)): face 2 = BC 306 `0.5, 300.0, none, env_off.dat`
                  with INPUT/env_off.dat = two chemistry-off rows carrying the 303 line (hconv 300,
                  T_edge 2800, p_w 101325), and the emissivity and absorptivity columns of the copied
                  virgin.dat and char.dat set to 0, so that the 306 routine hands the face to the BC 303
                  kernel and every output but the wall file matches the unmodified case byte for byte.

Python 3.6, no third-party modules.
"""

import argparse
import json
import math
import os
import re
import shutil
import stat
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
DONOR_FUSS_SH = os.path.join(REPO, "test", "moving-mesh", "freestream", "FUSS.sh")

# ---------------------------------------------------------------- common geometry
NY = 2            # transverse cells (ndir = 2 like every existing case)
DY = 0.1          # m
NK = 1
DZ = 0.1          # m

# ---------------------------------------------------------------- virgin TACOT diffusivity for the dt estimate
ALPHA_V = 0.4025 / (280.0 * 984.0)

TABLE_FILES = ("pyrolysis.dat", "virgin.dat", "char.dat", "gas.dat")

KINDS = {
    "isothermal": {"canonical": "test/pyrolysis/kinetics_isothermal/T%04d", "T0": None,
                   "t_end": 30.0, "sol_dtime": 1.0, "sol_overwrite": "false"},
    "energy_audit": {"canonical": "test/pyrolysis/energy_audit", "T0": 300.0,
                     "t_end": 10.0, "sol_dtime": 10.0, "sol_overwrite": "true"},
    "atc1": {"canonical": "test/pyrolysis/ATC1", "T0": 298.0,
             "t_end": 60.0, "sol_dtime": 10.0, "sol_overwrite": "false"},
    "atc2off": {"canonical": "test/pyrolysis/ATC2-chemistryOff", "T0": 300.0,
                "t_end": 60.0, "sol_dtime": 10.0, "sol_overwrite": "false"},
    "atc2x306": {"canonical": "test/pyrolysis/ablation/ATC21", "T0": 300.0,
                 "t_end": 120.0, "sol_dtime": 0.1, "sol_overwrite": "false"},
    "conduction": {"canonical": "test/pyrolysis/conduction", "T0": 298.0,
                   "t_end": 60.0, "sol_dtime": 10.0, "sol_overwrite": "false"},
    "adiabatic0d": {"canonical": "test/pyrolysis/adiabatic0d/T%04d", "T0": None,
                    "t_end": 20.0, "sol_dtime": 0.5, "sol_overwrite": "false"},
}

ATC_L = 0.05
ATC_N = 100              # default --ncells
ATC_RATIO = 0.1          # last/first cell, PATO blockMeshDict simpleGrading (1 0.1 1)
ATC1_DEPTHS_MM = [1, 2, 4, 8, 16]
ATC2_DEPTHS_MM = [1, 2, 4, 8, 12, 16, 24]
GRADED_KINDS = ("atc1", "atc2off", "atc2x306", "conduction")     # L = ATC_L, --ncells, bracketing probes, probes.json

TW_TIME = {
    "energy_audit": [(0.0, 300.0), (1.0, 1644.0), (60.0, 1644.0)],   # 1 s ramp: resolved by the coarsest run (dt ~ 0.17 s); a 0.1 s ramp aliased gate 4.3 (iv), measured 2026-10-07
    "atc1": [(0.0, 298.0), (0.1, 1644.0), (60.0, 1644.0), (60.1, 298.0), (120.0, 298.0)],
    "conduction": [(0.0, 1644.0), (60.0, 1644.0)],   # a step at t = 0: update() returns var(1) at t = time(1) and clamps beyond the last pair
}

BC303_LINE = (300.0, 0.0, 2800.0)      # hconv, qrad, Tref

# ---------------------------------------------------------------- plan 08: the ablating wall (BC 306)
BPRIME_NAME = "bprime_tacot_earth.dat"
DEFAULT_BPRIME = os.path.join(REPO, "test", "pyrolysis", "ablation", "tables", BPRIME_NAME)
TINF = 300.0                           # T_inf of the 306 line: PATO's Tbackground
# The env table of a 306 face, columns t p_w rhoeUeCH0 h_r q_rad chem hconv T_edge: PATO's
# tutorials/1D/AblationTestCase_2.x/constant/porousMat/BoundaryConditions (t p rhoUeCH h_r chemistryOn)
# completed with the constants of its origin.0/porousMat/Ta (qRad 0, hconv 0, Tedge 300). 2.3 is 2.2
# with h_r = 2.5e7 on the heating rows (ATC2 v2.8, no shipped PATO variant).
def env_atc2x(case, t_off=60.0):
    """t_off moves the switch-off (60 s in PATO) for the short closure cases; the rows keep PATO's shape."""
    hr = 2.5e7 if case == "2.3" else 1.5e6
    return [(0.0, 101325.0, 0.003, 0.0, 0.0, 1.0, 0.0, 300.0),
            (0.1, 101325.0, 0.3, hr, 0.0, 1.0, 0.0, 300.0),
            (t_off, 101325.0, 0.3, hr, 0.0, 1.0, 0.0, 300.0),
            (round(t_off + 0.1, 12), 101325.0, 0.003, 0.0, 0.0, 0.0, 0.0, 300.0),
            (max(120.0, 2.0 * t_off), 101325.0, 0.003, 0.0, 0.0, 0.0, 0.0, 300.0)]
# gate 3.1 (a): the 303 line as two chemistry-off rows (hconv 300, T_edge 2800, p_w the type-2 p-wall)
ENV_OFF = [(0.0, 101325.0, 0.0, 0.0, 0.0, 0.0, BC303_LINE[0], BC303_LINE[2]),
           (60.0, 101325.0, 0.0, 0.0, 0.0, 0.0, BC303_LINE[0], BC303_LINE[2])]
ENV_HEADER = "// t p_w rhoeUeCH0 h_r q_rad chem hconv T_edge   (SI; plan 08 section 2.1, linear in t, clamped)"


def die(msg):
    sys.stderr.write("[ERROR] make_case: %s\n" % msg)
    sys.exit(1)


def fmt(v):
    return " %.15E" % v


# ---------------------------------------------------------------- mesh / ic (generate.py formats)
def write_ic(path, xs, zone_name, with_vars, T0, qvol=None, ny=NY):
    """ic.tec (with_vars: T and matID per cell), mesh.tec (nodes only) or, with qvol, the volumetric
    source file st.tec: the same nodal blocks and ZONE line with [4]=CELLCENTERED, one qvol per cell
    (Read_vtk_tec reads INPUT/st.tec beside ic.tec and takes its first cell-centred variable as qvol).
    Returns the number of lines written."""
    ni, nj, nk = len(xs), ny + 1, NK + 1
    nx = ni - 1
    ys = [DY * j for j in range(nj)]
    zs = [DZ * k for k in range(nk)]

    lines = []
    if qvol is not None:
        lines.append(' VARIABLES ="x" "y" "z" qvol')
        varloc = "VARLOCATION=([1-3]=NODAL,[4]=CELLCENTERED)"
    elif with_vars:
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
    for coord, axis in ((xs, 0), (ys, 1), (zs, 2)):
        for k in range(nk):
            for j in range(nj):
                for i in range(ni):
                    v = (xs[i], ys[j], zs[k])[axis]
                    lines.append(fmt(v))
    ncell = NK * ny * nx
    if qvol is not None:
        for _ in range(ncell):
            lines.append(fmt(qvol))
    elif with_vars:
        for _ in range(ncell):
            lines.append(fmt(T0))
        for _ in range(ncell):
            lines.append(fmt(1.0))
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return len(lines)


# ---------------------------------------------------------------- bc.txt
VAL_FMT = "    %13.6E,"                                  # generate.py's constant value line
VAL3_FMT = "    %13.6E, %13.6E, %13.6E,"                 # BC 303: hconv, qrad, Tref
TIMEFILE_LINE = " Tw_time.dat,"                          # BC 302/301 time file: the bare name


def write_bc(path, nx, face1, face2, ny=NY):
    """face1/face2: (type, value-line or None). Faces 3/4 = 300, 5/6 = 0 (no value line).

    Record count = 2 nj nk + 2 ni nk + 2 nj ni (Allocate_BC); face order 1..6 as generate.py.
    """
    out = []

    def rec(i, j, k, f, typ, line=None):
        out.append("%8d%8d%8d%8d%8d%8d" % (1, i, j, k, f, typ))
        if line is not None:
            out.append(line)

    for j in range(1, ny + 1):
        for k in range(1, NK + 1):
            rec(1, j, k, 1, face1[0], face1[1])
    for j in range(1, ny + 1):
        for k in range(1, NK + 1):
            rec(nx, j, k, 2, face2[0], face2[1])
    for i in range(1, nx + 1):
        for k in range(1, NK + 1):
            rec(i, 1, k, 3, 300)
    for i in range(1, nx + 1):
        for k in range(1, NK + 1):
            rec(i, ny, k, 4, 300)
    for j in range(1, ny + 1):
        for i in range(1, nx + 1):
            rec(i, j, 1, 5, 0)
    for j in range(1, ny + 1):
        for i in range(1, nx + 1):
            rec(i, j, 1, 6, 0)
    with open(path, "w") as fh:
        fh.write("\n".join(out) + "\n")
    nrec = 2 * ny * NK + 2 * nx * NK + 2 * ny * nx
    return nrec


# ---------------------------------------------------------------- tables
def read_table_zone1(path, nvar):
    """First zone of a FUSS pyrolysis table: {T: [values]}, T integer-valued."""
    rows = {}
    nzone = 0
    with open(path) as fh:
        for line in fh:
            if "ZONE" in line or "Zone" in line:
                nzone += 1
                continue
            if nzone != 1:
                continue
            toks = line.split()
            try:
                vals = [float(t) for t in toks]
            except ValueError:
                continue
            if len(vals) != nvar + 1:
                die("%s: row with %d values, expected %d" % (path, len(vals), nvar + 1))
            rows[int(round(vals[0]))] = vals[1:]
    if not rows:
        die("%s: no rows in zone 1" % path)
    return rows


def read_pyrolysis_dat(path):
    """material name, rho_v = sum eps_I rho_I, nPyroReac list; validates the grammar loosely."""
    name = None
    keys = {}
    with open(path) as fh:
        for ln, line in enumerate(fh, 1):
            s = line.split("//")[0].strip()
            if not s:
                continue
            toks = s.split(None, 1)
            if len(toks) != 2:
                die("%s line %d: expected `key value`: %r" % (path, ln, line.rstrip()))
            if toks[0] not in ("K_v", "K_c", "B_v", "B_c") and len(toks[1].split()) != 1:
                die("%s line %d: expected `key value`: %r" % (path, ln, line.rstrip()))
            if toks[0] == "material":
                if name is not None:
                    die("%s: more than one material block (phase 2 cases carry one material)" % path)
                name = toks[1]
                continue
            if toks[0] in keys:
                die("%s line %d: key %s appears twice" % (path, ln, toks[0]))
            keys[toks[0]] = toks[1]
    if name is None:
        die("%s: no `material <name>` line" % path)
    nph = int(keys["nSolidPhases"])
    rho_v = 0.0
    nreac = []
    for i in range(1, nph + 1):
        rho_v += float(keys["epsI[%d]" % i]) * float(keys["rhoI[%d]" % i])
        nreac.append(int(keys["nPyroReac[%d]" % i]))
    return name, rho_v, nreac


def set_permeability(path, K):
    """--permeability K: isotropic K_v = K_c = K [m2] in the copied pyrolysis.dat (plan 07 section 3.1),
    replacing the K_v / K_c lines when present, appending them otherwise."""
    line = {"K_v": "K_v %g 0 0 0 %g 0 0 0 %g\n" % (K, K, K), "K_c": "K_c %g 0 0 0 %g 0 0 0 %g\n" % (K, K, K)}
    out, seen = [], set()
    with open(path) as fh:
        for s in fh:
            key = s.split("//", 1)[0].split()
            k = key[0] if key else ""
            if k in line:
                out.append(line[k])
                seen.add(k)
            else:
                out.append(s)
    for k in ("K_v", "K_c"):
        if k not in seen:
            out.append(line[k])
    with open(path, "w") as fh:
        fh.write("".join(out))


def strip_kinetics(path):
    """--no-kinetics: nPyroReac[i] 0 for every phase, reaction lines removed."""
    reac_re = re.compile(r"^\s*(F|A|E|m|n|T|h)\[\d+\]\[\d+\]\s")
    npr_re = re.compile(r"^(\s*nPyroReac\[\d+\])\s+\S+(.*)$")
    out = []
    with open(path) as fh:
        for line in fh:
            if reac_re.match(line):
                continue
            m = npr_re.match(line.rstrip("\n"))
            if m:
                line = "%s 0%s\n" % (m.group(1), m.group(2))
            out.append(line)
    out.insert(0, "// --no-kinetics: every nPyroReac[i] set to 0 and the reaction lines removed by make_case.py\n")
    with open(path, "w") as fh:
        fh.write("".join(out))


def write_properties_dat(path, name, virgin, rho_v):
    """Legacy properties.dat: T cp rho k h for T = 1..4000, cp and k clamped below the table,
    h = rho_v * int_0^T cp dT' (trapezoid on the 1 K grid; monotone, > 0 at 1 K)."""
    Ts = sorted(virgin)
    Tlo, Thi = Ts[0], Ts[-1]
    TMAX = 4000

    def cp_at(T):
        return virgin[min(max(T, Tlo), Thi)][0]

    def k_at(T):
        return virgin[min(max(T, Tlo), Thi)][2]

    lines = ['TITLE = "Mass Thermodynamic Properties"',
             'VARIABLES = "Temperature", "Cp", "Density", "Conductivity", "Energy"',
             'ZONE T="%s"' % name,
             "I=%d, F=POINT" % TMAX]
    h = 0.0
    cp_prev = cp_at(0)
    for T in range(1, TMAX + 1):
        cp = cp_at(T)
        h += rho_v * 0.5 * (cp_prev + cp)
        cp_prev = cp
        lines.append("%.1f %f %f %f %f" % (T, cp, rho_v, k_at(T), h))
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


# ---------------------------------------------------------------- grading / probes
def atc_nodes(n):
    """n cells on [0, ATC_L] with last/first = ATC_RATIO (geometric), the last node exactly ATC_L."""
    r = ATC_RATIO ** (1.0 / (n - 1))
    d1 = ATC_L * (1.0 - r) / (1.0 - r ** n)
    xs = [d1 * (1.0 - r ** i) / (1.0 - r) for i in range(n + 1)]
    xs[0] = 0.0
    xs[-1] = ATC_L
    return xs, r, d1


def uniform_nodes(n):
    """n equal cells on [0, ATC_L]; returns (xs, r = 1, delta_1 = L / n) like atc_nodes."""
    xs = [ATC_L * i / n for i in range(n + 1)]
    xs[0] = 0.0
    xs[-1] = ATC_L
    return xs, 1.0, ATC_L / n


def bracketing_probes(xs, depths_mm):
    """Two probes per depth at the cells whose centres bracket x = L - d, plus `back` at i = 1."""
    nx = len(xs) - 1
    xc = [0.5 * (xs[i - 1] + xs[i]) for i in range(1, nx + 1)]     # xc[i-1] = centre of cell i
    probes = []
    pairs = {}
    for d in depths_mm:
        xt = ATC_L - d * 1e-3
        lo = None
        for i in range(1, nx + 1):
            if xc[i - 1] <= xt:
                lo = i
        if lo is None or lo >= nx:
            die("depth %g mm has no bracketing cell pair" % d)
        hi = lo + 1
        # Probe names at most 5 characters: FUSS builds the file name as OUTPUT/<name>.txt in a
        # 16-character field (io_probes_t%file, clen) and silently truncates longer ones (.tx).
        n_lo, n_hi = "d%02dlo" % d, "d%02dhi" % d
        probes.append({"name": n_lo, "i": lo, "x": xc[lo - 1], "depth_mm": d, "role": "lo"})
        probes.append({"name": n_hi, "i": hi, "x": xc[hi - 1], "depth_mm": d, "role": "hi"})
        pairs[str(d)] = [n_lo, n_hi]
    probes.append({"name": "back", "i": 1, "x": xc[0], "depth_mm": ATC_L * 1e3, "role": "back"})
    return probes, pairs, xc


# ---------------------------------------------------------------- input.ini
INI_HEAD = {
    "isothermal": [
        "Isothermal kinetics on the --isothermal test tables  (plan 06, gate 4.2)",
        "",
        "4 cells x 1 mm, both i-faces held at T0 by BC 302, T(0) at T0, 30 s with a field file",
        "every second. With cp, h and the gas h_pg all linear in T and equal, T stays at T0 and the",
        "reaction advancements follow the closed-form propagator; verify.py compares chi and T.",
    ],
    "energy_audit": [
        "Energy and mass audit on the real TACOT tables  (plan 06, gate 4.3)",
        "",
        "40 cells x 1 mm, face 2 ramped 300 -> 1644 K in 1 s through Tw_time.dat (BC 302), face 1",
        "adiabatic (BC 301), 10 s. verify.py closes the energy and solid-mass balances from",
        "energy-budget.txt and the field file, and compares the probe series across vnn.",
    ],
    "atc1": [
        "Ablation Test Case 1 on the real TACOT tables  (plan 06, gate 4.4; oracle FIAT)",
        "",
        "L of 0.05 m, %(ncells)d cells graded 0.1 toward the heated face 2, which follows Tw_time.dat",
        "(298 -> 1644 K in 0.1 s, held to 60 s); face 1 adiabatic; T(0) at 298 K. Two probes bracket",
        "each FIAT thermocouple depth (see probes.json); verify.py interpolates in x and in time.",
    ],
    "atc2off": [
        "Ablation Test Case 2, chemistry off  (plan 06, gate 4.5; oracle a local PATO run)",
        "",
        "Mesh and tables as ATC1; T(0) at 300 K; face 2 is the closed-form Robin wall BC 303 with",
        "hconv 300 W/m2/K, qrad 0, Tref 2800 K from t equal 0 (PATO run with zero emissivity).",
        "Two probes bracket each PATO probe depth (see probes.json).",
    ],
    "atc2off-306": [
        "Ablation Test Case 2, chemistry off, face 2 as the ablating wall BC 306  (plan 08, gate 3.1 (a))",
        "",
        "The atc2off case with two edits: the emissivity and absorptivity columns (7 and 8) of",
        "virgin.dat and char.dat set to 0, and face 2 a BC 306 face whose env table INPUT/env_off.dat",
        "keeps the chemistry off with hconv 300 W/m2/K and T_edge 2800 K. The 306 routine then hands",
        "the face to the BC 303 kernel: every output but the wall file matches the unmodified case.",
    ],
    "atc2x306": [
        "Ablation Test Case %(atc_case)s, chemistry on  (plan 08 section 3.3; oracle a local PATO run)",
        "",
        "Mesh and tables as ATC2-chemistryOff; T(0) at 300 K; face 2 is the ablating wall BC 306, the",
        "surface energy balance with the equilibrium B' table INPUT/bprime_tacot_earth.dat and PATO's",
        "boundary conditions as the time table INPUT/%(env)s; lambda %(lam)s, T_inf 300 K.",
        "%(mesh_note)s",
        "%(probe_note)s %(io_note)s",
    ],
    # No `=` inside a comment line (hydra INI rule: ATLAS's FiNeR reads such a line as an option).
    "conduction": [
        "Pure conduction on the pyrolysis path, constant properties, step-heated face 2 and insulated",
        "face 1; oracle: the exact slab series (verify.py)",
        "",
        "L of 0.05 m, %(ncells)d cells %(grading)s.",
        "Face 2 held at 1644 K from t equal 0 through Tw_time.dat (BC 302, a step); face 1 adiabatic",
        "(BC 301); T(0) at %(T0)g K; the kinetics stripped from pyrolysis.dat (every nPyroReac 0).",
        "Two probes bracket each depth of ATC1 (see probes.json); verify.py interpolates in x and",
        "in time.",
    ],
    "adiabatic0d": [
        "Uniform adiabatic slab: the state follows the exact type-1 ODE, the gas leaving at h_pg(T)",
        "(verify.py)",
        "",
        "4 cells x 1 mm, both i-faces adiabatic (BC 301, q 0), T(0) at %(T0)g K, tables unchanged,",
        "20 s with a field file every 0.5 s, no probes. Nothing varies along i, so every cell",
        "follows the zero-dimensional decomposition ODE that verify.py integrates.",
    ],
    "adiabatic0d-qvol": [
        "Uniform slab heated by a volumetric source q: the state follows the exact type-1 ODE",
        "(verify.py)",
        "",
        "4 cells x 1 mm, both i-faces adiabatic (BC 301, zero flux), T(0) at %(T0)g K, tables",
        "unchanged, qvol %(qvol)g W/m3 in every cell from INPUT/st.tec, %(tend)g s with a field file",
        "every %(sol_dtime)g s, no probes. Nothing varies along i, so every cell follows the",
        "zero-dimensional decomposition ODE with the source that verify.py integrates.",
    ],
}


def fmt_time(v):
    """`%.1f` as the first cases were written (30.0, 10.0, 0.5); the exact repr when that would round (0.05)."""
    s = "%.1f" % v
    return s if float(s) == v else repr(v)


def write_ini(path, head_key, cfg, vnn, probes, subst, type2=False, open_codes=None, extra=None):
    head = ["; " + "=" * 77]
    for s in INI_HEAD[head_key]:
        if "%(" in s:
            s = s % subst
        head.append(";  " + s if s else ";")
    if type2:
        head.append(";")
        head.append(";  Type-2 gas closure (plan 07): the gas pressure is solved ([FUSS-Darcy]), the gas storage")
        head.append(";  eps_pg (rho_g h_g - p) is in E, the gas leaves through the open faces at p-wall.")
    head.append(";")
    head.append(";  Generated by test/pyrolysis/tools/make_case.py -- edit that, not this.")
    head.append("; " + "=" * 77)
    if type2:
        pyro = ["[FUSS-Pyrolysis]", "gas-closure = type2", "", "[FUSS-Darcy]", "solver = tridiag",
                "p-init = 101325.0", "p-wall = 101325.0"]
        if open_codes:
            pyro.append("open-bc-types = %s" % open_codes)
    else:
        pyro = ["[FUSS-Pyrolysis]", "gas-closure = type1", "outflow-direction = +i", "p-gas = 101325.0"]
    body = head + [
        "",
        "[FUSS-Parameters]",
        "newrun = true",
        "res-threshold = 1e-30",
        "time-threshold = %s" % fmt_time(cfg["t_end"]),
        "iter-threshold = 100000000",
        "",
        "[FUSS-Numerics]",
        "time-scheme = RK3",
        "vnn = %s" % repr(vnn),
        "time-accurate = true",
        "integration-variables = cons",
        "irs = false",
        "",
    ] + pyro + (extra or []) + [
        "",
        "[FUSS-IO]",
        "sol-dtime     = %s" % fmt_time(cfg["sol_dtime"]),
        "sol-overwrite = %s" % cfg["sol_overwrite"],
        "sol-format    = tecplot ascii",
        "ini-diter     = 100000000",
    ]
    if probes:
        body += ["", "[FUSS-Probes]"]
        for n, p in enumerate(probes, 1):
            body.append("probe%d = %s" % (n, p["name"]))
        for p in probes:
            body += ["", "[%s]" % p["name"],
                     "index-position = 1 %d 1 1" % p["i"],
                     "variables = T p_pg" if type2 else "variables = T",
                     "dtime = 0.1"]
    with open(path, "w") as fh:
        fh.write("\n".join(body) + "\n")


def write_env_file(path, rows):
    """The env table of a 306 face (plan 08 section 2.1): one `//` header, then the rows (repr: exact)."""
    with open(path, "w") as fh:
        fh.write(ENV_HEADER + "\n")
        for r in rows:
            fh.write(" ".join(repr(float(v)) for v in r) + "\n")


def zero_optics(path):
    """Set the emissivity and absorptivity columns (the last two of 8) of a FUSS pyrolysis table to 0,
    every other byte unchanged."""
    out = []
    n = 0
    with open(path) as fh:
        for line in fh:
            toks = line.split()
            try:
                [float(t) for t in toks]
                is_row = len(toks) == 8
            except ValueError:
                is_row = False
            if is_row:
                k = line.rstrip("\n")
                parts = k.split(" ")
                if len(parts) != 8:
                    die("%s: a data row not separated by single blanks: %r" % (path, k))
                parts[6] = parts[7] = "0.000000000E+00"
                line = " ".join(parts) + "\n"
                n += 1
            out.append(line)
    if n == 0:
        die("%s: no 8-column rows" % path)
    with open(path, "w") as fh:
        fh.write("".join(out))
    return n


def write_time_file(path, pairs):
    with open(path, "w") as fh:
        fh.write("\n".join("%g %g" % (t, T) for t, T in pairs) + "\n")


def write_fuss_sh(dst, donor, depth):
    with open(donor) as fh:
        sh = fh.read()
    old = "MASTERDIR=../../.."
    if old not in sh:
        die("%s has no `%s` line to rewrite" % (donor, old))
    new = "MASTERDIR=" + "/".join([".."] * depth)
    with open(dst, "w") as fh:
        fh.write(sh.replace(old, new, 1))
    os.chmod(dst, os.stat(dst).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return new


def case_depth(case_dir, kind, T0, override):
    if override is not None:
        return override, "--depth"
    real = os.path.realpath(case_dir)
    test_root = os.path.join(os.path.realpath(REPO), "test") + os.sep
    if real.startswith(test_root):
        rel = os.path.relpath(real, os.path.realpath(REPO))
        return len(rel.split(os.sep)), "path " + rel
    canon = KINDS[kind]["canonical"]
    if "%" in canon:
        canon = canon % int(round(T0))
    return len(canon.split("/")), "canonical path " + canon


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case_dir")
    ap.add_argument("--kind", required=True, choices=sorted(KINDS))
    ap.add_argument("--tables", required=True, help="directory holding pyrolysis.dat virgin.dat char.dat gas.dat")
    ap.add_argument("--T0", type=float, default=None, help="initial temperature [K] (isothermal: 700/900/1200; adiabatic0d: required)")
    ap.add_argument("--vnn", type=float, default=0.5)
    ap.add_argument("--no-kinetics", action="store_true")
    ap.add_argument("--ncells", type=int, default=None, help="cells along i for atc1, atc2off and conduction (default %d)" % ATC_N)
    ap.add_argument("--uniform", action="store_true", help="conduction only: equal cells instead of the 0.1 grading")
    ap.add_argument("--qvol", type=float, default=None,
                    help="adiabatic0d only: uniform volumetric source Q > 0 [W/m3] written to INPUT/st.tec (then 3 s, a field file every 0.05 s)")
    ap.add_argument("--tend", type=float, default=None, help="override the kind's time-threshold [s]")
    ap.add_argument("--sol-dtime", type=float, default=None, help="override the kind's sol-dtime [s]")
    ap.add_argument("--fuss-sh", default=DONOR_FUSS_SH, help="FUSS.sh to copy (default test/moving-mesh/freestream/FUSS.sh)")
    ap.add_argument("--depth", type=int, default=None, help="MASTERDIR depth below the repo root")
    ap.add_argument("--type2", action="store_true",
                    help="plan 07: gas-closure = type2 with a [FUSS-Darcy] block (p-init = p-wall = 101325; the 301 faces of"
                         " adiabatic0d and the 302 faces of atc1/conduction/energy_audit/isothermal opened, atc2off's 303 by"
                         " default), probes record T and p_pg")
    ap.add_argument("--permeability", type=float, default=None,
                    help="isotropic K_v = K_c [m2] written into the copied pyrolysis.dat (replaces the tables' tensors)")
    ap.add_argument("--case", dest="atc_case", default=None, choices=("2.1", "2.2", "2.3"),
                    help="atc2x306 only: the Ablation Test Case variant (default 2.1)")
    ap.add_argument("--lambda", dest="lam", type=float, default=None,
                    help="atc2x306 only: blowing-correction lambda on the 306 line (default 0.5)")
    ap.add_argument("--bprime", default=None, help="atc2x306 only: the B' table to copy (default %s)" % DEFAULT_BPRIME)
    ap.add_argument("--heat-until", type=float, default=None, help="atc2x306 only: the switch-off time of the env table (default 60)")
    ap.add_argument("--sol-overwrite", choices=("true", "false"), default=None, help="override the kind's sol-overwrite")
    ap.add_argument("--bc306-chemoff", action="store_true",
                    help="atc2off only: face 2 = BC 306 with a chemistry-off env table carrying the 303 line, optics zeroed")
    args = ap.parse_args()

    kind = args.kind
    cfg = dict(KINDS[kind])                    # the kind's defaults; --qvol / --tend / --sol-dtime adjust the copy
    if kind != "atc2x306" and (args.atc_case or args.lam is not None or args.bprime or args.heat_until is not None):
        die("--case, --lambda, --bprime and --heat-until apply to --kind atc2x306 only")
    if args.sol_overwrite is not None:
        cfg["sol_overwrite"] = args.sol_overwrite
    t_off = 60.0 if args.heat_until is None else args.heat_until
    if t_off <= 0.1:
        die("--heat-until must be > 0.1 (the heating ramp ends at 0.1 s)")
    if args.bc306_chemoff and kind != "atc2off":
        die("--bc306-chemoff applies to --kind atc2off only")
    atc_case = args.atc_case or "2.1"
    lam = 0.5 if args.lam is None else args.lam
    bprime = args.bprime or DEFAULT_BPRIME
    if kind == "atc2x306":
        if atc_case not in ("2.1", "2.2", "2.3"):
            die("--case must be 2.1, 2.2 or 2.3")
        if lam < 0:
            die("--lambda must be >= 0")
        if not os.path.isfile(bprime):
            die("B' table not found: %s" % bprime)
        args.type2 = True                      # PATO's model: the gas pressure is solved
    receding = kind == "atc2x306" and atc_case != "2.1"      # plan 08 stage C: law nodeshift
    T0 = args.T0 if args.T0 is not None else cfg["T0"]
    if T0 is None:
        die("--T0 is required for --kind %s%s" % (kind, " (700, 900 or 1200)" if kind == "isothermal" else ""))
    if kind == "isothermal" and int(round(T0)) not in (700, 900, 1200):
        print("[NOTE] --T0 %g is not one of the gate-4.2 values 700/900/1200" % T0)
    if args.vnn <= 0:
        die("--vnn must be > 0")
    if args.ncells is not None and kind not in GRADED_KINDS:
        die("--ncells applies to the kinds %s only" % ", ".join(GRADED_KINDS))
    ncells = args.ncells if args.ncells is not None else ATC_N
    if ncells < 2:
        die("--ncells must be at least 2")
    if args.uniform and kind != "conduction":
        die("--uniform applies to --kind conduction only")
    if args.qvol is not None:
        if kind != "adiabatic0d":
            die("--qvol applies to --kind adiabatic0d only")
        if args.qvol <= 0:
            die("--qvol must be > 0 (omit it for the adiabatic case)")
        cfg["t_end"], cfg["sol_dtime"] = 3.0, 0.05
    if args.tend is not None:
        if args.tend <= 0:
            die("--tend must be > 0")
        cfg["t_end"] = args.tend
    if args.sol_dtime is not None:
        if args.sol_dtime <= 0:
            die("--sol-dtime must be > 0")
        cfg["sol_dtime"] = args.sol_dtime
    strip = args.no_kinetics or kind == "conduction"
    for f in TABLE_FILES:
        if not os.path.isfile(os.path.join(args.tables, f)):
            die("%s not found in --tables %s" % (f, args.tables))
    if not os.path.isfile(args.fuss_sh):
        die("FUSS.sh donor not found: %s" % args.fuss_sh)

    case = os.path.abspath(args.case_dir)
    os.makedirs(os.path.join(case, "INPUT"), exist_ok=True)
    os.makedirs(os.path.join(case, "MESH"), exist_ok=True)

    # tables
    for f in TABLE_FILES:
        shutil.copy(os.path.join(args.tables, f), os.path.join(case, "INPUT", f))
    pyro = os.path.join(case, "INPUT", "pyrolysis.dat")
    if strip:
        strip_kinetics(pyro)
    if args.permeability is not None:
        if args.permeability <= 0:
            die("--permeability must be > 0")
        set_permeability(pyro, args.permeability)
    name, rho_v, nreac = read_pyrolysis_dat(pyro)
    if len(name) > 16:
        die("material name %r longer than 16 characters" % name)
    virgin = read_table_zone1(os.path.join(case, "INPUT", "virgin.dat"), 7)
    if kind == "isothermal":
        if any(abs(v[0] - 1000.0) > 0 or abs(v[2] - 0.4) > 0 for v in virgin.values()):
            die("kind isothermal needs the --isothermal tables (virgin cp 1000, k 0.4); %s is not that set" % args.tables)
    elif kind == "conduction":
        if len(set(v[0] for v in virgin.values())) > 1 or len(set(v[2] for v in virgin.values())) > 1:
            print("[WARNING] kind conduction on tables whose virgin cp or k_i varies with T: the exact slab"
                  " series assumes constant properties (use the --isothermal set, not %s)" % args.tables)
    elif all(abs(v[0] - 1000.0) == 0 for v in virgin.values()):
        print("[WARNING] kind %s on tables whose virgin cp is constant 1000: is %s the real TACOT set?" % (kind, args.tables))

    # mesh
    ratio, d1 = None, None
    if kind in ("isothermal", "adiabatic0d"):
        nx = 4
        xs = [1e-3 * i for i in range(nx + 1)]
    elif kind == "energy_audit":
        nx = 40
        xs = [1e-3 * i for i in range(nx + 1)]
    elif args.uniform:
        nx = ncells
        xs, ratio, d1 = uniform_nodes(nx)
    else:
        nx = ncells
        xs, ratio, d1 = atc_nodes(nx)
    write_ic(os.path.join(case, "INPUT", "ic.tec"), xs, "B1-SP", True, T0)
    write_ic(os.path.join(case, "MESH", "mesh.tec"), xs, "Block1", False, T0)
    nst = None
    if args.qvol is not None:
        nst = write_ic(os.path.join(case, "INPUT", "st.tec"), xs, "B1-SP", False, T0, qvol=args.qvol)

    # boundary conditions
    if kind == "isothermal":
        face1 = (302, VAL_FMT % T0)
        face2 = (302, VAL_FMT % T0)
    elif kind == "adiabatic0d":
        face1 = (301, VAL_FMT % 0.0)
        face2 = (301, VAL_FMT % 0.0)
    elif kind in TW_TIME:                      # energy_audit, atc1, conduction
        face1 = (301, VAL_FMT % 0.0)
        face2 = (302, TIMEFILE_LINE)
        write_time_file(os.path.join(case, "Tw_time.dat"), TW_TIME[kind])
    elif kind == "atc2x306":
        env_name = "env_atc%s.dat" % atc_case.replace(".", "")
        face1 = (301, VAL_FMT % 0.0)
        face2 = (306, "     %s, %s, %s, %s" % (repr(lam), repr(TINF), BPRIME_NAME, env_name))
        write_env_file(os.path.join(case, "INPUT", env_name), env_atc2x(atc_case, t_off))
        shutil.copy(bprime, os.path.join(case, "INPUT", BPRIME_NAME))
    elif args.bc306_chemoff:                   # atc2off, gate 3.1 (a)
        face1 = (301, VAL_FMT % 0.0)
        face2 = (306, "     0.5, %s, none, env_off.dat" % repr(TINF))
        write_env_file(os.path.join(case, "INPUT", "env_off.dat"), ENV_OFF)
        nz = [zero_optics(os.path.join(case, "INPUT", f)) for f in ("virgin.dat", "char.dat")]
    else:                                      # atc2off
        face1 = (301, VAL_FMT % 0.0)
        face2 = (303, VAL3_FMT % BC303_LINE)
    nrec = write_bc(os.path.join(case, "INPUT", "bc.txt"), nx, face1, face2)

    # phase.txt, legacy properties.dat
    with open(os.path.join(case, "INPUT", "phase.txt"), "w") as fh:
        fh.write("solid phase\n%s 1\n" % name)
    write_properties_dat(os.path.join(case, "INPUT", "properties.dat"), name, virgin, rho_v)

    # probes
    probes = []
    pairs = {}
    if kind == "energy_audit":
        probes = [{"name": "p05mm", "i": 40}, {"name": "p15mm", "i": 39},
                  {"name": "p35mm", "i": 37}, {"name": "p75mm", "i": 33}]
        for p in probes:
            p["x"] = 0.5 * (xs[p["i"] - 1] + xs[p["i"]])
    elif kind in GRADED_KINDS and not receding:
        depths = ATC2_DEPTHS_MM if kind in ("atc2off", "atc2x306") else ATC1_DEPTHS_MM
        probes, pairs, xc = bracketing_probes(xs, depths)
        with open(os.path.join(case, "probes.json"), "w") as fh:
            json.dump({"kind": kind, "L": ATC_L, "n_cells": nx, "ratio_last_first": 1.0 if args.uniform else ATC_RATIO,
                       "r": ratio, "delta_1": d1, "delta_n": xs[-1] - xs[-2],
                       "depths_mm": depths, "pairs": pairs,
                       "note": "x is the cell-centre coordinate of cell i (nodes in MESH/mesh.tec); "
                               "depth below the heated face 2 is L - x",
                       "probes": probes}, fh, indent=1, sort_keys=True)
            fh.write("\n")
    for p in probes:
        if len(p["name"]) > 16:
            die("probe name %r longer than 16 characters" % p["name"])

    subst = {"ncells": nx, "T0": T0, "qvol": args.qvol, "tend": cfg["t_end"], "sol_dtime": cfg["sol_dtime"],
             "grading": "of equal size" if args.uniform else "graded 0.1 toward the heated face 2"}
    head_key = "adiabatic0d-qvol" if args.qvol is not None else kind
    extra = None
    if args.bc306_chemoff:
        head_key = "atc2off-306"
    if kind == "atc2x306":
        if receding:
            mesh_note = "Receding mesh, law nodeshift: the char is consumed at the rate of the B' table (PATO with movingMesh yes)."
            probe_note = "No index probes (a probe is a cell and every cell moves): verify.py samples the field files."
            extra = ["", "[FUSS-MeshMotion]", "law = nodeshift"]
        else:
            mesh_note = "Static mesh, no-char true: the char is not consumed (PATO with movingMesh no)."
            probe_note = "Two probes bracket each PATO probe depth (see probes.json)."
            extra = ["", "[FUSS-Ablation]", "no-char = true"]
        subst.update({"atc_case": atc_case, "env": env_name, "lam": repr(lam), "probe_note": probe_note,
                      "mesh_note": mesh_note + ("" if t_off == 60.0 else " Heating until %g s." % t_off),
                      "io_note": "A field file every %g s (the fronts)." % cfg["sol_dtime"] if cfg["sol_overwrite"] == "false"
                                 else "seb.txt rows every %g s." % cfg["sol_dtime"]})
    open_codes = None
    if args.type2:
        # the faces the gas leaves through (plan 07 section 3.4): atc2off's 303 is open by default, a 306
        # face always (plan 08 P-D6)
        open_codes = {"adiabatic0d": "301", "atc2off": None, "atc2x306": None}.get(kind, "302")
    write_ini(os.path.join(case, "input.ini"), head_key, cfg, args.vnn, probes, subst, type2=args.type2,
              open_codes=open_codes, extra=extra)

    depth, how = case_depth(case, kind, T0, args.depth)
    masterdir = write_fuss_sh(os.path.join(case, "FUSS.sh"), args.fuss_sh, depth)

    # summary
    dxs = [xs[i + 1] - xs[i] for i in range(nx)]
    dx_min = min(dxs)
    inv_l2 = 1.0 / dx_min ** 2 + 1.0 / DY ** 2
    dt = 0.5 * args.vnn / (ALPHA_V * inv_l2)
    print("make_case.py: %s -> %s" % (kind, case))
    print("  material %s (rho_v = %g), tables from %s%s" % (name, rho_v, args.tables,
          ", kinetics stripped" if strip else ", nPyroReac = %s" % nreac))
    print("  cells %d x %d x %d, dx = %.6g .. %.6g m, dy = %g, dz = %g; bc.txt records = %d"
          % (nx, NY, NK, dx_min, max(dxs), DY, DZ, nrec))
    if kind in GRADED_KINDS and args.uniform:
        print("  uniform cells, dx = %.6f mm" % (d1 * 1e3))
    elif kind in GRADED_KINDS:
        print("  grading r = %.6f, delta_1 = %.6f mm (back face), delta_%d = %.6f mm (heated face)"
              % (ratio, d1 * 1e3, nx, dxs[-1] * 1e3))
    print("  T0 = %g K, t_end = %g s, vnn = %g, dt estimate (VNN, virgin alpha %.4e) = %.4e s -> ~%d steps"
          % (T0, cfg["t_end"], args.vnn, ALPHA_V, dt, int(cfg["t_end"] / dt)))
    print("  faces: 1 = %d, 2 = %d, 3/4 = 300, 5/6 = 0%s" % (face1[0], face2[0],
          "; Tw_time.dat at the case root" if kind in TW_TIME else ""))
    if face2[0] == 306:
        print("  BC 306 line 2: %s" % face2[1].strip())
    if args.bc306_chemoff:
        print("  emissivity and absorptivity zeroed: %d rows of virgin.dat, %d of char.dat" % tuple(nz))
    if args.type2:
        print("  gas closure: type2 ([FUSS-Darcy], open-bc-types %s%s)" % (open_codes or "default 303 305 306",
              ", K = %g m2 isotropic" % args.permeability if args.permeability is not None else ""))
    if nst is not None:
        print("  volumetric source qvol = %g W/m3 in every cell -> INPUT/st.tec (%d lines: 2 + 3 x %d nodes + %d cells)"
              % (args.qvol, nst, (nx + 1) * (NY + 1) * (NK + 1), nx * NY * NK))
    if probes:
        print("  probes (name, i, x [m], depth below face 2 [mm]):")
        for p in probes:
            print("    %-8s i=%3d  x=%.6e  d=%.4f" % (p["name"], p["i"], p["x"], (xs[-1] - p["x"]) * 1e3))
    print("  FUSS.sh: %s (depth %d from %s)" % (masterdir, depth, how))


if __name__ == "__main__":
    main()
