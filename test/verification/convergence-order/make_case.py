#!/usr/bin/env python3
"""Decks and verifier of the verification cases with an exact solution (ConvergenceOrder, MultigridSolution, VnnStability,
the cases of MPIBitIdentity). Standard library only (python 3.6), no ATLAS.

A unit box (2-D: one cell thick) with constant properties (k 52, rho 7850, cp 460), walls 302 (T = exact solution at the
face centre) on every outer face, mapped by a smooth, sheared, random or stretched lattice, optionally cut into px x py x pz
blocks joined by connections (101) and with 1, 2, 3 ... multigrid levels (bc<m>.txt written for every level: the wall value
is sampled at the face centre of ITS level). The solution is not linear, so the observed order is measurable:

  harm2  T = 1000 + 100 sin(pi x) sinh(pi y) / sinh(pi)                     harmonic (no source), 2-D
  pois2  T = 1000 + 100 sin(pi x) sin(pi y),  qvol = 2 pi^2 k 100 sin sin    Poisson, source in INPUT/st.tec, 2-D
  harm3  T = 1000 + 100 sin(pi x) sin(pi y) sinh(sqrt(2) pi z)/sinh(sqrt(2) pi)   harmonic, 3-D
  pois3  T = 1000 + 100 sin(pi x) sin(pi y) sin(pi z), qvol = 3 pi^2 k 100 sin sin sin, 3-D
  const  T = 1000 (the vnn cases: no steady gradient; --ic noise puts a random field on top)

Meshes (--mesh): cart; shear (2-D x = s + 0.5 t; 3-D x = s + 0.3 t + 0.1 u, y = t + 0.2 u, z = u + 0.25 s); rand (interior
nodes moved by 0.25 of the cell size in 2-D, 0.15 in 3-D, the same hash-seeded map at every resolution, scaled with h); stretch (the cells grow
1:5 geometrically in y (2-D) or z (3-D)); wavy (smooth sinusoidal map); jump (x cells 4 times wider for s > 0.5, sheared
x += 0.3 t: with --blocks 2,1,1 the jump sits on the connection); skew115 (rhombohedral lattice, edges 115 degrees apart in
all three pairs, as test/numerics-features/skewed_3D).

Commands
  make_case.py case  <dir> --sol S --mesh M --n N [options]      writes input.ini, MESH/, INPUT/, case.json
  make_case.py error <dir> [--no-walls]                          one line "ERR key=value ..." from OUTPUT/ of a finished run
  make_case.py diff  <dirA> <dirB>                               max |T_A - T_B| over the cells of the two OUTPUT/field.tec
  make_case.py order <name> <N:err ...>                          observed orders of consecutive resolutions
  make_case.py orders <name> <min L2> <min Linf> <N:L2:Linf ...> orders of the finest pair against the two minima (exit 1 below)
  make_case.py classify <dir> <exit status> <cap>                how a vnn run ended: converged, cap, diverged or failed
  make_case.py selftest                                          generator + verifier on a synthetic field (no solver)
Options of `case`: --blocks px,py,pz  --metric stencil|mean  --scheme euler|RK2|RK3  --vnn x  --vars cons|prim  --irs beta
  --levels L  --cycle none|v|w  --pre p --post q --coarse c  --cap n (iteration budget of the finest level)  --lviter n
  (budget of every coarser level; default: cap with cycle none, 1 with a cycle)  --res thr  --resditer n  --shellditer n
  --ic uniform|noise  --noise A  --k2 K2 (second material, i > N/2)  --probe  --decomp FILE  --dim 2|3.
"""
import argparse
import json
import math
import os
import re
import sys

K, RHO, CP, T0 = 52.0, 7850.0, 460.0, 1000.0
LZ = 0.01                               # thickness of the 2-D cases
PI = math.pi
SEED = 20261028
AMP = 0.25                              # random node displacement of a 2-D mesh, in cells
AMP3 = 0.15                             # ... of a 3-D mesh (0.25 in three directions folds a cell of the 32^3 lattice: cellcheck.py, 5 cells)
M64 = (1 << 64) - 1


# ----------------------------------------------------------------------------------------------------------- numbers
def e(v):
    return ' %.15E' % v


def iv(*v):
    return ''.join('%8d' % x for x in v)


def rv(*v):
    return ''.join('    %.15E,' % x for x in v)


def urand(seed, i, j, k, c):
    """Uniform in (-1, 1) from a hash of the node index: the same map in every block layout and at every level."""
    x = (seed * 0x9E3779B97F4A7C15 + i * 0xBF58476D1CE4E5B9 + j * 0x94D049BB133111EB + k * 0xD6E8FEB86659FD93
         + c * 0xCA5A826395121157 + 0x1234567) & M64
    x ^= x >> 30
    x = (x * 0xBF58476D1CE4E5B9) & M64
    x ^= x >> 27
    x = (x * 0x94D049BB133111EB) & M64
    x ^= x >> 31
    return 2.0 * (x / 18446744073709551616.0) - 1.0


# --------------------------------------------------------------------------------------------------- exact solutions
def solution(name):
    """(dim, T(x,y,z), grad T, qvol)."""
    if name == 'harm2':
        sh = math.sinh(PI)
        return (2, lambda x, y, z: T0 + 100.0 * math.sin(PI * x) * math.sinh(PI * y) / sh,
                lambda x, y, z: (100.0 * PI * math.cos(PI * x) * math.sinh(PI * y) / sh,
                                 100.0 * PI * math.sin(PI * x) * math.cosh(PI * y) / sh, 0.0),
                lambda x, y, z: 0.0)
    if name == 'pois2':
        return (2, lambda x, y, z: T0 + 100.0 * math.sin(PI * x) * math.sin(PI * y),
                lambda x, y, z: (100.0 * PI * math.cos(PI * x) * math.sin(PI * y),
                                 100.0 * PI * math.sin(PI * x) * math.cos(PI * y), 0.0),
                lambda x, y, z: 2.0 * PI * PI * K * 100.0 * math.sin(PI * x) * math.sin(PI * y))
    if name == 'harm3':
        c = math.sqrt(2.0) * PI
        sh = math.sinh(c)
        return (3, lambda x, y, z: T0 + 100.0 * math.sin(PI * x) * math.sin(PI * y) * math.sinh(c * z) / sh,
                lambda x, y, z: (100.0 * PI * math.cos(PI * x) * math.sin(PI * y) * math.sinh(c * z) / sh,
                                 100.0 * PI * math.sin(PI * x) * math.cos(PI * y) * math.sinh(c * z) / sh,
                                 100.0 * c * math.sin(PI * x) * math.sin(PI * y) * math.cosh(c * z) / sh),
                lambda x, y, z: 0.0)
    if name == 'pois3':
        return (3, lambda x, y, z: T0 + 100.0 * math.sin(PI * x) * math.sin(PI * y) * math.sin(PI * z),
                lambda x, y, z: (100.0 * PI * math.cos(PI * x) * math.sin(PI * y) * math.sin(PI * z),
                                 100.0 * PI * math.sin(PI * x) * math.cos(PI * y) * math.sin(PI * z),
                                 100.0 * PI * math.sin(PI * x) * math.sin(PI * y) * math.cos(PI * z)),
                lambda x, y, z: 3.0 * PI * PI * K * 100.0 * math.sin(PI * x) * math.sin(PI * y) * math.sin(PI * z))
    if name == 'const':
        return (0, lambda x, y, z: T0, lambda x, y, z: (0.0, 0.0, 0.0), lambda x, y, z: 0.0)
    sys.exit('unknown solution ' + name)


# ----------------------------------------------------------------------------------------------------------- meshes
def skew115_edges():
    c = math.cos(math.radians(115.0))
    s2 = math.sqrt(1 - c * c)
    y3 = (c - c * c) / s2
    return [(1.0, 0.0, 0.0), (c, s2, 0.0), (c, y3, math.sqrt(1 - c * c - y3 * y3))]


def base_map(mesh, dim, s, t, u):
    if mesh in ('cart', 'rand'):
        return s, t, u
    if mesh == 'shear':
        if dim == 2:
            return s + 0.5 * t, t, u
        return s + 0.3 * t + 0.1 * u, t + 0.2 * u, u + 0.25 * s
    if mesh == 'stretch':
        a = math.log(5.0)
        f = lambda w: (math.exp(a * w) - 1.0) / (math.exp(a) - 1.0)
        return (s, f(t), u) if dim == 2 else (s, t, f(u))
    if mesh == 'wavy':
        if dim == 2:
            return s + 0.06 * math.sin(2 * PI * t), t + 0.06 * math.sin(2 * PI * s), u
        return (s + 0.04 * math.sin(2 * PI * t), t + 0.04 * math.sin(2 * PI * u), u + 0.04 * math.sin(2 * PI * s))
    if mesh == 'jump':
        x = 0.4 * s if s <= 0.5 else 0.2 + 1.6 * (s - 0.5)
        return x + 0.3 * t, t, u
    sys.exit('unknown mesh ' + mesh)


def lattice_node(mesh, dim, N, i, j, k):
    NX, NY, NZ = N
    if mesh == 'skew115':
        ed = skew115_edges()
        return tuple((i * ed[0][d] / NX + j * ed[1][d] / NY + k * ed[2][d] / NZ) for d in range(3))
    s, t = i / float(NX), j / float(NY)
    u = k / float(NZ) if dim == 3 else 0.0
    x, y, z = base_map(mesh, dim, s, t, u)
    if dim == 2:
        z = (k - 0.5) * LZ
    if mesh == 'rand' and 0 < i < NX and 0 < j < NY and (dim == 2 or 0 < k < NZ):
        kk = 0 if dim == 2 else k
        amp = AMP if dim == 2 else AMP3
        x += amp / NX * urand(SEED, i, j, kk, 0)
        y += amp / NY * urand(SEED, i, j, kk, 1)
        if dim == 3:
            z += amp / NZ * urand(SEED, i, j, kk, 2)
    return x, y, z


class Blk(object):
    """A structured block: nx x ny x nz cells, nodes in the order of a Tecplot BLOCK zone (i fastest)."""

    def __init__(self, nx, ny, nz, X, Y, Z):
        self.nx, self.ny, self.nz, self.X, self.Y, self.Z = nx, ny, nz, X, Y, Z

    def nid(self, i, j, k):
        return (k * (self.ny + 1) + j) * (self.nx + 1) + i

    def pt(self, i, j, k):
        n = self.nid(i, j, k)
        return (self.X[n], self.Y[n], self.Z[n])

    def face_cells(self, f):
        nx, ny, nz = self.nx, self.ny, self.nz
        if f in (1, 2):
            return [(1 if f == 1 else nx, j, k) for k in range(1, nz + 1) for j in range(1, ny + 1)]
        if f in (3, 4):
            return [(i, 1 if f == 3 else ny, k) for k in range(1, nz + 1) for i in range(1, nx + 1)]
        return [(i, j, 1 if f == 5 else nz) for j in range(1, ny + 1) for i in range(1, nx + 1)]

    def corners(self, i, j, k):
        """The eight corners of cell (i, j, k) (1-based) as P[c][b][a]."""
        return [[[self.pt(i - 1 + a, j - 1 + b, k - 1 + c) for a in (0, 1)] for b in (0, 1)] for c in (0, 1)]

    def face_pts(self, f, i, j, k):
        i0, j0, k0 = i - 1, j - 1, k - 1
        q = {1: [(i0, j0, k0), (i0, j0 + 1, k0), (i0, j0 + 1, k0 + 1), (i0, j0, k0 + 1)],
             2: [(i0 + 1, j0, k0), (i0 + 1, j0 + 1, k0), (i0 + 1, j0 + 1, k0 + 1), (i0 + 1, j0, k0 + 1)],
             3: [(i0, j0, k0), (i0 + 1, j0, k0), (i0 + 1, j0, k0 + 1), (i0, j0, k0 + 1)],
             4: [(i0, j0 + 1, k0), (i0 + 1, j0 + 1, k0), (i0 + 1, j0 + 1, k0 + 1), (i0, j0 + 1, k0 + 1)],
             5: [(i0, j0, k0), (i0 + 1, j0, k0), (i0 + 1, j0 + 1, k0), (i0, j0 + 1, k0)],
             6: [(i0, j0, k0 + 1), (i0 + 1, j0, k0 + 1), (i0 + 1, j0 + 1, k0 + 1), (i0, j0 + 1, k0 + 1)]}[f]
        return [self.pt(*c) for c in q]

    def face_centre(self, f, i, j, k):
        p = self.face_pts(f, i, j, k)
        return tuple(sum(q[d] for q in p) / 4.0 for d in range(3))


def cell_centre_of(P):
    pts = [P[c][b][a] for c in (0, 1) for b in (0, 1) for a in (0, 1)]
    return tuple(sum(q[d] for q in pts) / 8.0 for d in range(3))


def vec_area(p):
    """Area vector (1/2 d1 x d2) of the quadrilateral p[0..3] (in order around it)."""
    d1 = [p[2][d] - p[0][d] for d in range(3)]
    d2 = [p[3][d] - p[1][d] for d in range(3)]
    return [0.5 * (d1[1] * d2[2] - d1[2] * d2[1]), 0.5 * (d1[2] * d2[0] - d1[0] * d2[2]), 0.5 * (d1[0] * d2[1] - d1[1] * d2[0])]


def hex_volume(blk, i, j, k):
    """Volume of cell (i, j, k): (1/3) sum over the six faces of (area vector . face centre), oriented outward."""
    P = blk.corners(i, j, k)
    cc = cell_centre_of(P)
    v = 0.0
    for f in range(1, 7):
        p = blk.face_pts(f, i, j, k)
        S = vec_area(p)
        c = [sum(q[d] for q in p) / 4.0 for d in range(3)]
        if sum(S[d] * (c[d] - cc[d]) for d in range(3)) < 0:
            S = [-s for s in S]
        v += sum(S[d] * c[d] for d in range(3))
    return v / 3.0


def build_blocks(mesh, dim, n, layout):
    """Level-1 blocks of the lattice n = (NX, NY, NZ) cut into layout = (px, py, pz) blocks (block number a + px (b + py c) + 1)."""
    NX, NY, NZ = n
    px, py, pz = layout
    if NX % px or NY % py or NZ % pz:
        sys.exit('the cells %s are not divisible by the blocks %s' % (n, layout))
    mx, my, mz = NX // px, NY // py, NZ // pz
    blocks = []
    for c in range(pz):
        for b in range(py):
            for a in range(px):
                X, Y, Z = [], [], []
                for k in range(mz + 1):
                    for j in range(my + 1):
                        for i in range(mx + 1):
                            x, y, z = lattice_node(mesh, dim, n, a * mx + i, b * my + j, c * mz + k)
                            X.append(x); Y.append(y); Z.append(z)
                blocks.append(Blk(mx, my, mz, X, Y, Z))
    return blocks


def coarsen(blocks, m):
    """Level m (1 = finest): every 2^(m-1)-th node in i, j and (3-D) k, as Coarse_Grid does."""
    s = 2 ** (m - 1)
    out = []
    for b in blocks:
        sk = s if b.nz > 1 else 1
        nx, ny, nz = b.nx // s, b.ny // s, b.nz // sk
        X, Y, Z = [], [], []
        for k in range(nz + 1):
            for j in range(ny + 1):
                for i in range(nx + 1):
                    n = b.nid(i * s, j * s, k * sk)
                    X.append(b.X[n]); Y.append(b.Y[n]); Z.append(b.Z[n])
        out.append(Blk(nx, ny, nz, X, Y, Z))
    return out


# ------------------------------------------------------------------------------------------------------- the files
def zone_header(title, b, varloc):
    return (' ZONE  T = %s, I=%d, J=%d, K=%d, DATAPACKING=BLOCK, VARLOCATION=(%s), SOLUTIONTIME=0.0, STRANDID = 0\n'
            % (title, b.nx + 1, b.ny + 1, b.nz + 1, varloc))


def nodes_flat(b):
    return b.X + b.Y + b.Z


def cells_of(b):
    return [(i, j, k) for k in range(1, b.nz + 1) for j in range(1, b.ny + 1) for i in range(1, b.nx + 1)]


def bc_lines(blocks, layout, dim, Tfun):
    px, py, pz = layout
    lines = []
    for bn, b in enumerate(blocks, start=1):
        a = (bn - 1) % px
        bb = ((bn - 1) // px) % py
        c = (bn - 1) // (px * py)
        for f in range(1, 7):
            nbr = None
            if f == 1 and a > 0:
                nbr = (bn - 1, 2)
            elif f == 2 and a < px - 1:
                nbr = (bn + 1, 1)
            elif f == 3 and bb > 0:
                nbr = (bn - px, 4)
            elif f == 4 and bb < py - 1:
                nbr = (bn + px, 3)
            elif f == 5 and c > 0:
                nbr = (bn - px * py, 6)
            elif f == 6 and c < pz - 1:
                nbr = (bn + px * py, 5)
            for (i, j, k) in b.face_cells(f):
                if dim == 2 and f in (5, 6):
                    lines.append(iv(bn, i, j, k, f, 0))
                elif nbr:
                    ob, of = nbr
                    o = blocks[ob - 1]
                    oi, oj, ok = i, j, k
                    if of == 1:
                        oi = 1
                    elif of == 2:
                        oi = o.nx
                    elif of == 3:
                        oj = 1
                    elif of == 4:
                        oj = o.ny
                    elif of == 5:
                        ok = 1
                    else:
                        ok = o.nz
                    lines += [iv(bn, i, j, k, f, 101), iv(ob, oi, oj, ok, of, 1, 0, 0, 1)]
                else:
                    lines += [iv(bn, i, j, k, f, 302), rv(Tfun(*b.face_centre(f, i, j, k)))]
    return lines


def ini_text(a, dim):
    L = a.levels
    cycle = a.cycle if L > 1 else 'none'
    cap = a.cap
    lines = ['[FUSS-Parameters]', 'newrun = true', 'res-threshold = %s' % a.res]
    if L == 1:
        lines.append('iter-threshold = %d' % cap)
    lines += ['', '[FUSS-Numerics]', 'time-scheme = %s' % a.scheme, 'vnn = %s' % a.vnn, 'time-accurate = false',
              'integration-variables = %s' % a.vars]
    if a.metric != 'stencil':
        lines.append('diffusive-metric = %s' % a.metric)
    if a.irs > 0:
        lines.append('irs-beta = %s' % a.irs)
    if L > 1:
        lv = a.lviter if a.lviter is not None else (cap if cycle == 'none' else 1)
        lines += ['', '[FUSS-Multigrid]', 'levels = %d' % L, 'cycle = %s' % cycle, 'level1-iter = %d' % cap]
        lines += ['level%d-iter = %d' % (m, lv) for m in range(2, L + 1)]
        if cycle != 'none':
            lines += ['pre-iter = %d' % a.pre, 'post-iter = %d' % a.post, 'coarse-iter = %d' % a.coarse]
    shell = a.shellditer if a.shellditer else (100 if cycle == 'none' else 1)
    lines += ['', '[FUSS-IO]', 'sol-diter = 1000000000', 'sol-overwrite = true', 'sol-format = tecplot ascii',
              'res-diter = %d' % a.resditer, 'shell-diter = %d' % shell, 'ini-diter = 1000000000']
    if a.decomp:
        lines.append('decomposition-map = %s' % a.decomp)
    if a.probe:
        lines += ['', '[FUSS-Probes]', 'probe1 = p1', '', '[p1]', 'index-position = 1 2 2 1', 'variables = T',
                  'diter = %d' % a.probediter]
    return '\n'.join(lines) + '\n'


def write_case(a):
    solname = a.sol
    sdim, Tex, Gex, Qex = solution(solname)
    dim = a.dim or sdim
    if dim not in (2, 3):
        sys.exit('--dim 2|3 needed')
    n = (a.n, a.n, a.n if dim == 3 else 1)
    layout = tuple(int(v) for v in a.blocks.split(','))
    L = a.levels
    out = a.dir
    for sub in ('MESH', 'INPUT', 'OUTPUT'):
        os.makedirs(os.path.join(out, sub), exist_ok=True)
    fine = build_blocks(a.mesh, dim, n, layout)
    rap = 2 ** (L - 1)
    for b in fine:
        for d, nn in enumerate((b.nx, b.ny, b.nz)):
            if nn > 1 and (nn % rap or nn // rap < 2):
                sys.exit('block of %dx%dx%d cells: %d levels need every direction divisible by %d with 2 cells on the coarsest level'
                         % (b.nx, b.ny, b.nz, L, rap))
    mesh_txt, ic_txt, st_txt = [' VARIABLES ="x" "y" "z"\n'], [' VARIABLES ="x" "y" "z" T matID\n'], [' VARIABLES ="x" "y" "z"  qvol\n']
    anysrc = False
    for bn, b in enumerate(fine, start=1):
        flat = nodes_flat(b)
        mesh_txt += [zone_header('Block%d' % bn, b, '[1-3]=NODAL'), '\n'.join(e(v) for v in flat) + '\n']
        cells = cells_of(b)
        T, mat, q = [], [], []
        for (i, j, k) in cells:
            P = b.corners(i, j, k)
            cx = cell_centre_of(P)
            if a.ic == 'noise':
                u = urand(SEED + 7, (bn << 20) + i, j, k, 3)
                T.append(T0 + a.noise * u)
            elif a.ic == 'exact':
                T.append(Tex(*cx))
            else:
                T.append(T0)
            mat.append(2.0 if (a.k2 and (i - 1 + (bn - 1) % layout[0] * b.nx) >= n[0] // 2) else 1.0)
            qq = Qex(*cx)
            q.append(qq)
            anysrc = anysrc or qq != 0.0
        ic_txt += [zone_header('B%d-SP' % bn, b, '[1-3]=NODAL,[4-5]=CELLCENTERED'), '\n'.join(e(v) for v in flat + T + mat) + '\n']
        st_txt += [zone_header('B%d-ST' % bn, b, '[1-3]=NODAL,[4]=CELLCENTERED'), '\n'.join(e(v) for v in flat + q) + '\n']
    open(os.path.join(out, 'MESH', 'mesh.tec'), 'w').write(''.join(mesh_txt))
    open(os.path.join(out, 'INPUT', 'ic.tec'), 'w').write(''.join(ic_txt))
    if anysrc:
        open(os.path.join(out, 'INPUT', 'st.tec'), 'w').write(''.join(st_txt))
    for m in range(1, L + 1):
        blocks = fine if m == 1 else coarsen(fine, m)
        name = 'bc.txt' if m == 1 else 'bc%d.txt' % m
        open(os.path.join(out, 'INPUT', name), 'w').write('\n'.join(bc_lines(blocks, layout, dim, Tex)) + '\n')
    mats = [('mat1', K), ('mat2', a.k2)] if a.k2 else [('mat', K)]
    open(os.path.join(out, 'INPUT', 'phase.txt'), 'w').write('solid-bulk\n' + ''.join('%s 1\n' % nm for nm, _ in mats))
    with open(os.path.join(out, 'INPUT', 'properties.dat'), 'w') as fh:
        fh.write('TITLE = "Mass Thermodynamic Properties"\nVARIABLES = "Temperature", "Cp", "Density", "Conductivity", "Energy"\n')
        for nm, kk in mats:
            fh.write('ZONE T="%s"\nI=2500, F=POINT\n' % nm)
            fh.writelines('%.1f %.6f %.6f %.6f %.6f\n' % (t, CP, RHO, kk, RHO * CP * t) for t in range(1, 2501))
    open(os.path.join(out, 'input.ini'), 'w').write(ini_text(a, dim))
    meta = {'sol': solname, 'mesh': a.mesh, 'N': a.n, 'dim': dim, 'layout': list(layout), 'metric': a.metric,
            'levels': L, 'cycle': a.cycle if L > 1 else 'none', 'vars': a.vars, 'irs': a.irs, 'k2': a.k2}
    json.dump(meta, open(os.path.join(out, 'case.json'), 'w'))
    print('%s %s N=%d blocks=%s levels=%d %s: %d cells' % (solname, a.mesh, a.n, a.blocks, L, a.cycle if L > 1 else 'single',
                                                          sum(b.nx * b.ny * b.nz for b in fine)))


# ------------------------------------------------------------------------------------------------------- verifier
def read_tec(path):
    """Tecplot ASCII BLOCK file -> list of (title, dims, {variable: values}) (as fuss-stress/tools/fusstec.read_tec)."""
    text = open(path).read().split('\n')
    names, zones, i = [], [], 0
    while i < len(text):
        line = text[i].strip()
        if re.match(r'(?i)variables', line):
            names = [x or y for x, y in re.findall(r'"([^"]*)"|([^\s,"]+)', line.split('=', 1)[1])]
        elif re.match(r'(?i)zone', line):
            title = re.search(r'T\s*=\s*"?([^",\s]+)', line).group(1)
            ni, nj, nk = (int(re.search(r'\b%s\s*=\s*(\d+)' % c, line).group(1)) for c in 'IJK')
            cc = re.search(r'\[(\d+)-(\d+)\]\s*=\s*CELLCENTERED', line, re.I)
            sizes = [(max(ni - 1, 1) * max(nj - 1, 1) * max(nk - 1, 1))
                     if cc and int(cc.group(1)) <= v <= int(cc.group(2)) else ni * nj * nk for v in range(1, len(names) + 1)]
            vals, i = [], i + 1
            while len(vals) < sum(sizes) and i < len(text):
                s = text[i].strip()
                if s and not s[0].isalpha():
                    vals += [float(x) for x in s.split()]
                elif re.match(r'(?i)(zone|dataset|text)', s):
                    break
                i += 1
            zone, off = {}, 0
            for nm, sz in zip(names, sizes):
                zone[nm] = vals[off:off + sz]
                off += sz
            zones.append((title, (ni, nj, nk), zone))
            continue
        i += 1
    return zones


def zone_block(dims, z):
    ni, nj, nk = dims
    return Blk(ni - 1, nj - 1, nk - 1, z['x'], z['y'], z['z'])


def error_line(rundir, with_walls=True):
    meta = json.load(open(os.path.join(rundir, 'case.json')))
    dim, Tex, Gex, Qex = solution(meta['sol'])
    zones = read_tec(os.path.join(rundir, 'OUTPUT', 'field.tec'))
    sv = se = 0.0
    linf, at = 0.0, None
    ncell = 0
    for bn, (title, dims, z) in enumerate(zones, start=1):
        b = zone_block(dims, z)
        T = z['T']
        for idx, (i, j, k) in enumerate(cells_of(b)):
            cx = cell_centre_of(b.corners(i, j, k))
            err = T[idx] - Tex(*cx)
            v = hex_volume(b, i, j, k)
            sv += v
            se += v * err * err
            if abs(err) > linf:
                linf, at = abs(err), (bn, i, j, k)
            ncell += 1
    out = {'cells': ncell, 'l2': math.sqrt(se / sv), 'linf': linf, 'at': '%d:%d,%d,%d' % at if at else '-'}
    wpath = os.path.join(rundir, 'OUTPUT', 'wall.tec')
    if with_walls and os.path.exists(wpath):
        blocks = [zone_block(d, z) for (_, d, z) in zones]
        sa = sw = sai = swi = sq = 0.0
        wlinf = wlinfi = 0.0
        wat = '-'
        tang = {1: (1, 2), 2: (1, 2), 3: (0, 2), 4: (0, 2), 5: (0, 1), 6: (0, 1)}
        for title, dims, z in read_tec(wpath):
            m = re.match(r'B(\d+)F(\d+)', title)
            if not m or 'qw' not in z:
                continue
            bn, f = int(m.group(1)), int(m.group(2))
            b = blocks[bn - 1]
            nn = (b.nx, b.ny, b.nz)
            for idx, (i, j, k) in enumerate(b.face_cells(f)):
                p = b.face_pts(f, i, j, k)
                S = vec_area(p)
                c = [sum(q[d] for q in p) / 4.0 for d in range(3)]
                cc = cell_centre_of(b.corners(i, j, k))
                if sum(S[d] * (c[d] - cc[d]) for d in range(3)) < 0:
                    S = [-s for s in S]
                A = math.sqrt(sum(s * s for s in S))
                g = Gex(*c)
                qex = K * sum(g[d] * S[d] for d in range(3)) / A
                d = z['qw'][idx] - qex
                ijk = (i, j, k)
                edge = any(nn[t] > 1 and ijk[t] in (1, nn[t]) for t in tang[f])
                sa += A
                sw += A * d * d
                sq += A * qex * qex
                if abs(d) > wlinf:
                    wlinf, wat = abs(d), '%d:F%d:%d,%d,%d%s' % (bn, f, i, j, k, ':edge' if edge else ':inner')
                if not edge:
                    sai += A
                    swi += A * d * d
                    wlinfi = max(wlinfi, abs(d))
        if sa > 0:
            out.update({'wl2': math.sqrt(sw / sa), 'wlinf': wlinf, 'qrms': math.sqrt(sq / sa), 'wat': wat})
            if sai > 0:
                out.update({'wl2i': math.sqrt(swi / sai), 'wlinfi': wlinfi})
    return meta, out


def print_error(rundir, with_walls=True):
    meta, o = error_line(rundir, with_walls)
    print('ERR ' + ' '.join('%s=%s' % (k, ('%.6e' % v) if isinstance(v, float) else v) for k, v in sorted(o.items())))
    return 0


def global_cells(rundir):
    """{(I, J, K): T} of a run with the cells numbered over the whole lattice (the blocks of case.json's layout put side by side)."""
    meta = json.load(open(os.path.join(rundir, 'case.json')))
    px, py, pz = meta['layout']
    out = {}
    for bn, (title, dims, z) in enumerate(read_tec(os.path.join(rundir, 'OUTPUT', 'field.tec')), start=1):
        a, b, c = (bn - 1) % px, ((bn - 1) // px) % py, (bn - 1) // (px * py)
        nx, ny, nz = dims[0] - 1, dims[1] - 1, dims[2] - 1
        T = z['T']
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    out[(a * nx + i, b * ny + j, c * nz + k)] = T[(k * ny + j) * nx + i]
    return out


def field_diff(da, db):
    """max |T_A - T_B| over the cells of two runs of the same lattice, whatever their block layouts."""
    A, B = global_cells(da), global_cells(db)
    if set(A) != set(B):
        return float('nan')
    return max(abs(A[k] - B[k]) for k in A)


def order(name, pairs):
    rows = []
    for p in pairs:
        n, v = p.split(':')
        rows.append((int(n), float(v)))
    rows.sort()
    out = []
    for (n1, v1), (n2, v2) in zip(rows, rows[1:]):
        out.append('%d->%d: %.3f' % (n1, n2, math.log(v1 / v2) / math.log(n2 / float(n1))) if v1 > 0 and v2 > 0 else '%d->%d: -' % (n1, n2))
    print('%s  %s' % (name, '  '.join(out)))


def classify(rundir, ex, cap):
    """How a run of the const solution ended, from its exit status, residual history, log and field. One line
    "RUN status=S iterations=N last-residual=R max-dev=D finite=F msg=...", S being
      cap        exit 0 and the budget `cap` used up: the residual did not reach the threshold (a stable, not converging field)
      converged  exit 0 before the budget: the residual reached the threshold (the field ended at its fixed point)
      diverged   exit != 0 with "Integration failed" in the log (the run stopped itself on a non-physical state)
      failed     any other exit != 0, or no history
    max-dev = max |T - T0| over the cells of OUTPUT/field.tec (empty when the run wrote none), finite = all T finite."""
    it, last = 0, float('nan')
    hp = os.path.join(rundir, 'OUTPUT', 'residual-history.dat')
    if os.path.exists(hp):
        for line in open(hp, errors='replace'):
            p = line.split()
            if len(p) >= 3:
                try:
                    it, last = int(p[0]), float(p[2])
                except ValueError:
                    pass
    msg, integ = '-', False
    lp = os.path.join(rundir, 'run.log')
    if os.path.exists(lp):
        for line in open(lp, errors='replace'):
            if 'Integration failed' in line:
                integ = True
            if msg == '-' and re.search(r'\[ERROR\]|STOP:|forrtl|Fatal', line):
                msg = line.strip()[:120].replace(' ', '_')
    dev, finite = '', ''
    fld = os.path.join(rundir, 'OUTPUT', 'field.tec')
    if os.path.exists(fld):
        t = [v for (_, _, z) in read_tec(fld) for v in z['T']]
        finite = str(all(math.isfinite(v) for v in t))
        dev = '%.6e' % max(abs(v - T0) for v in t if math.isfinite(v)) if t else ''
    if int(ex) != 0:
        status = 'diverged' if integ else 'failed'
    elif it >= cap:
        status = 'cap'
    elif it > 0:
        status = 'converged'
    else:
        status = 'failed'
    print('RUN status=%s iterations=%d last-residual=%.4e max-dev=%s finite=%s msg=%s' % (status, it, last, dev or '-', finite or '-', msg))
    return 0


def orders(name, min_l2, min_linf, triples):
    """Observed orders of consecutive resolutions in L2 and Linf (triples N:l2:linf); the finest pair must reach the minima."""
    rows = []
    for t in triples:
        n, a, b = t.split(':')
        rows.append((int(n), float(a), float(b)))
    rows.sort()
    pl2, pli = [], []
    for (n1, a1, b1), (n2, a2, b2) in zip(rows, rows[1:]):
        r = math.log(n2 / float(n1))
        pl2.append(math.log(a1 / a2) / r)
        pli.append(math.log(b1 / b2) / r)
    ok = bool(pl2) and pl2[-1] >= min_l2 and pli[-1] >= min_linf
    print('%-22s L2 %s | Linf %s | finest pair %d->%d: p(L2) %.2f (>= %.1f), p(Linf) %.2f (>= %.1f)  %s'
          % (name, ' '.join('%.2f' % p for p in pl2), ' '.join('%.2f' % p for p in pli), rows[-2][0], rows[-1][0],
             pl2[-1], min_l2, pli[-1], min_linf, 'ok' if ok else 'FAIL'))
    return 0 if ok else 1


# ----------------------------------------------------------------------------------------------------------- self-test
def selftest():
    import shutil
    import tempfile
    tmp = tempfile.mkdtemp()
    ok = True
    for sol, mesh, n, lay in (('harm2', 'rand', 16, '2,2,1'), ('harm3', 'shear', 16, '2,1,1'), ('pois2', 'jump', 16, '2,1,1')):
        d = os.path.join(tmp, '%s-%s' % (sol, mesh))
        args = parse(['case', d, '--sol', sol, '--mesh', mesh, '--n', str(n), '--blocks', lay, '--levels', '3', '--cycle', 'w'])
        write_case(args)
        # fabricate OUTPUT/field.tec = exact + 1e-3 (constant offset): the error must be 1e-3 in both norms
        sdim, Tex, Gex, Qex = solution(sol)
        zones = []
        b = None
        mt = read_tec(os.path.join(d, 'MESH', 'mesh.tec'))
        txt = [' VARIABLES ="x" "y" "z" "T"\n']
        for title, dims, z in mt:
            b = zone_block(dims, z)
            T = [Tex(*cell_centre_of(b.corners(i, j, k))) + 1e-3 for (i, j, k) in cells_of(b)]
            txt += [' ZONE  T = %s, I=%d, J=%d, K=%d, DATAPACKING=BLOCK, VARLOCATION=([1-3]=NODAL,[4]=CELLCENTERED)\n'
                    % (title, dims[0], dims[1], dims[2]), '\n'.join(e(v) for v in z['x'] + z['y'] + z['z'] + T) + '\n']
        open(os.path.join(d, 'OUTPUT', 'field.tec'), 'w').write(''.join(txt))
        meta, o = error_line(d)
        good = abs(o['l2'] - 1e-3) < 1e-9 and abs(o['linf'] - 1e-3) < 1e-9 and o['cells'] > 0
        print('selftest %-12s cells %d l2 %.6e linf %.6e %s' % (sol + '-' + mesh, o['cells'], o['l2'], o['linf'], 'ok' if good else 'FAIL'))
        ok = ok and good
        for f in sorted(os.listdir(os.path.join(d, 'INPUT'))):
            print('   INPUT/%s %d bytes' % (f, os.path.getsize(os.path.join(d, 'INPUT', f))))
    print(open(os.path.join(tmp, 'harm2-rand', 'input.ini')).read())
    shutil.rmtree(tmp)
    return 0 if ok else 1


def parse(argv):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('cmd', choices=['case'])
    p.add_argument('dir')
    p.add_argument('--sol', default='harm2')
    p.add_argument('--mesh', default='cart')
    p.add_argument('--n', type=int, default=16)
    p.add_argument('--dim', type=int, default=0)
    p.add_argument('--blocks', default='1,1,1')
    p.add_argument('--metric', default='stencil', choices=['stencil', 'mean'])
    p.add_argument('--scheme', default='RK3')
    p.add_argument('--vnn', default='0.4')
    p.add_argument('--vars', default='cons', choices=['cons', 'prim'])
    p.add_argument('--irs', type=float, default=0.0)
    p.add_argument('--levels', type=int, default=1)
    p.add_argument('--cycle', default='none', choices=['none', 'v', 'w'])
    p.add_argument('--pre', type=int, default=2)
    p.add_argument('--post', type=int, default=1)
    p.add_argument('--coarse', type=int, default=10)
    p.add_argument('--cap', type=int, default=0)
    p.add_argument('--lviter', type=int, default=None)
    p.add_argument('--res', default='1e-10')
    p.add_argument('--resditer', type=int, default=1)
    p.add_argument('--shellditer', type=int, default=0)
    p.add_argument('--ic', default='uniform', choices=['uniform', 'noise', 'exact'])
    p.add_argument('--noise', type=float, default=50.0)
    p.add_argument('--k2', type=float, default=0.0)
    p.add_argument('--probe', action='store_true')
    p.add_argument('--probediter', type=int, default=10)
    p.add_argument('--decomp', default='')
    a = p.parse_args(argv)
    if not a.cap:
        a.cap = 100 * a.n * a.n
    return a


if __name__ == '__main__':
    if len(sys.argv) >= 3 and sys.argv[1] == 'error':
        sys.exit(print_error(sys.argv[2], '--no-walls' not in sys.argv[3:]))
    if len(sys.argv) >= 4 and sys.argv[1] == 'diff':
        print('%.6e' % field_diff(sys.argv[2], sys.argv[3]))
        sys.exit(0)
    if len(sys.argv) >= 3 and sys.argv[1] == 'order':
        order(sys.argv[2], sys.argv[3:])
        sys.exit(0)
    if len(sys.argv) >= 5 and sys.argv[1] == 'classify':
        sys.exit(classify(sys.argv[2], sys.argv[3], int(sys.argv[4])))
    if len(sys.argv) >= 6 and sys.argv[1] == 'orders':
        sys.exit(orders(sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), sys.argv[5:]))
    if len(sys.argv) >= 2 and sys.argv[1] == 'selftest':
        sys.exit(selftest())
    if len(sys.argv) >= 2 and sys.argv[1] == 'case':
        write_case(parse(sys.argv[1:]))
        sys.exit(0)
    sys.exit(__doc__)
