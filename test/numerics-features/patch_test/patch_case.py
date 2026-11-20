#!/usr/bin/env python3
"""Linear-field patch test of the FUSS diffusive fluxes (the conduction analogue of MOSE's test/fast/viscous-patch).

A linear field T = A + g.x has a uniform heat flux, so the net conduction of every cell is zero: a face gradient that is
exact for linear fields gives a zero residual to round-off on any mesh. One forward-Euler step from the exact field, with
the wall values of the exact field (302: T at the face centre; 301: q = k g.n_out into the solid; 303: Tref = Tw + q/h;
305: Tref from h (Tref - Tw) + sigma eps (Tref^4 - Tw^4) = q), and OUTPUT/diagnostic.tec holds the residual E of every cell.
The check reports max |E| / (k |g| A_cell), A_cell the total face area of the cell, for the interior cells, the cells next to
block connections (101), to symmetry faces (300) and to each wall type, and max |qw - exact| / (k |g|) of the wall output.

Families (2-D cases are one cell thick, k faces of type 0; nodes moved at random with a fixed seed where said):
  para     one block 16x16 of parallelogram cells (shear 20 deg), walls 302: the control, every metric is exact here
  rand     one block 16x16, interior nodes moved by up to 0.25 of the cell size, walls 302
  jump     one block 16x16, cell sizes alternating in the ratio 1.5 in x and y, sheared by x += 0.5 y, walls 302
  fold2    two blocks 8x8 joined at x = 0.5 by a 101 connection, the j lines of block 2 rising as y = t + 0.25 s t
           (a fold at the interface), nodes moved by 0.25, walls 302
  sym      one block 16x16 sheared by x += 0.35 y on a symmetry plane (300) at y = 0, field T = A + B x (no y gradient),
           nodes moved by 0.15, other walls 302
  wallmix  as rand with walls 302 (face 1), 303 (face 2), 301 (face 3), 305 (face 4): the wall formulation on random cells
  box3d    two blocks 6x6x6 joined at x = 0.5, the j lines of block 2 rising as y = t + 0.2 s t, nodes moved by 0.15,
           walls 302 on the six outer faces (3-D: the k faces are walls too)
  step     one block 16x16, Cartesian: 8 cells of width 0.05 then 8 of 0.075 in x (one jump of the cell size, ratio 1.5),
           y uniform, walls 302. Exact under the stencil metric; under diffusive-metric = mean the cells at the jump
           plane carry 6.4201e-3 (left) and 5.2528e-3 (right) of k|g|A_cell (the (1+r)^2/(4r) - 1 bias of that metric,
           invisible to the alternating sizes of `jump`)

Families added with the verification of October 2026 (ported from MOSE test/fast/viscous-patch/make_case.py and extended; the
battery of the sync, fuss-stress/tools/patch_case.py, has the same ones and a few more):
  jump2    as jump with the sizes 1 and 2 (the first cells next to faces 1 and 3 half as wide), dyadic coordinates (unit 2^-7 m),
           walls 302/301/302/303 on faces 1..4 (MOSE: jump2)
  fan      8 x 8 cells, the j spacing grows 7/4 across the first cell at face 1 (1/16 at x = 0, 7/64 after), dyadic coordinates,
           walls 302/301/302/303 (MOSE: fan)
  fan4     as fan, spacing 1/4 after the first cell: both ghost layers behind face 1 turned over (MOSE: fan4)
  star3    a triangle cut into 3 blocks of 8 x 8 cells meeting at a point shared by 3 blocks; connections 1 <-> 3 (same parity,
  star5    index axes exchanged); a pentagon in 5 blocks; outer walls by block (302 301 303 305), cells at the shared point
           form the group `centre`; nodes moved by 0.15 of the cell size (MOSE: star3, star5)
  rough3d  3 blocks 6x6x6 along x, z stretched 1:5, folds at both interfaces (slopes +0.2, -0.3), nodes moved by 0.25, walls 302
  wall3dA  one skewed 6x6x6 block (shear in the three planes, nodes moved by 0.15), faces 1..6 = 302 301 303 305 302 301
  wall3dB  the same with the types shifted: 301 303 305 302 303 305 (each type on at least two different faces over A and B)
  sym3d    skewed block with symmetry plane (300) at z = 0 (x, y sheared with z, the plane flat), other faces 302 301 303 305 302
  conn14   2 blocks 8x8 joined by a 101, cell width ratio 1:4 across it, fold (block 2 j lines rise) and shear 0.3
  conn14s  the same with the block order swapped (the large cells are block 1, the connection is 1 <-> 2 seen from the other side)
  conn14p  control: conn14 without fold or moved nodes (same metric in every cell of a block)
  flipA    2 blocks, block 2 with both index axes reversed: face 2 <-> face 2 with j flipped (d = -1 0 0 1 form), nodes moved
  flipB    2 blocks, block 2 with the axes exchanged: face 2 <-> face 3, i reversed on the partner, nodes moved
  conn14c  conn14 in `conserve` mode: random T field (1000 +- 100 K), adiabatic walls (301, q = 0): the sum of the cell
           residuals is 0, so no flux is seen differently by the two sides of a connection (`check` prints CONSERVE)
  starc3   star3 in conserve mode

Usage:  patch_case.py case <family> <case dir> [--metric mean]
                                                     writes input.ini, MESH/mesh.tec, INPUT/{ic.tec,bc.txt,phase.txt,
                                                     properties.dat}, patch.json; --metric mean sets
                                                     diffusive-metric = mean in input.ini (default: the key is absent)
        patch_case.py check <case dir> [tol]         reads OUTPUT/diagnostic.tec and OUTPUT/wall.tec of the run in the case
                                                     dir; prints one line per group; exit 1 when a group exceeds tol (1e-11)
"""
import json
import math
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fusstec import read_tec  # noqa: E402

A, G = 1000.0, (200.0, 300.0, 150.0)       # T = A + g . x  [K], x in metres (the meshes are unit boxes)
K, RHO, CP = 52.0, 7850.0, 460.0
H, EPS, SIGMA = 100.0, 0.8, 5.67e-8
LZ = 0.01                                   # thickness of the 2-D cases
SEED = 20261023

INI = """[FUSS-Parameters]
newrun = true
res-threshold = 1e-30
iter-threshold = 1

[FUSS-Numerics]
time-scheme = euler
vnn = 0.5
time-accurate = false
integration-variables = cons

[FUSS-IO]
sol-diter     = 1
sol-overwrite = true
sol-format    = tecplot ascii
res-diter = 1
shell-diter = 1
ini-diter = 100000000
"""


def e(v):
    return ' %.15E' % v


def iv(*v):
    return ''.join('%8d' % x for x in v)


def rv(*v):
    return ''.join('    %.15E,' % x for x in v)


def exact(p, g):
    return A + g[0] * p[0] + g[1] * p[1] + g[2] * p[2]


def tref_305(tw, q):
    x = tw
    for _ in range(80):
        x -= (H * (x - tw) + SIGMA * EPS * (x ** 4 - tw ** 4) - q) / (H + 4 * SIGMA * EPS * x ** 3)
    return x


# ----------------------------------------------------------------------------------------------------------- meshes
class Block:
    """Structured block: nodes[k][j][i] = (x, y, z), cells nx x ny x nz; faces: dict face -> (type, data) or 101 target."""

    def __init__(self, name, nx, ny, nz, fn):
        self.name, self.nx, self.ny, self.nz = name, nx, ny, nz
        self.nodes = [[[list(fn(i, j, k)) for i in range(nx + 1)] for j in range(ny + 1)] for k in range(nz + 1)]
        self.faces = {}   # face -> ('wall', type) | ('sym',) | ('null',) | ('conn', other block number, other face)
        self.partner = {}  # (face, i, j, k) -> (other block, i, j, k, other face, d11, d12, d21, d22) where the default does not hold

    def perturb(self, amp, rng, fixed_i=(), fixed_j=(), fixed_k=()):
        """Move the interior nodes by up to amp of the local cell size (boundary and listed index planes fixed)."""
        nx, ny, nz = self.nx, self.ny, self.nz
        # a 2-D case moves its two node planes together, once (a loop over both planes moved every node twice, up to
        # 0.51 of the cell size, and made a few cells non-convex: found by the metrics chat, 2026-10-23)
        for k in (range(nz + 1) if nz > 1 else (0,)):
            for j in range(ny + 1):
                for i in range(nx + 1):
                    if i in (0, nx) or i in fixed_i or j in (0, ny) or j in fixed_j:
                        continue
                    if nz > 1 and (k in (0, nz) or k in fixed_k):
                        continue
                    n = self.nodes[k][j][i]
                    hx = 0.5 * (self.nodes[k][j][i + 1][0] - self.nodes[k][j][i - 1][0])
                    hy = 0.5 * (self.nodes[k][j + 1][i][1] - self.nodes[k][j - 1][i][1])
                    dx, dy = amp * hx * rng.uniform(-1, 1), amp * hy * rng.uniform(-1, 1)
                    dz = 0.0
                    if nz > 1:
                        hz = 0.5 * (self.nodes[k + 1][j][i][2] - self.nodes[k - 1][j][i][2])
                        dz = amp * hz * rng.uniform(-1, 1)
                    if nz == 1:   # the two node planes of a 2-D case move together
                        for kk in (0, 1):
                            self.nodes[kk][j][i][0] += dx
                            self.nodes[kk][j][i][1] += dy
                    else:
                        n[0] += dx; n[1] += dy; n[2] += dz

    def face_cells(self, f):
        """Cells (i, j, k) 1-based on face f, in FUSS record order."""
        nx, ny, nz = self.nx, self.ny, self.nz
        if f in (1, 2):
            return [(1 if f == 1 else nx, j, k) for k in range(1, nz + 1) for j in range(1, ny + 1)]
        if f in (3, 4):
            return [(i, 1 if f == 3 else ny, k) for k in range(1, nz + 1) for i in range(1, nx + 1)]
        return [(i, j, 1 if f == 5 else nz) for j in range(1, ny + 1) for i in range(1, nx + 1)]

    def face_nodes(self, f, i, j, k):
        """The four nodes of face f of cell (i, j, k) (1-based cell), node indices 0-based."""
        i0, j0, k0 = i - 1, j - 1, k - 1
        if f == 1:
            q = [(i0, j0, k0), (i0, j0 + 1, k0), (i0, j0 + 1, k0 + 1), (i0, j0, k0 + 1)]
        elif f == 2:
            q = [(i0 + 1, j0, k0), (i0 + 1, j0 + 1, k0), (i0 + 1, j0 + 1, k0 + 1), (i0 + 1, j0, k0 + 1)]
        elif f == 3:
            q = [(i0, j0, k0), (i0 + 1, j0, k0), (i0 + 1, j0, k0 + 1), (i0, j0, k0 + 1)]
        elif f == 4:
            q = [(i0, j0 + 1, k0), (i0 + 1, j0 + 1, k0), (i0 + 1, j0 + 1, k0 + 1), (i0, j0 + 1, k0 + 1)]
        elif f == 5:
            q = [(i0, j0, k0), (i0 + 1, j0, k0), (i0 + 1, j0 + 1, k0), (i0, j0 + 1, k0)]
        else:
            q = [(i0, j0, k0 + 1), (i0 + 1, j0, k0 + 1), (i0 + 1, j0 + 1, k0 + 1), (i0, j0 + 1, k0 + 1)]
        return [self.nodes[c][b][a] for a, b, c in q]

    def face_centre_normal(self, f, i, j, k):
        """Centre, outward unit normal and area of face f of cell (i, j, k)."""
        p = self.face_nodes(f, i, j, k)
        c = [sum(q[d] for q in p) / 4 for d in range(3)]
        d1 = [p[2][d] - p[0][d] for d in range(3)]
        d2 = [p[3][d] - p[1][d] for d in range(3)]
        n = [d1[1] * d2[2] - d1[2] * d2[1], d1[2] * d2[0] - d1[0] * d2[2], d1[0] * d2[1] - d1[1] * d2[0]]
        area = 0.5 * math.sqrt(sum(v * v for v in n))
        n = [v / (2 * area) for v in n]
        # orient outward: towards the face from the cell centre
        cc = self.cell_centre(i, j, k)
        if sum((c[d] - cc[d]) * n[d] for d in range(3)) < 0:
            n = [-v for v in n]
        return c, n, area

    def cell_centre(self, i, j, k):
        pts = [self.nodes[k - 1 + c][j - 1 + b][i - 1 + a] for c in (0, 1) for b in (0, 1) for a in (0, 1)]
        return [sum(q[d] for q in pts) / 8 for d in range(3)]

    def cell_area(self, i, j, k, faces):
        return sum(self.face_centre_normal(f, i, j, k)[2] for f in faces)

    def cell_volume(self, i, j, k):
        """Volume from the six faces: V = (1/3) sum_f A_f (c_f . n_f), the divergence of x over the closed surface."""
        v = 0.0
        for f in range(1, 7):
            c, n, a = self.face_centre_normal(f, i, j, k)
            v += a * sum(c[d] * n[d] for d in range(3))
        return v / 3.0

    def nodes_flat(self):
        return [n[d] for d in range(3) for plane in self.nodes for row in plane for n in row]


# ---------------------------------------------------------------------------------------- Part 2 families (2026-10-28)
EXTRA = {}          # family -> {'mode': 'linear' (default) | 'conserve'}
LREF2 = 2.0 ** -7   # dyadic unit of jump2, fan, fan4


def mesh_star(m, n):
    """MOSE mesh_star: a regular m-gon in m quadrilateral blocks meeting at its centre C; block b has (xi, eta) = (0, 0) at C,
    (1, 0) at the middle of the side before vertex b, (1, 1) at vertex b, (0, 1) at the middle of the side after it, so that the
    eta = 0 side of block b+1 is the xi = 0 side of block b (the same nodes). Returns a list of (n+1)x(n+1) arrays of (x, y)."""
    ang = [2.0 * math.pi * b / m + 0.3 for b in range(m)]
    V = [(math.cos(a), math.sin(a)) for a in ang]
    E = [(0.5 * (V[b][0] + V[(b + 1) % m][0]), 0.5 * (V[b][1] + V[(b + 1) % m][1])) for b in range(m)]
    t = [i / float(n) for i in range(n + 1)]
    spoke = [[(tt * E[b][0], tt * E[b][1]) for tt in t] for b in range(m)]      # C -> E_b, shared by blocks b and b+1
    out = []
    for b in range(m):
        Eb, Em = E[b], E[(b - 1) % m]
        XY = [[None] * (n + 1) for _ in range(n + 1)]
        for i in range(n + 1):
            for j in range(n + 1):
                x, y = t[i], t[j]
                XY[i][j] = tuple(x * (1 - y) * Em[d] + x * y * V[b][d] + (1 - x) * y * Eb[d] for d in (0, 1))
        for q in range(n + 1):
            XY[q][0] = spoke[(b - 1) % m][q]       # eta = 0: C -> E_{b-1}
            XY[0][q] = spoke[b][q]                 # xi = 0: C -> E_b
        out.append(XY)
    return out


def _key(nodes):
    return tuple(sorted(tuple(round(c, 9) for c in p) for p in nodes))


def tangential(face):
    return {1: (1, 2), 2: (1, 2), 3: (0, 2), 4: (0, 2), 5: (0, 1), 6: (0, 1)}[face]


def find_connections(blocks):
    """Fill Block.partner for every cell face whose kind is 'conn' by node coincidence, with the index map d of MOSE
    `orientation`: the step of the source cell in the two tangential directions of the source face, for a step of +1 of the
    master cell along its first and second tangential direction (identity (1 0 0 1) on the pairs 2 <-> 1 and 4 <-> 3)."""
    seen = {}
    for bn, b in enumerate(blocks, start=1):
        for f in range(1, 7):
            if b.faces[f][0] != 'conn':
                continue
            for (i, j, k) in b.face_cells(f):
                seen.setdefault(_key(b.face_nodes(f, i, j, k)), []).append((bn, f, i, j, k))
    conn = {}
    for v in seen.values():
        if len(v) == 2:
            conn[v[0]] = v[1]; conn[v[1]] = v[0]
        elif len(v) > 2:
            sys.exit('face shared by more than two cells')
        else:
            sys.exit('connection face without partner: %s' % (v,))
    for rec, (bs, fs, is_, js, ks) in conn.items():
        bn, f, i, j, k = rec
        b = blocks[bn - 1]
        dims = (b.nx, b.ny, b.nz)
        out = []
        for t in tangential(f):
            idx = [i, j, k]
            lim = dims[t]
            if lim == 1:
                out.append((0, 1) if fs in (1, 2, 3, 4) else (1, 0))
                continue
            step = 1 if idx[t] < lim else -1
            idx[t] += step
            nb = conn.get((bn, f, idx[0], idx[1], idx[2]))
            if nb is None or nb[0] != bs or nb[1] != fs:
                sys.exit('the connection does not continue along a tangential direction')
            ds = [nb[2] - is_, nb[3] - js, nb[4] - ks]
            ts = tangential(fs)
            out.append((step * ds[ts[0]], step * ds[ts[1]]))
        b.partner[(f, i, j, k)] = (bs, is_, js, ks, fs, out[0][0], out[0][1], out[1][0], out[1][1])


def geom_stretch(n, ratio):
    """n+1 nodes on [0, 1], the last cell `ratio` times the first, geometric."""
    r = ratio ** (1.0 / (n - 1)) if n > 1 else 1.0
    w = [r ** m for m in range(n)]
    cum = [0.0]
    for v in w:
        cum.append(cum[-1] + v)
    return [c / cum[-1] for c in cum]


def build_part2(family, rng):
    """Blocks and gradient of the Part 2 families (None when `family` is not one)."""
    g = [G[0], G[1], 0.0]
    EXTRA[family] = {'mode': 'linear'}
    if family == 'jump2':
        n = 8
        sizes = [1.0 if m % 2 == 0 else 2.0 for m in range(n)]
        cum = [0.0]
        for v in sizes:
            cum.append(cum[-1] + v)
        a = [c / 16.0 for c in cum]                                   # units of 1/16, dyadic
        b = Block('Block1', n, n, 1, lambda i, j, k: (LREF2 * (a[i] + 0.5 * a[j]), LREF2 * a[j], LREF2 * (k - 0.5) * 0.0625))
        for f, t in {1: 302, 2: 301, 3: 302, 4: 303}.items():
            b.faces[f] = ('wall', t)
        b.faces[5] = b.faces[6] = ('null',)
        return [b], [v * 128.0 for v in g]
    if family in ('fan', 'fan4'):
        n, top = 8, (7.0 / 64.0 if family == 'fan' else 0.25)
        b = Block('Block1', n, n, 1, lambda i, j, k: (LREF2 * i / float(n), LREF2 * j * (1.0 / 16.0 if i == 0 else top),
                                                       LREF2 * (k - 0.5) * 0.0625))
        for f, t in {1: 302, 2: 301, 3: 302, 4: 303}.items():
            b.faces[f] = ('wall', t)
        b.faces[5] = b.faces[6] = ('null',)
        return [b], [v * 128.0 for v in g]
    if family in ('star3', 'star5', 'starc3'):
        m, n = (5, 8) if family == 'star5' else (3, 8)
        meshes = mesh_star(m, n)
        blocks = []
        for bi, XY in enumerate(meshes):
            blocks.append(Block('Block%d' % (bi + 1), n, n, 1, lambda i, j, k, XY=XY: (XY[i][j][0], XY[i][j][1], (k - 0.5) * 0.02)))
        srng = random.Random(SEED + 5)
        for b in blocks:
            for j in range(1, n):
                for i in range(1, n):
                    dx, dy = 0.15 * (0.5 / n) * srng.uniform(-1, 1), 0.15 * (0.5 / n) * srng.uniform(-1, 1)
                    for kk in (0, 1):
                        b.nodes[kk][j][i][0] += dx
                        b.nodes[kk][j][i][1] += dy
        mix = {0: (302, 302), 1: (301, 302), 2: (303, 301), 3: (302, 305), 4: (301, 301)}
        for bi, b in enumerate(blocks):
            if family == 'starc3':
                b.faces[2] = ('wall', 301); b.faces[4] = ('wall', 301)
            else:
                b.faces[2] = ('wall', mix[bi][0]); b.faces[4] = ('wall', mix[bi][1])
            b.faces[1] = ('conn', (bi + 1) % m + 1, 3)       # xi = 0 side = eta = 0 side of the next block
            b.faces[3] = ('conn', (bi - 1) % m + 1, 1)
            b.faces[5] = b.faces[6] = ('null',)
        find_connections(blocks)
        if family == 'starc3':
            EXTRA[family] = {'mode': 'conserve'}
        return blocks, g
    if family == 'rough3d':
        n = 6
        zs = geom_stretch(n, 5.0)
        slope = [(0.0, 0.0), (0.2, 0.0), (-0.3, 1.2)]            # y = t * (c0 + a s), per block
        blocks = []
        for bi in range(3):
            a, c0 = slope[bi]
            c0 = (1.0 if bi == 0 else (1.0 if bi == 1 else 1.2))
            blocks.append(Block('Block%d' % (bi + 1), n, n, n,
                                lambda i, j, k, bi=bi, a=a, c0=c0: ((bi + i / float(n)) / 3.0, (j / float(n)) * (c0 + a * i / float(n)), zs[k])))
        for b in blocks:
            b.perturb(0.25, rng)
            b.faces.update({3: ('wall', 302), 4: ('wall', 302), 5: ('wall', 302), 6: ('wall', 302)})
        blocks[0].faces[1] = ('wall', 302); blocks[0].faces[2] = ('conn', 2, 1)
        blocks[1].faces[1] = ('conn', 1, 2); blocks[1].faces[2] = ('conn', 3, 1)
        blocks[2].faces[1] = ('conn', 2, 2); blocks[2].faces[2] = ('wall', 302)
        return blocks, list(G)
    if family in ('wall3dA', 'wall3dB', 'sym3d'):
        n = 6
        sh = (0.3, 0.2, 0.25)
        if family == 'sym3d':
            fn = lambda i, j, k: (i / float(n) + sh[0] * k / float(n), j / float(n) + sh[1] * k / float(n), k / float(n))
        else:
            fn = lambda i, j, k: (i / float(n) + sh[0] * j / float(n) + 0.1 * k / float(n),
                                  j / float(n) + sh[1] * k / float(n), k / float(n) + sh[2] * i / float(n))
        b = Block('Block1', n, n, n, fn)
        b.perturb(0.15, rng)
        types = {'wall3dA': [302, 301, 303, 305, 302, 301], 'wall3dB': [301, 303, 305, 302, 303, 305],
                 'sym3d': [302, 301, 303, 305, 'sym', 302]}[family]
        for f in range(1, 7):
            b.faces[f] = ('sym',) if types[f - 1] == 'sym' else ('wall', types[f - 1])
        return [b], ([G[0], G[1], 0.0] if family == 'sym3d' else list(G))
    if family in ('conn14', 'conn14s', 'conn14c', 'conn14p'):
        n = 8
        shear = 0.3
        xw = 0.2                                                   # block of the small cells: x in [0, 0.2], 8 cells of 0.025
        small = lambda i, j, k: (xw * i / n, j / float(n), (k - 0.5) * LZ)
        large = lambda i, j, k: (xw + (1 - xw) * i / n,
                                 (j / float(n)) * (1 + 0.25 * i / n), (k - 0.5) * LZ)   # fold at the interface: j lines rise
        # the interface must coincide: at i = n of `small` and i = 0 of `large` both have y = j / n and x = xw + shear y
        bs_, bl_ = Block('Block1', n, n, 1, small), Block('Block2', n, n, 1, large)
        if family == 'conn14p':     # control: no fold, nodes not moved: every cell of a block has the same metric
            large = lambda i, j, k: (xw + (1 - xw) * i / n, j / float(n), (k - 0.5) * LZ)
            bl_ = Block('Block2', n, n, 1, large)
        else:
            bs_.perturb(0.25, rng); bl_.perturb(0.25, rng)
        for b in (bs_, bl_):                                       # shear x += 0.3 y after the nodes are moved (same map in both)
            for plane in b.nodes:
                for row in plane:
                    for nd in row:
                        nd[0] += shear * nd[1]
        for b in (bs_, bl_):
            b.faces.update({3: ('wall', 302), 4: ('wall', 302), 5: ('null',), 6: ('null',)})
        bs_.faces[1] = ('wall', 302); bs_.faces[2] = ('conn', 2, 1)
        bl_.faces[1] = ('conn', 1, 2); bl_.faces[2] = ('wall', 302)
        blocks = [bs_, bl_]
        if family == 'conn14s':
            # swapped: block 1 = large cells, block 2 = small ones; the large block keeps its geometry, face 1 <-> face 2
            blocks = [bl_, bs_]
            bl_.name, bs_.name = 'Block1', 'Block2'
            bl_.faces[1] = ('conn', 2, 2); bs_.faces[2] = ('conn', 1, 1)
            bl_.faces[2] = ('wall', 302); bs_.faces[1] = ('wall', 302)
        else:
            bs_.name, bl_.name = 'Block1', 'Block2'
        if family == 'conn14c':
            EXTRA[family] = {'mode': 'conserve'}
            for b in blocks:
                for f in (1, 2, 3, 4):
                    if b.faces[f][0] == 'wall':
                        b.faces[f] = ('wall', 301)
        find_connections(blocks)
        return blocks, g
    if family in ('flipA', 'flipB'):
        n = 8
        b1 = Block('Block1', n, n, 1, lambda i, j, k: (0.5 * i / n, j / float(n), (k - 0.5) * LZ))
        if family == 'flipA':      # block 2 with both index axes reversed: face 2 <-> face 2, j flipped
            b2 = Block('Block2', n, n, 1, lambda i, j, k: (1.0 - 0.5 * i / n, 1.0 - j / float(n), (k - 0.5) * LZ))
            b1.faces.update({1: ('wall', 302), 2: ('conn', 2, 2), 3: ('wall', 302), 4: ('wall', 302)})
            b2.faces.update({1: ('wall', 302), 2: ('conn', 1, 2), 3: ('wall', 302), 4: ('wall', 302)})
        else:                      # block 2 with the axes exchanged and i reversed: face 2 <-> face 3
            b2 = Block('Block2', n, n, 1, lambda i, j, k: (0.5 + 0.5 * j / n, 1.0 - i / float(n), (k - 0.5) * LZ))
            b1.faces.update({1: ('wall', 302), 2: ('conn', 2, 3), 3: ('wall', 302), 4: ('wall', 302)})
            b2.faces.update({1: ('wall', 302), 2: ('wall', 302), 3: ('conn', 1, 2), 4: ('wall', 302)})
        amp = 0.2 * 0.5 / n
        for b in (b1, b2):
            for j in range(1, n):
                for i in range(1, n):
                    dx, dy = amp * rng.uniform(-1, 1), amp * rng.uniform(-1, 1)
                    for kk in (0, 1):
                        b.nodes[kk][j][i][0] += dx
                        b.nodes[kk][j][i][1] += dy
            b.faces[5] = b.faces[6] = ('null',)
        find_connections([b1, b2])
        return [b1, b2], g
    return None


PART2 = ('jump2', 'fan', 'fan4', 'star3', 'star5', 'starc3', 'rough3d', 'wall3dA', 'wall3dB', 'sym3d', 'conn14', 'conn14s',
         'conn14c', 'conn14p', 'flipA', 'flipB')


def build(family):
    rng = random.Random(SEED)
    g = list(G)
    blocks = []
    if family in ('para', 'rand', 'jump', 'wallmix', 'sym', 'step'):
        n = 16
        if family == 'step':
            xs = [0.05 * i if i <= 8 else 0.4 + 0.075 * (i - 8) for i in range(n + 1)]
            fn = lambda i, j, k: (xs[i], j / n, (k - 0.5) * LZ)
        elif family == 'para':
            s = math.tan(math.radians(20.0))
            fn = lambda i, j, k: (i / n + s * j / n, j / n, (k - 0.5) * LZ)
        elif family in ('rand', 'wallmix'):
            fn = lambda i, j, k: (i / n, j / n, (k - 0.5) * LZ)
        elif family == 'jump':
            w = [1.0 if m % 2 == 0 else 1.5 for m in range(n)]
            cum = [0.0]
            for v in w:
                cum.append(cum[-1] + v)
            cum = [v / cum[-1] for v in cum]
            fn = lambda i, j, k: (cum[i] + 0.5 * cum[j], cum[j], (k - 0.5) * LZ)
        else:   # sym: plane y = 0 is a symmetry plane, the field has no y gradient
            g = [G[0], 0.0, 0.0]
            fn = lambda i, j, k: (i / n + 0.35 * j / n, j / n, (k - 0.5) * LZ)
        b = Block('Block1', n, n, 1, fn)
        if family in ('rand', 'wallmix'):
            b.perturb(0.25, rng)
        if family == 'sym':
            b.perturb(0.15, rng)
        walls = {1: 302, 2: 302, 3: 302, 4: 302}
        if family == 'wallmix':
            walls = {1: 302, 2: 303, 3: 301, 4: 305}
        for f, t in walls.items():
            b.faces[f] = ('wall', t)
        if family == 'sym':
            b.faces[3] = ('sym',)
        b.faces[5] = b.faces[6] = ('null',)
        blocks = [b]
        g = g if family == 'sym' else [G[0], G[1], 0.0]   # step: the 2-D field of para
    elif family == 'fold2':
        n = 8
        b1 = Block('Block1', n, n, 1, lambda i, j, k: (0.5 * i / n, j / n, (k - 0.5) * LZ))
        b2 = Block('Block2', n, n, 1, lambda i, j, k: (0.5 + 0.5 * i / n, (j / n) * (1 + 0.25 * i / n), (k - 0.5) * LZ))
        b1.perturb(0.25, rng); b2.perturb(0.25, rng)
        for b in (b1, b2):
            b.faces.update({3: ('wall', 302), 4: ('wall', 302), 5: ('null',), 6: ('null',)})
        b1.faces[1] = ('wall', 302); b1.faces[2] = ('conn', 2, 1)
        b2.faces[1] = ('conn', 1, 2); b2.faces[2] = ('wall', 302)
        blocks = [b1, b2]
        g = [G[0], G[1], 0.0]
    elif family == 'box3d':
        n = 6
        b1 = Block('Block1', n, n, n, lambda i, j, k: (0.5 * i / n, j / n, k / n))
        b2 = Block('Block2', n, n, n, lambda i, j, k: (0.5 + 0.5 * i / n, (j / n) * (1 + 0.2 * i / n), k / n))
        b1.perturb(0.15, rng); b2.perturb(0.15, rng)
        for b in (b1, b2):
            b.faces.update({3: ('wall', 302), 4: ('wall', 302), 5: ('wall', 302), 6: ('wall', 302)})
        b1.faces[1] = ('wall', 302); b1.faces[2] = ('conn', 2, 1)
        b2.faces[1] = ('conn', 1, 2); b2.faces[2] = ('wall', 302)
        blocks = [b1, b2]
    elif family in PART2:
        blocks, g = build_part2(family, rng)
    else:
        sys.exit('unknown family ' + family)
    return blocks, g


def tfield(p, g):
    """Initial and wall temperature at a point: the exact linear field."""
    return exact(p, g)


def write_case(family, out, metric=None):
    blocks, g = build(family)
    extra = EXTRA.get(family, {'mode': 'linear'})
    mode = extra['mode']
    trng = random.Random(SEED + 99)
    for sub in ('MESH', 'INPUT', 'OUTPUT'):
        os.makedirs(os.path.join(out, sub), exist_ok=True)
    mesh, ic, bc = [' VARIABLES ="x" "y" "z"\n'], [' VARIABLES ="x" "y" "z" T matID\n'], []
    for bn, b in enumerate(blocks, start=1):
        kk = b.nz + 1
        zone = (' ZONE  T = %s, I=%d, J=%d, K=%d, DATAPACKING=BLOCK, VARLOCATION=([1-3]=NODAL), SOLUTIONTIME=0.0, STRANDID = 0\n'
                % (b.name, b.nx + 1, b.ny + 1, kk))
        flat = b.nodes_flat()
        mesh += [zone, '\n'.join(e(v) for v in flat) + '\n']
        if mode == 'conserve':
            tcell = [1000.0 + 100.0 * trng.uniform(-1, 1) for k in range(1, b.nz + 1) for j in range(1, b.ny + 1)
                     for i in range(1, b.nx + 1)]
        else:
            tcell = [tfield(b.cell_centre(i, j, k), g) for k in range(1, b.nz + 1) for j in range(1, b.ny + 1)
                     for i in range(1, b.nx + 1)]
        ic += [zone.replace(b.name, 'B%d-SP' % bn).replace('VARLOCATION=([1-3]=NODAL)',
                                                           'VARLOCATION=([1-3]=NODAL,[4-5]=CELLCENTERED)'),
               '\n'.join(e(v) for v in flat + tcell + [1.0] * len(tcell)) + '\n']
        for f in range(1, 7):
            kind = b.faces[f]
            for (i, j, k) in b.face_cells(f):
                if kind[0] == 'null':
                    bc.append(iv(bn, i, j, k, f, 0))
                elif kind[0] == 'sym':
                    bc.append(iv(bn, i, j, k, f, 300))
                elif kind[0] == 'conn':
                    ob, of = kind[1], kind[2]
                    if (f, i, j, k) in b.partner:
                        bc += [iv(bn, i, j, k, f, 101), iv(*b.partner[(f, i, j, k)])]
                    else:
                        oi = 1 if of == 1 else blocks[ob - 1].nx
                        bc += [iv(bn, i, j, k, f, 101), iv(ob, oi, j, k, of, 1, 0, 0, 1)]
                else:
                    t = kind[1]
                    c, nrm, _ = b.face_centre_normal(f, i, j, k)
                    tw = tfield(c, g)
                    q = K * sum(g[d] * nrm[d] for d in range(3)) if mode != 'conserve' else 0.0   # flux into the solid
                    data = {302: [tw], 301: [q], 303: [H, 0.0, tw + q / H], 305: [H, EPS, tref_305(tw, q)]}[t]
                    bc += [iv(bn, i, j, k, f, t), rv(*data)]
    open(os.path.join(out, 'MESH', 'mesh.tec'), 'w').write(''.join(mesh))
    open(os.path.join(out, 'INPUT', 'ic.tec'), 'w').write(''.join(ic))
    open(os.path.join(out, 'INPUT', 'bc.txt'), 'w').write('\n'.join(bc) + '\n')
    open(os.path.join(out, 'INPUT', 'phase.txt'), 'w').write('solid-bulk\nmat 1\n')
    with open(os.path.join(out, 'INPUT', 'properties.dat'), 'w') as fh:
        fh.write('TITLE = "Mass Thermodynamic Properties"\n'
                 'VARIABLES = "Temperature", "Cp", "Density", "Conductivity", "Energy"\nZONE T="mat"\nI=2500, F=POINT\n')
        fh.writelines('%.1f %.6f %.6f %.6f %.6f\n' % (t, CP, RHO, K, RHO * CP * t) for t in range(1, 2501))
    ini = INI if metric is None else INI.replace('integration-variables = cons\n',
                                                 'integration-variables = cons\ndiffusive-metric = %s\n' % metric)
    open(os.path.join(out, 'input.ini'), 'w').write(ini)
    meta = {'family': family, 'g': g, 'k': K, 'mode': mode, 'metric': metric or 'stencil',
            'blocks': [{'name': b.name, 'n': [b.nx, b.ny, b.nz],
             'faces': {str(f): list(b.faces[f]) for f in range(1, 7)}} for b in blocks]}
    json.dump(meta, open(os.path.join(out, 'patch.json'), 'w'), indent=1)
    print('%s: %d block(s), %d boundary records' % (family, len(blocks), sum(1 for l in bc if len(l.split()) == 6)))


def check_conserve(meta, blocks, diag, tol):
    """Random T, adiabatic walls: the sum over all cells of the flux residual R = -E V / dt must vanish (every internal flux,
    the connection ones included, enters two cells with opposite signs)."""
    tot, absum, n = 0.0, 0.0, 0
    per = []
    for bn, (b, (title, z)) in enumerate(zip(blocks, diag.items()), start=1):
        E, DT = z['E'], z['dt']
        bt = 0.0
        for k in range(1, b.nz + 1):
            for j in range(1, b.ny + 1):
                for i in range(1, b.nx + 1):
                    idx = ((k - 1) * b.ny + (j - 1)) * b.nx + (i - 1)
                    r = -E[idx] * b.cell_volume(i, j, k) / DT[idx]
                    tot += r; absum += abs(r); bt += r; n += 1
        per.append(bt)
    rel = abs(tot) / absum
    ok = rel <= tol
    print('CONSERVE %-8s sum R / sum |R| %.6e  %d cells, per-block sums %s  %s'
          % (meta['family'], rel, n, ' '.join('%.3e' % v for v in per), 'ok' if ok else 'FAIL'))
    return 0 if ok else 1


def check_case(case, tol):
    meta = json.load(open(os.path.join(case, 'patch.json')))
    blocks, g = build(meta['family'])
    gnorm = math.sqrt(sum(v * v for v in g))
    diag = read_tec(os.path.join(case, 'OUTPUT', 'diagnostic.tec'))
    groups = {}
    status = 0
    if meta.get('mode') == 'conserve':
        return check_conserve(meta, blocks, diag, tol)
    # diagnostic.tec holds, after the Euler step, the scaled update E = (-R/V + qvol) dt of every cell (Lib_Newstate:
    # Scale_Residual_Cons), with R the sum of the face fluxes [W]: the flux sum is |E| V / dt (qvol = 0 here)
    for bn, (b, (title, z)) in enumerate(zip(blocks, diag.items()), start=1):
        E, DT = z['E'], z['dt']
        active = [1, 2, 3, 4] + ([5, 6] if b.nz > 1 else [])
        for k in range(1, b.nz + 1):
            for j in range(1, b.ny + 1):
                for i in range(1, b.nx + 1):
                    tags = []
                    for f in active:
                        on = (f == 1 and i == 1) or (f == 2 and i == b.nx) or (f == 3 and j == 1) or \
                             (f == 4 and j == b.ny) or (f == 5 and k == 1) or (f == 6 and k == b.nz)
                        if on:
                            kind = b.faces[f]
                            tags.append({'wall': 'wall%d' % kind[1] if kind[0] == 'wall' else '',
                                         'conn': 'conn', 'sym': 'sym'}.get(kind[0], kind[0]) or 'wall')
                    grp = 'interior' if not tags else '+'.join(sorted(set(tags)))
                    if meta['family'].startswith('star') and i == 1 and j == 1:
                        grp = 'centre'       # the cells at the point shared by all blocks
                    idx = ((k - 1) * b.ny + (j - 1)) * b.nx + (i - 1)
                    flux_sum = abs(E[idx]) * b.cell_volume(i, j, k) / DT[idx]
                    val = flux_sum / (K * gnorm * b.cell_area(i, j, k, active))
                    groups.setdefault(grp, [0.0, 0])
                    groups[grp][0] = max(groups[grp][0], val)
                    groups[grp][1] += 1
    for grp in sorted(groups):
        v, n = groups[grp]
        ok = v <= tol
        status |= (not ok)
        print('RESID %-8s %-22s %.6e  %4d cells  %s' % (meta['family'], grp, v, n, 'ok' if ok else 'FAIL'))
    # wall output: qw of every wall face against k g . n_out
    wall_path = os.path.join(case, 'OUTPUT', 'wall.tec')
    if os.path.exists(wall_path):
        wz = read_tec(wall_path)
        worst = {}
        for title, z in wz.items():
            bn = int(title[1:title.index('F')]); f = int(title[title.index('F') + 1:])
            b = blocks[bn - 1]
            cells = b.face_cells(f)
            if b.faces[f][0] != 'wall':
                continue
            t = b.faces[f][1]
            qw = z['qw']
            # wall zone cells: for faces 1/2 index (1, j, k) -> j fastest; 3/4 (i, 1, k) -> i fastest; 5/6 (i, j, 1)
            for m, (i, j, k) in enumerate(cells):
                c, nrm, _ = b.face_centre_normal(f, i, j, k)
                q = K * sum(g[d] * nrm[d] for d in range(3))
                d = abs(qw[m] - q) / (K * gnorm)
                key = 'wall%d' % t
                worst[key] = max(worst.get(key, 0.0), d)
        for key in sorted(worst):
            ok = worst[key] <= tol
            status |= (not ok)
            print('QWALL %-8s %-22s %.6e  %s' % (meta['family'], key, worst[key], 'ok' if ok else 'FAIL'))
    return status


if __name__ == '__main__':
    if len(sys.argv) >= 4 and sys.argv[1] == 'case':
        write_case(sys.argv[2], sys.argv[3], sys.argv[5] if len(sys.argv) >= 6 and sys.argv[4] == '--metric' else None)
    elif len(sys.argv) >= 3 and sys.argv[1] == 'check':
        sys.exit(check_case(sys.argv[2], float(sys.argv[3]) if len(sys.argv) > 3 else 1e-11))
    else:
        sys.exit(__doc__)
