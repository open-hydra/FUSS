#!/usr/bin/env python3
"""PeriodicSector: periodic faces (201) against the full ring they stand for.

  make_sector.py case <kind> <case dir> [--metric mean]

The geometry is a sector of an annular wall about the x axis, NS = 6 sectors in the full ring, twisted into a helix
(theta = theta0 + ALPHA x, as the side faces of the channels of a cooling jacket), with the interior nodes moved at
random (fixed seed): 6 x 6 x 4 cells (x, theta, r) per sector, cut in two blocks of 3 cells in theta (A and B). The
nodes of the two side faces of a sector (theta = 0 and theta = 360/NS deg) are moved only along the face, by the same
amounts, so that the two faces are images of each other by the rotation of one sector about x.

Every sector has the same boundary values, functions of the position in its sector (a pattern that repeats from sector
to sector and is not symmetric about the middle plane of the sector): walls 302 at x = 0 and x = L, a convective and
radiative hot wall 303 at the inner radius, a wall 301 at the outer radius; the initial field is a periodic,
non-symmetric function of the same kind. RK3, ITERS iterations.

Kinds:
  ring     the full ring, 2 NS blocks (A and B of every sector) joined by block connections (101), the last sector to the
           first: the reference of `sector`
  ring1    the full ring, NS blocks (one per sector) joined by 101: the reference of `sector1`
  sector   one sector, blocks A and B joined by a 101, face 3 of A periodic with face 4 of B (201): the side faces
  sector1  one sector as one block of 6 cells in theta, its faces 3 and 4 periodic with each other (201)
  sym      as sector with symmetry (300) on its side faces instead: not the ring (the pattern is not symmetric)
  sym1     as sector1 with symmetry on its side faces
  flip     as sector with the cells of the periodic faces paired in reverse order along x: the faces are not images of
           each other with that index map, and the run must stop at set-up

A periodic face gives the ghost cells the cells of the other face, moved by the rotation of the pair: the fields of
`sector` and `sector1` must be those of every sector of `ring` and `ring1` to round-off (compare_sector.py), for both
metrics. Each is compared with the ring cut in the same blocks: on this mesh a sector cut in two blocks is not the same
discretisation as one block, with or without periodic faces (sym against sym1: 4.8e-4 of T at convergence, the same
with the code before 201; the largest differences next to the walls at the cut, where the wall ghost cells are filled
from differences inside the block).
"""
import math
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'patch_test'))
import patch_case as pc  # noqa: E402

NS = 6                      # sectors of the full ring
NX, NT, NR = 6, 6, 4        # cells per sector along x, theta, r
NTA = 3                     # cells in theta of block A (B takes the others)
L, R0, R1 = 0.05, 0.05, 0.06
ALPHA = 3.0                 # helix: theta grows by ALPHA x [rad/m]
DTH = 2.0 * math.pi / NS
SEED = 20261116
ITERS = 100

INI = """[FUSS-Parameters]
newrun = true
res-threshold = 1e-30
iter-threshold = {iters}

[FUSS-Numerics]
time-scheme = RK3
vnn = 0.5
time-accurate = false
integration-variables = cons
{metric}
[FUSS-IO]
sol-diter     = {iters}
sol-overwrite = true
sol-format    = tecplot ascii
res-diter = 1
shell-diter = 10
ini-diter = 100000000
"""


def base_nodes():
    """Nodes of one sector in base coordinates (x, r, s), s in [0, 1] the fraction of the sector: [k][j][i]."""
    rng = random.Random(SEED)
    nodes = [[[None] * (NX + 1) for _ in range(NT + 1)] for _ in range(NR + 1)]
    side = {}     # the moves of the nodes of the side faces: the same on s = 0 and s = 1
    for k in range(NR + 1):
        for j in range(NT + 1):
            for i in range(NX + 1):
                x = L * i / NX
                r = R0 + (R1 - R0) * k / NR
                s = j / NT
                ix, kr = 0 < i < NX, 0 < k < NR
                if j in (0, NT):
                    if (i, k) not in side:
                        side[(i, k)] = ((rng.uniform(-0.2, 0.2) if ix else 0.0), (rng.uniform(-0.2, 0.2) if kr else 0.0))
                    dx, dr = side[(i, k)]
                    ds = 0.0
                else:
                    dx = rng.uniform(-0.2, 0.2) if ix else 0.0
                    dr = rng.uniform(-0.2, 0.2) if kr else 0.0
                    ds = rng.uniform(-0.2, 0.2)
                nodes[k][j][i] = (x + dx * L / NX, r + dr * (R1 - R0) / NR, s + ds / NT)
    return nodes


def cart(p, sector):
    x, r, s = p
    th = (sector + s) * DTH + ALPHA * x
    return (x, r * math.cos(th), -r * math.sin(th))  # ours: z -> -z, a right-handed block (MORPH refuses a left-handed one)


# Boundary and initial values: functions of the base position (x, r, s) of a face or cell centre, periodic in s
def t_end(x, r, s):
    return 450.0 + 80.0 * math.sin(2 * math.pi * s) + 30.0 * math.cos(2 * math.pi * s) + 1000.0 * (r - R0)


def tref_hot(x, r, s):
    return 1500.0 + 300.0 * math.sin(2 * math.pi * s + 0.7) + 2000.0 * x


def t_init(x, r, s):
    return 500.0 + 60.0 * math.sin(2 * math.pi * s + 0.3) * math.cos(40.0 * x) + 3000.0 * (r - R0)


def mean(pts):
    return tuple(sum(p[c] for p in pts) / len(pts) for c in range(3))


class Blk:
    def __init__(self, name, sector, j0, nt, base):
        self.name, self.sector, self.j0, self.nt = name, sector, j0, nt
        self.base = [[[base[k][j0 + j][i] for i in range(NX + 1)] for j in range(nt + 1)] for k in range(NR + 1)]
        self.faces = {}

    def cell_base(self, i, j, k):      # 1-based cell indices
        return mean([self.base[k - 1 + c][j - 1 + b][i - 1 + a] for a in (0, 1) for b in (0, 1) for c in (0, 1)])

    def face_base(self, f, i, j, k):
        if f in (1, 2):
            a = i - 1 if f == 1 else i
            pts = [self.base[k - 1 + c][j - 1 + b][a] for b in (0, 1) for c in (0, 1)]
        elif f in (3, 4):
            b = j - 1 if f == 3 else j
            pts = [self.base[k - 1 + c][b][i - 1 + a] for a in (0, 1) for c in (0, 1)]
        else:
            c = k - 1 if f == 5 else k
            pts = [self.base[c][j - 1 + b][i - 1 + a] for a in (0, 1) for b in (0, 1)]
        return mean(pts)

    def face_cells(self, f):
        n = (NX, self.nt, NR)
        rng = [range(1, n[0] + 1), range(1, n[1] + 1), range(1, n[2] + 1)]
        d = (f - 1) // 2
        rng[d] = [1 if f % 2 else n[d]]
        return [(i, j, k) for k in rng[2] for j in rng[1] for i in rng[0]]


def blocks_of(kind):
    base = base_nodes()
    if kind == 'ring':
        blks = []
        for s in range(NS):
            blks += [Blk('S%dA' % s, s, 0, NTA, base), Blk('S%dB' % s, s, NTA, NT - NTA, base)]
        for s in range(NS):
            a, b = 2 * s + 1, 2 * s + 2
            nxt = (2 * ((s + 1) % NS)) + 1
            blks[a - 1].faces[4] = ('conn', b, 3, 101)
            blks[b - 1].faces[3] = ('conn', a, 4, 101)
            blks[b - 1].faces[4] = ('conn', nxt, 3, 101)
            blks[nxt - 1].faces[3] = ('conn', b, 4, 101)
    elif kind == 'ring1':
        blks = [Blk('S%d' % s, s, 0, NT, base) for s in range(NS)]
        for s in range(NS):
            nxt = (s + 1) % NS
            blks[s].faces[4] = ('conn', nxt + 1, 3, 101)
            blks[nxt].faces[3] = ('conn', s + 1, 4, 101)
    elif kind in ('sector', 'sym', 'flip'):
        blks = [Blk('A', 0, 0, NTA, base), Blk('B', 0, NTA, NT - NTA, base)]
        blks[0].faces[4] = ('conn', 2, 3, 101)
        blks[1].faces[3] = ('conn', 1, 4, 101)
        if kind in ('sector', 'flip'):
            blks[0].faces[3] = ('conn', 2, 4, 201)
            blks[1].faces[4] = ('conn', 1, 3, 201)
        else:
            blks[0].faces[3] = blks[1].faces[4] = ('sym',)
    elif kind in ('sector1', 'sym1'):
        blks = [Blk('A', 0, 0, NT, base)]
        if kind == 'sector1':
            blks[0].faces[3] = ('conn', 1, 4, 201)
            blks[0].faces[4] = ('conn', 1, 3, 201)
        else:
            blks[0].faces[3] = blks[0].faces[4] = ('sym',)
    else:
        sys.exit(__doc__)
    for b in blks:
        b.faces[1] = b.faces[2] = ('end',)
        b.faces[5] = ('hot',)
        b.faces[6] = ('outer',)
    return blks


def write_case(kind, out, metric):
    blks = blocks_of(kind)
    for sub in ('MESH', 'INPUT', 'OUTPUT'):
        os.makedirs(os.path.join(out, sub), exist_ok=True)
    mesh, ic, bc = [' VARIABLES ="x" "y" "z"\n'], [' VARIABLES ="x" "y" "z" T matID\n'], []
    for bn, b in enumerate(blks, start=1):
        zone = (' ZONE  T = %s, I=%d, J=%d, K=%d, DATAPACKING=BLOCK, VARLOCATION=([1-3]=NODAL), SOLUTIONTIME=0.0, STRANDID = 0\n'
                % (b.name, NX + 1, b.nt + 1, NR + 1))
        xyz = [[cart(b.base[k][j][i], b.sector)[c] for k in range(NR + 1) for j in range(b.nt + 1) for i in range(NX + 1)]
               for c in range(3)]
        flat = xyz[0] + xyz[1] + xyz[2]
        mesh += [zone, '\n'.join(pc.e(v) for v in flat) + '\n']
        tcell = [t_init(*b.cell_base(i, j, k)) for k in range(1, NR + 1) for j in range(1, b.nt + 1) for i in range(1, NX + 1)]
        ic += [zone.replace(' T = %s,' % b.name, ' T = B%d-SP,' % bn)
                   .replace('VARLOCATION=([1-3]=NODAL)', 'VARLOCATION=([1-3]=NODAL,[4-5]=CELLCENTERED)'),
               '\n'.join(pc.e(v) for v in flat + tcell + [1.0] * len(tcell)) + '\n']
        for f in range(1, 7):
            kf = b.faces[f]
            for (i, j, k) in b.face_cells(f):
                fb = b.face_base(f, i, j, k)
                if kf[0] == 'conn':
                    ob, of, code = kf[1], kf[2], kf[3]
                    oj = 1 if of == 3 else blks[ob - 1].nt
                    if kind == 'flip' and code == 201:     # i reversed on the other face: d11 = -1
                        bc += [pc.iv(bn, i, j, k, f, code), pc.iv(ob, NX + 1 - i, oj, k, of, -1, 0, 0, 1)]
                    else:
                        bc += [pc.iv(bn, i, j, k, f, code), pc.iv(ob, i, oj, k, of, 1, 0, 0, 1)]
                elif kf[0] == 'sym':
                    bc.append(pc.iv(bn, i, j, k, f, 300))
                elif kf[0] == 'end':
                    bc += [pc.iv(bn, i, j, k, f, 302), pc.rv(t_end(*fb))]
                elif kf[0] == 'hot':
                    bc += [pc.iv(bn, i, j, k, f, 303), pc.rv(800.0, 2.0e5, tref_hot(*fb))]
                else:
                    bc += [pc.iv(bn, i, j, k, f, 301), pc.rv(-3.0e4)]
    open(os.path.join(out, 'MESH', 'mesh.tec'), 'w').write(''.join(mesh))
    open(os.path.join(out, 'INPUT', 'ic.tec'), 'w').write(''.join(ic))
    open(os.path.join(out, 'INPUT', 'bc.txt'), 'w').write('\n'.join(bc) + '\n')
    open(os.path.join(out, 'INPUT', 'phase.txt'), 'w').write('solid-bulk\nmat 1\n')
    with open(os.path.join(out, 'INPUT', 'properties.dat'), 'w') as fh:
        fh.write('TITLE = "Mass Thermodynamic Properties"\n'
                 'VARIABLES = "Temperature", "Cp", "Density", "Conductivity", "Energy"\nZONE T="mat"\nI=3000, F=POINT\n')
        fh.writelines('%.1f %.6f %.6f %.6f %.6f\n' % (t, pc.CP, pc.RHO, pc.K, pc.RHO * pc.CP * t) for t in range(1, 3001))
    open(os.path.join(out, 'input.ini'), 'w').write(
        INI.format(iters=ITERS, metric='diffusive-metric = mean\n' if metric == 'mean' else ''))


if __name__ == '__main__':
    args = sys.argv[1:]
    metric = 'stencil'
    if '--metric' in args:
        n = args.index('--metric')
        metric = args[n + 1]
        del args[n:n + 2]
    if len(args) == 3 and args[0] == 'case':
        write_case(args[1], args[2], metric)
    else:
        sys.exit(__doc__)
