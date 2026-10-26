#!/usr/bin/env python3
"""Decks of the MetricChecks test (no ATLAS): make_case.py <part> <dir> writes input.ini, MESH/mesh.tec and INPUT/ of one part.
Standard library only. Every deck is a few hundred cells of constant properties (k 52, rho 7850, cp 460), T0 1000 K.

  thin3d    two Cartesian 3-D blocks, block 1 one cell thick in k, joined to block 2 on its k face; walls 302 with T linear in z
            (the lateral faces of the two blocks at the face centre value), steady: the field is exactly linear in z.
  mixed2d   a plane block and an axisymmetric wedge block (k = 1 plane turned by 5 degrees about the x axis) in one 2-D mesh
  fold2d    a 2-D block with one node moved across the next cell
  fold3d    a 3-D block with one node moved across the next cell
  chimturn  two plane blocks joined by chimera faces (102, one donor per ghost cell); the edges along the face of block 1 grow
            4 times across its first cell at the chimera face, so that the first ghost layer turns over
  chimwall  the block 1 of chimturn alone, a wall (302) where the chimera face was
"""
import math
import os
import sys

K, RHO, CP, T0 = 52.0, 7850.0, 460.0, 1000.0
LZ = 0.01


def e(v):
    return ' %.15E' % v


def iv(*v):
    return ''.join('%8d' % x for x in v)


def rv(*v):
    return ''.join('    %.15E,' % x for x in v)


class Block:
    def __init__(self, name, nx, ny, nz, fn):
        self.name, self.nx, self.ny, self.nz = name, nx, ny, nz
        self.nodes = [[[list(fn(i, j, k)) for i in range(nx + 1)] for j in range(ny + 1)] for k in range(nz + 1)]
        self.faces = {}      # face -> ('wall', fn(face centre) -> T) | ('null',) | ('conn', other block, other face) | ('chim', other block)

    def face_cells(self, f):
        nx, ny, nz = self.nx, self.ny, self.nz
        if f in (1, 2):
            return [(1 if f == 1 else nx, j, k) for k in range(1, nz + 1) for j in range(1, ny + 1)]
        if f in (3, 4):
            return [(i, 1 if f == 3 else ny, k) for k in range(1, nz + 1) for i in range(1, nx + 1)]
        return [(i, j, 1 if f == 5 else nz) for j in range(1, ny + 1) for i in range(1, nx + 1)]

    def face_centre(self, f, i, j, k):
        i0, j0, k0 = i - 1, j - 1, k - 1
        q = {1: [(i0, j0, k0), (i0, j0 + 1, k0), (i0, j0 + 1, k0 + 1), (i0, j0, k0 + 1)],
             2: [(i0 + 1, j0, k0), (i0 + 1, j0 + 1, k0), (i0 + 1, j0 + 1, k0 + 1), (i0 + 1, j0, k0 + 1)],
             3: [(i0, j0, k0), (i0 + 1, j0, k0), (i0 + 1, j0, k0 + 1), (i0, j0, k0 + 1)],
             4: [(i0, j0 + 1, k0), (i0 + 1, j0 + 1, k0), (i0 + 1, j0 + 1, k0 + 1), (i0, j0 + 1, k0 + 1)],
             5: [(i0, j0, k0), (i0 + 1, j0, k0), (i0 + 1, j0 + 1, k0), (i0, j0 + 1, k0)],
             6: [(i0, j0, k0 + 1), (i0 + 1, j0, k0 + 1), (i0 + 1, j0 + 1, k0 + 1), (i0, j0 + 1, k0 + 1)]}[f]
        p = [self.nodes[c][b][a] for a, b, c in q]
        return [sum(v[d] for v in p) / 4 for d in range(3)]

    def cell_centre(self, i, j, k):
        pts = [self.nodes[k - 1 + c][j - 1 + b][i - 1 + a] for c in (0, 1) for b in (0, 1) for a in (0, 1)]
        return [sum(q[d] for q in pts) / 8 for d in range(3)]


def write_case(out, blocks, ini, tfield=None):
    for sub in ('MESH', 'INPUT', 'OUTPUT'):
        os.makedirs(os.path.join(out, sub), exist_ok=True)
    mesh, ic, bc = [' VARIABLES ="x" "y" "z"\n'], [' VARIABLES ="x" "y" "z" T matID\n'], []
    for bn, b in enumerate(blocks, start=1):
        zone = (' ZONE  T = %s, I=%d, J=%d, K=%d, DATAPACKING=BLOCK, VARLOCATION=([1-3]=NODAL), SOLUTIONTIME=0.0, STRANDID = 0\n'
                % (b.name, b.nx + 1, b.ny + 1, b.nz + 1))
        flat = [n[d] for d in range(3) for plane in b.nodes for row in plane for n in row]
        mesh += [zone, '\n'.join(e(v) for v in flat) + '\n']
        nc = b.nx * b.ny * b.nz
        tcell = [T0] * nc
        ic += [zone.replace(b.name, 'B%d-SP' % bn).replace('VARLOCATION=([1-3]=NODAL)',
                                                           'VARLOCATION=([1-3]=NODAL,[4-5]=CELLCENTERED)'),
               '\n'.join(e(v) for v in flat + tcell + [1.0] * nc) + '\n']
        for f in range(1, 7):
            kind = b.faces.get(f, ('null',))
            for (i, j, k) in b.face_cells(f):
                if kind[0] == 'null':
                    bc.append(iv(bn, i, j, k, f, 0))
                elif kind[0] == 'wall':
                    bc += [iv(bn, i, j, k, f, 302), rv(kind[1](b.face_centre(f, i, j, k)))]
                elif kind[0] == 'conn':
                    ob, of = kind[1], kind[2]
                    o = blocks[ob - 1]
                    oi, oj, ok = i, j, k
                    if of in (1, 2): oi = 1 if of == 1 else o.nx
                    if of in (3, 4): oj = 1 if of == 3 else o.ny
                    if of in (5, 6): ok = 1 if of == 5 else o.nz
                    bc += [iv(bn, i, j, k, f, 101), iv(ob, oi, oj, ok, of, 1, 0, 0, 1)]
                elif kind[0] == 'chim':
                    ob = kind[1]
                    o = blocks[ob - 1]
                    d1, d2 = (1, 2) if f == 2 else (o.nx, o.nx - 1)     # donor cells across: first and second layer
                    bc += [iv(bn, i, j, k, f, 102), iv(1, 1), iv(ob, d1, j, k) + '    %.10E' % 1.0,
                           iv(ob, d2, j, k) + '    %.10E' % 1.0]
    open(os.path.join(out, 'MESH', 'mesh.tec'), 'w').write(''.join(mesh))
    open(os.path.join(out, 'INPUT', 'ic.tec'), 'w').write(''.join(ic))
    open(os.path.join(out, 'INPUT', 'bc.txt'), 'w').write('\n'.join(bc) + '\n')
    open(os.path.join(out, 'INPUT', 'phase.txt'), 'w').write('solid-bulk\nmat 1\n')
    with open(os.path.join(out, 'INPUT', 'properties.dat'), 'w') as fh:
        fh.write('TITLE = "Mass Thermodynamic Properties"\n'
                 'VARIABLES = "Temperature", "Cp", "Density", "Conductivity", "Energy"\nZONE T="mat"\nI=2500, F=POINT\n')
        fh.writelines('%.1f %.6f %.6f %.6f %.6f\n' % (t, CP, RHO, K, RHO * CP * t) for t in range(1, 2501))
    open(os.path.join(out, 'input.ini'), 'w').write(ini)


def ini_text(iters, steady=True, vnn=0.25, scheme='RK3', extra=''):
    return f"""[FUSS-Parameters]
newrun = true
res-threshold = 1e-30
iter-threshold = {iters}

[FUSS-Numerics]
time-scheme = {scheme}
vnn = {vnn}
time-accurate = false
integration-variables = cons
{extra}
[FUSS-IO]
sol-diter     = {iters}
sol-overwrite = true
sol-format    = tecplot ascii
res-diter = 1
shell-diter = 100000
ini-diter = 100000000
"""


def walls(blocks, which, fn):
    for b in blocks:
        for f in which:
            b.faces.setdefault(f, ('wall', fn))


# --------------------------------------------------------------------------------------------------------------- parts
GZ, A = 400.0, 800.0                     # thin3d: T = A + GZ z


def thin3d():
    n = 4
    b1 = Block('Block1', n, n, 1, lambda i, j, k: (i / n, j / n, 0.2 * k))
    b2 = Block('Block2', n, n, n, lambda i, j, k: (i / n, j / n, 0.2 + 0.8 * k / n))
    tz = lambda c: A + GZ * c[2]
    for b in (b1, b2):
        for f in (1, 2, 3, 4):
            b.faces[f] = ('wall', tz)
    b1.faces[5] = ('wall', tz)
    b1.faces[6] = ('conn', 2, 5)
    b2.faces[5] = ('conn', 1, 6)
    b2.faces[6] = ('wall', tz)
    return [b1, b2], ini_text(40000, vnn=0.5, extra='irs-beta = 0.5\n')


def plane_fn(n, ny, h, y0):
    return lambda i, j, k: (i / n, y0 + h * j / ny, (k - 0.5) * LZ)


def mixed2d():
    n = 4
    b1 = Block('Block1', n, n, 1, lambda i, j, k: (i / n, 0.5 + j / n, (k - 0.5) * LZ))
    phi = math.radians(2.5)
    b2 = Block('Block2', n, n, 1, lambda i, j, k: (1.0 + i / n, (0.5 + j / n) * math.cos(phi),
                                                    (0.5 + j / n) * math.sin(phi) * (1 if k else -1)))
    walls([b1, b2], (1, 2, 3, 4), lambda c: T0)
    return [b1, b2], ini_text(5)


def fold_node(n):
    return lambda i, j, k: ((i + 1.5 if (i, j) == (3, 3) else i) / n, j / n, (k - 0.5) * LZ)


def fold2d():
    n = 6
    b1 = Block('Block1', n, n, 1, fold_node(n))
    walls([b1], (1, 2, 3, 4), lambda c: T0)
    return [b1], ini_text(5)


def fold3d():
    n = 4
    b1 = Block('Block1', n, n, n, lambda i, j, k: ((i + 1.5 if (i, j, k) == (2, 2, 2) else i) / n, j / n, k / n))
    walls([b1], (1, 2, 3, 4, 5, 6), lambda c: T0)
    return [b1], ini_text(5)


def tapered(n, ny):
    """Block 1 of chimturn: x in [0, 1], the height of the last node plane (the chimera face) a quarter of the one before."""
    return lambda i, j, k: (i / n, (j / ny) * (0.25 if i == n else 1.0), (k - 0.5) * LZ)


def chimturn():
    n, ny = 8, 6
    b1 = Block('Block1', n, ny, 1, tapered(n, ny))
    b2 = Block('Block2', n, ny, 1, lambda i, j, k: (1.0 + i / n, 0.25 * j / ny, (k - 0.5) * LZ))
    walls([b1, b2], (3, 4), lambda c: T0)
    b1.faces[1] = ('wall', lambda c: T0)
    b1.faces[2] = ('chim', 2)
    b2.faces[1] = ('chim', 1)
    b2.faces[2] = ('wall', lambda c: T0)
    return [b1, b2], ini_text(5)


def chimwall():
    n, ny = 8, 6
    b1 = Block('Block1', n, ny, 1, tapered(n, ny))
    walls([b1], (1, 2, 3, 4), lambda c: T0)
    return [b1], ini_text(5)


PARTS = {'thin3d': thin3d, 'mixed2d': mixed2d, 'fold2d': fold2d, 'fold3d': fold3d, 'chimturn': chimturn, 'chimwall': chimwall}

if __name__ == '__main__':
    if len(sys.argv) != 3 or sys.argv[1] not in PARTS:
        sys.exit(__doc__)
    blocks, ini = PARTS[sys.argv[1]]()
    write_case(sys.argv[2], blocks, ini)
    print('%s: %d block(s)' % (sys.argv[1], len(blocks)))
