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
    else:
        sys.exit('unknown family ' + family)
    return blocks, g


def write_case(family, out, metric=None):
    blocks, g = build(family)
    for sub in ('MESH', 'INPUT', 'OUTPUT'):
        os.makedirs(os.path.join(out, sub), exist_ok=True)
    mesh, ic, bc = [' VARIABLES ="x" "y" "z"\n'], [' VARIABLES ="x" "y" "z" T matID\n'], []
    for bn, b in enumerate(blocks, start=1):
        kk = b.nz + 1
        zone = (' ZONE  T = %s, I=%d, J=%d, K=%d, DATAPACKING=BLOCK, VARLOCATION=([1-3]=NODAL), SOLUTIONTIME=0.0, STRANDID = 0\n'
                % (b.name, b.nx + 1, b.ny + 1, kk))
        flat = b.nodes_flat()
        mesh += [zone, '\n'.join(e(v) for v in flat) + '\n']
        tcell = [exact(b.cell_centre(i, j, k), g) for k in range(1, b.nz + 1) for j in range(1, b.ny + 1)
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
                    oi = 1 if of == 1 else blocks[ob - 1].nx
                    bc += [iv(bn, i, j, k, f, 101), iv(ob, oi, j, k, of, 1, 0, 0, 1)]
                else:
                    t = kind[1]
                    c, nrm, _ = b.face_centre_normal(f, i, j, k)
                    tw = exact(c, g)
                    q = K * sum(g[d] * nrm[d] for d in range(3))     # flux into the solid
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
    meta = {'family': family, 'g': g, 'k': K, 'blocks': [{'name': b.name, 'n': [b.nx, b.ny, b.nz],
             'faces': {str(f): list(b.faces[f]) for f in range(1, 7)}} for b in blocks]}
    json.dump(meta, open(os.path.join(out, 'patch.json'), 'w'), indent=1)
    print('%s: %d block(s), %d boundary records' % (family, len(blocks), sum(1 for l in bc if len(l.split()) == 6)))


def check_case(case, tol):
    meta = json.load(open(os.path.join(case, 'patch.json')))
    blocks, g = build(meta['family'])
    gnorm = math.sqrt(sum(v * v for v in g))
    diag = read_tec(os.path.join(case, 'OUTPUT', 'diagnostic.tec'))
    groups = {}
    status = 0
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
