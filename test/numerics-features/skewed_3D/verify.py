"""Skewed3D: the time step of every cell after one iteration against a replica of the formula
    dt = vnn rho cp / (f sum_d k_d |grad xi_d|^2),  f = (2 + m + 1/m)/4 on cells whose three directions are skewed in all
    three pairs (cos12^2 + cos13^2 + cos23^2 + 2 cos12 cos13 cos23 > 1; m + 1 the largest eigenvalue of the matrix of the
    cosines between the rows grad xi_d of the metric), 1 otherwise;
k_d the larger conductivity of the two faces of the cell in direction d: its own k on a wall face and between cells of one
material, the series (harmonic) mean 2 k1 k2/(k1 + k2) at the interface of the two materials (the cells are equal, so the
half widths are). The replica builds the metric from the nodes of MESH/mesh.tec as the solver does (rows of M = rows of the
inverse transpose of the mean edges) and reads dt from OUTPUT/diagnostic.tec. Standard library only.
Tolerance 1e-12 relative. The factor of this mesh is 1.0368 (red on a solver without it: 3.7 % off, and a factor 2 off when
vnn follows the convention of FUSS before October 2026)."""
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_inputs as mk  # noqa: E402


def read_block(path, key):
    """The numbers of the first zone of a Tecplot ASCII BLOCK file, as {variable: list}, and the dims."""
    text = open(path).read().split('\n')
    names, vals, dims = [], [], None
    for i, line in enumerate(text):
        s = line.strip()
        if re.match(r'(?i)variables', s):
            names = re.findall(r'"([^"]*)"|([^\s,"]+)', s.split('=', 1)[1])
            names = [a or b for a, b in names]
        elif re.match(r'(?i)zone', s):
            dims = tuple(int(re.search(r'\b%s\s*=\s*(\d+)' % c, s).group(1)) for c in 'IJK')
            for t in text[i + 1:]:
                t = t.strip()
                if t and not t[0].isalpha():
                    vals += [float(x) for x in t.split()]
            break
    return names, vals, dims


def inv3(a):
    c = [[a[1][1] * a[2][2] - a[1][2] * a[2][1], -(a[1][0] * a[2][2] - a[1][2] * a[2][0]), a[1][0] * a[2][1] - a[1][1] * a[2][0]],
         [-(a[0][1] * a[2][2] - a[0][2] * a[2][1]), a[0][0] * a[2][2] - a[0][2] * a[2][0], -(a[0][0] * a[2][1] - a[0][1] * a[2][0])],
         [a[0][1] * a[1][2] - a[0][2] * a[1][1], -(a[0][0] * a[1][2] - a[0][2] * a[1][0]), a[0][0] * a[1][1] - a[0][1] * a[1][0]]]
    det = sum(a[0][j] * c[0][j] for j in range(3))
    return [[c[i][j] / det for j in range(3)] for i in range(3)]      # (cofactor/det) = inverse transpose of a


def skew_factor(m):
    n = [math.sqrt(sum(v * v for v in r)) for r in m]
    dot = lambda p, q: sum(m[p][d] * m[q][d] for d in range(3)) / (n[p] * n[q])
    c12, c13, c23 = dot(0, 1), dot(0, 2), dot(1, 2)
    s, q = c12 ** 2 + c13 ** 2 + c23 ** 2, c12 * c13 * c23
    if s + 2 * q <= 1:
        return 1.0
    f = lambda x: x ** 3 - s * x - 2 * q          # m = largest root, eigenvalue of the cosine matrix minus 1
    lo, hi = 0.0, 2.0
    for _ in range(200):                           # f(0) = -2q < 0 for q > 0, f(2) > 0: bisection on the largest root
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if f(mid) < 0 else (lo, mid)
    mm = 0.5 * (lo + hi)
    return 0.25 * (2 + mm + 1 / mm)


def main():
    names, vals, dims = read_block(os.path.join(HERE, 'MESH', 'mesh.tec'), 'x')
    ni, nj, nk = dims
    nn = ni * nj * nk
    node = lambda i, j, k: tuple(vals[d * nn + (k * nj + j) * ni + i] for d in range(3))
    dnames, dvals, ddims = read_block(os.path.join(HERE, 'OUTPUT', 'diagnostic.tec'), 'dt')
    nx, ny, nz = ni - 1, nj - 1, nk - 1
    ncell = nx * ny * nz
    off = 3 * nn + dnames.index('dt') * 0
    # variables: x y z (nodal) then E dt (cell-centred)
    first_cell = 3 * nn
    dt = dvals[first_cell + ncell: first_cell + 2 * ncell] if dnames.index('dt') == 4 else dvals[first_cell + (dnames.index('dt') - 3) * ncell: first_cell + (dnames.index('dt') - 2) * ncell]
    worst, worst_plain, fmax = 0.0, 0.0, 1.0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                p = lambda a, b, c: node(i + a, j + b, k + c)
                a = [[0.25 * sum(p(*q)[d] for q in qs) for d in range(3)] for qs in ()]
                A = []
                for axis in range(3):
                    row = [0.0, 0.0, 0.0]
                    for b in (0, 1):
                        for c in (0, 1):
                            lo = [0, 0, 0]; hi = [0, 0, 0]
                            other = [x for x in range(3) if x != axis]
                            lo[other[0]] = b; hi[other[0]] = b; lo[other[1]] = c; hi[other[1]] = c
                            hi[axis] = 1
                            pl, ph = p(*lo), p(*hi)
                            for d in range(3):
                                row[d] += ph[d] - pl[d]
                    A.append([0.25 * v for v in row])
                M = inv3(A)
                G = [sum(v * v for v in M[d]) for d in range(3)]
                f = skew_factor(M)
                fmax = max(fmax, f)
                mat1 = i < mk.NX // 2
                kc = mk.K1 if mat1 else mk.K2
                kd = [kc, kc, kc]
                ks = 2 * mk.K1 * mk.K2 / (mk.K1 + mk.K2)
                if i == mk.NX // 2 - 1: kd[0] = max(kc, ks)
                if i == mk.NX // 2: kd[0] = max(kc, ks)
                ref = mk.VNN * mk.RHO * mk.CP / (f * sum(kd[d] * G[d] for d in range(3)))
                plain = mk.VNN * mk.RHO * mk.CP / sum(kd[d] * G[d] for d in range(3))
                got = dt[(k * ny + j) * nx + i]
                worst = max(worst, abs(got - ref) / ref)
                worst_plain = max(worst_plain, abs(got - plain) / plain)
    print('Skewed3D: %d cells, skew factor of the mesh %.6f, max |dt - replica| / replica = %.3e (without the factor: %.3e)'
          % (ncell, fmax, worst, worst_plain))
    if worst <= 1e-12:
        print('Skewed3D: PASS')
        return 0
    print('Skewed3D: FAIL (the time step is not vnn rho cp / (f sum_d k_d |grad xi_d|^2))')
    return 1


if __name__ == '__main__':
    sys.exit(main())
