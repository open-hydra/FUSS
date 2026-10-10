"""Steady 1-D conduction in a rod with a convective wall that also receives a prescribed flux (BC 303).

x = 0: T = T0 (BC 302). x = L: q = hconv (Tref - Tw) + qrad into the rod (BC 303, bc.txt: hconv, qrad, Tref).
Constant k: T is linear, T(x) = T0 + (TR - T0) x / L, and the flux k (TR - T0) / L through the rod is the
one entering at x = L, which gives TR = (k T0 / L + hconv Tref + qrad) / (k / L + hconv). The finite-volume
scheme is exact for a linear profile, so the converged field matches it to the convergence error.
Without qrad (the 303 flux used to pass a never-set component) TR would be 347.4 K instead of 366.3 K.
Reads OUTPUT/field.tec and OUTPUT/wall.tec: run ./FUSS.sh solve first.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
L, K = 0.1, 55.563                               # rod length [m], conductivity [W/(m K)] (properties.dat)
T0, HCONV, TREF, QRAD = 300.0, 500.0, 400.0, 2.0e4
TOL_T = 1e-3                                     # K
TOL_Q = 1e-4                                     # relative to the flux through the rod


def exact_TR(qrad):
    return (K * T0 / L + HCONV * TREF + qrad) / (K / L + HCONV)


def read_tec(path):
    """Zones of a Tecplot ASCII BLOCK file: {title: {variable: values}} (names quoted or bare)."""
    text = open(path).read().split('\n')
    names, zones, i = [], {}, 0
    while i < len(text):
        line = text[i].strip()
        if re.match(r'(?i)variables', line):
            names = [a or b for a, b in re.findall(r'"([^"]*)"|([^\s,"]+)', line.split('=', 1)[1])]
        elif re.match(r'(?i)zone', line):
            title = re.search(r'T\s*=\s*"?([^",\s]+)', line).group(1)
            ni, nj, nk = (int(re.search(r'\b%s\s*=\s*(\d+)' % c, line).group(1)) for c in 'IJK')
            cc = re.search(r'\[(\d+)-(\d+)\]\s*=\s*CELLCENTERED', line, re.I)
            sizes = [(max(ni - 1, 1) * max(nj - 1, 1) * max(nk - 1, 1)) if cc and int(cc.group(1)) <= v <= int(cc.group(2))
                     else ni * nj * nk for v in range(1, len(names) + 1)]
            vals, i = [], i + 1
            while len(vals) < sum(sizes):
                vals += [float(x) for x in text[i].split()]
                i += 1
            zone, off = {}, 0
            for n, s in zip(names, sizes):
                zone[n] = vals[off:off + s]
                off += s
            zones[title] = zone
            continue
        i += 1
    return zones


def main():
    field, wall = os.path.join(HERE, 'OUTPUT', 'field.tec'), os.path.join(HERE, 'OUTPUT', 'wall.tec')
    for f in (field, wall):
        if not os.path.isfile(f):
            print(f'ERROR: {f} not found: run ./FUSS.sh solve first')
            sys.exit(1)
    TR = exact_TR(QRAD)
    q = K * (TR - T0) / L

    # Field: cell centres from the nodes of the first cell row (x is the fastest index)
    blk = next(iter(read_tec(field).values()))
    T = blk['T']
    nx = len(T) // 2                                 # 2 cells across the rod, 1 in z
    x_nodes = blk['x'][:nx + 1]
    err_T = max(abs(T[j * nx + i] - (T0 + (TR - T0) * 0.5 * (x_nodes[i] + x_nodes[i + 1]) / L))
                for j in range(2) for i in range(nx))

    # Wall: B1F2 is the BC 303 face, B1F1 the 300 K wall; qw is the flux into the rod, so -q on B1F1
    walls = read_tec(wall)
    Tw, qw = walls['B1F2']['Tw'], walls['B1F2']['qw']
    err_Tw = max(abs(t - TR) for t in Tw)
    err_qw = max(abs(v - q) for v in qw) / abs(q)
    err_q1 = max(abs(v + q) for v in walls['B1F1']['qw']) / abs(q)

    print('-' * 72)
    print(' 1-D rod, convective wall with a prescribed flux (BC 303)')
    print('-' * 72)
    print(f'  exact TR {TR:.4f} K (without qrad {exact_TR(0.0):.4f} K), flux through the rod {q:.2f} W/m2')
    print(f'  field   max |T - T_exact|          {err_T:.3e} K    (tol {TOL_T:g} K)')
    print(f'  B1F2    max |Tw - TR|              {err_Tw:.3e} K    (tol {TOL_T:g} K)')
    print(f'  B1F2    max |qw - q| / q           {err_qw:.3e}      (tol {TOL_Q:g})')
    print(f'  B1F1    max |qw + q| / q           {err_q1:.3e}      (tol {TOL_Q:g})')
    ok = err_T < TOL_T and err_Tw < TOL_T and err_qw < TOL_Q and err_q1 < TOL_Q
    print('  Result:', 'PASS' if ok else 'FAIL')
    print('-' * 72)
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
