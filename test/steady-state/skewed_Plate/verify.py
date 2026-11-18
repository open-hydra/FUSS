"""Linear patch test on a 45 deg sheared plate with the four wall types 301, 302, 303, 305 (see make_inputs.py).

A linear field is reproduced exactly by the finite-volume scheme on any affine mesh, provided the wall
conditions are consistent with the interior scheme: the converged solution must match T = A + B x + C y
to the convergence error, in the field and in the wall output (Tw at the face centres, qw = heat flux
into the solid). Reads OUTPUT/field.tec and OUTPUT/wall.tec: run ./FUSS.sh solve.
"""
import os
import re
import sys

import numpy as np

from make_inputs import N, Q_IN, exact

HERE = os.path.dirname(os.path.abspath(__file__))
TOL_T = 1e-6                                     # K
TOL_Q = 1e-6                                     # relative to the flux


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
            zone['dims'] = (ni, nj, nk)
            zones[title] = zone
            continue
        i += 1
    return zones



def centres(zone):
    """Cell centres (x, y) of a zone: node coordinates averaged over each direction with more than one node."""
    ni, nj, nk = zone['dims']
    out = []
    for c in ('x', 'y'):
        a = np.array(zone[c]).reshape(nk, nj, ni)
        for ax, n in ((2, ni), (1, nj), (0, nk)):
            if n > 1:
                a = 0.5 * (a.take(range(n - 1), axis=ax) + a.take(range(1, n), axis=ax))
        out.append(a.ravel())
    return out


def main():
    field, wall = os.path.join(HERE, 'OUTPUT', 'field.tec'), os.path.join(HERE, 'OUTPUT', 'wall.tec')
    for f in (field, wall):
        if not os.path.isfile(f):
            print(f'ERROR: {f} not found: run ./FUSS.sh solve first')
            sys.exit(1)
    print('-' * 72)
    print(f' Linear patch test, {N} x {N} cells sheared by 45 deg, walls 302 / 303 / 301 / 305')
    print('-' * 72)
    blk = next(iter(read_tec(field).values()))
    x, y = centres(blk)
    err = np.max(np.abs(np.array(blk['T']) - exact(x, y)))
    print(f'  field   max |T - T_exact|      {err:.3e} K    (tol {TOL_T:g} K)')
    ok = err < TOL_T
    walls = read_tec(wall)
    for f, bc in ((1, 302), (2, 303), (3, 301), (4, 305)):
        z = walls[f'B1F{f}']
        x, y = centres(z)
        eT = np.max(np.abs(np.array(z['Tw']) - exact(x, y)))
        eq = np.max(np.abs(np.array(z['qw']) - Q_IN[f])) / abs(Q_IN[f])
        print(f'  B1F{f} ({bc})  max |Tw - T_exact| {eT:.3e} K   max |qw - q| / |q| {eq:.3e}   (q = {Q_IN[f]:.1f} W/m2)')
        ok = ok and eT < TOL_T and eq < TOL_Q
    print('  Result:', 'PASS' if ok else 'FAIL')
    print('-' * 72)
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
