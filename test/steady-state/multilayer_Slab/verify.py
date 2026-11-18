"""Three material layers in series inside one block (see make_inputs.py).

With k constant in each material, the steady field is linear in each layer and the heat flux is the
same through all of them. A face between two materials takes the series conductivity of its two half
cells, so the finite-volume scheme reproduces that field exactly, whatever the cell widths on the two
sides: the converged solution must match it to the convergence error, and the wall heat flux qw (into
the solid) must be +q at x = 0 and -q at x = L. Reads OUTPUT/field.tec and OUTPUT/wall.tec: run
./FUSS.sh solve.
"""
import os
import re
import sys

import numpy as np

from make_inputs import LAYERS, Q, exact

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


def centres_x(zone):
    """x of the cell centres of a zone: node x averaged over each direction with more than one node."""
    ni, nj, nk = zone['dims']
    a = np.array(zone['x']).reshape(nk, nj, ni)
    for ax, n in ((2, ni), (1, nj), (0, nk)):
        if n > 1:
            a = 0.5 * (a.take(range(n - 1), axis=ax) + a.take(range(1, n), axis=ax))
    return a.ravel()


def main():
    field, wall = os.path.join(HERE, 'OUTPUT', 'field.tec'), os.path.join(HERE, 'OUTPUT', 'wall.tec')
    for f in (field, wall):
        if not os.path.isfile(f):
            print(f'ERROR: {f} not found: run ./FUSS.sh solve first')
            sys.exit(1)
    print('-' * 72)
    print(' Three layers in series inside one block: k = '
          + ' / '.join(f'{k:g}' for _, _, _, k, _, _ in LAYERS) + f' W/(m K), q = {Q:.1f} W/m2')
    print('-' * 72)
    blk = next(iter(read_tec(field).values()))
    x = centres_x(blk)
    err = np.max(np.abs(np.array(blk['T']) - np.array([exact(v) for v in x])))
    print(f'  field   max |T - T_exact|      {err:.3e} K    (tol {TOL_T:g} K)')
    ok = err < TOL_T
    walls = read_tec(wall)
    for f, q in ((1, Q), (2, -Q)):
        eq = np.max(np.abs(np.array(walls[f'B1F{f}']['qw']) - q)) / abs(Q)
        print(f'  B1F{f} (302)  max |qw - q| / |q| {eq:.3e}   (q = {q:.1f} W/m2)')
        ok = ok and eq < TOL_Q
    print('  Result:', 'PASS' if ok else 'FAIL')
    print('-' * 72)
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
