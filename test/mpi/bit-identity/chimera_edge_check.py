#!/usr/bin/env python3
"""ChimeraEdge: a linear field is a solution of the conduction equation, so the fluxes of every cell sum to zero. On the
sheared chimera pair of make_case.py (edge) with walls consistent with the field, the state after a few iterations is
the initial one up to round-off, also at the cells at the ends of the chimera faces (edges of the blocks), where the ghost
row beyond the edge is extrapolated along the row (Fill_BC_Ghost_Chimera) and not copied.

Usage: chimera_edge_check.py <case dir> [tolerance in K]    reads <case dir>/OUTPUT/field.tec and INPUT/ic.tec
Prints the largest |T - T0| of each block, with its cell, and exits 1 above the tolerance.
"""
import re
import sys

import numpy as np


def read_tec(path):
    """Zones of a Tecplot ASCII BLOCK file: {title: {variable: values}}."""
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


def main():
    case = sys.argv[1]
    tol = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-6
    ic = read_tec(f'{case}/INPUT/ic.tec')
    out = read_tec(f'{case}/OUTPUT/field.tec')
    worst = 0.0
    for (t0name, z0), (tname, z1) in zip(sorted(ic.items()), sorted(out.items())):
        ni, nj, nk = (n - 1 for n in z0['dims'])
        d = np.abs(np.array(z1['T']) - np.array(z0['T'])).reshape(nk, nj, ni)
        k, j, i = np.unravel_index(np.argmax(d), d.shape)
        print(f'  {tname}: max |T - T0| = {d.max():.3e} K at cell ({i + 1},{j + 1},{k + 1})')
        worst = max(worst, d.max())
    print(f'ChimeraEdge: max |T - T0| = {worst:.3e} K (tolerance {tol:g})')
    sys.exit(0 if worst <= tol else 1)


if __name__ == '__main__':
    main()
