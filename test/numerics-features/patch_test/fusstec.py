"""Tecplot ASCII (BLOCK, nodal x y z + cell-centred variables) reader and small helpers shared by the FUSS stress
battery (patch_case.py, restart_cmp.py, balance.py). Python 3.6 compatible, no numpy needed.

read_tec(path) -> ordered dict  zone title -> {'dims': (I, J, K) nodes, 'x','y','z': node lists, <var>: cell lists}
cell_centres(zone) -> list of (xc, yc, zc) in cell order (i fastest, then j, then k)
"""
import re
from collections import OrderedDict


def read_tec(path):
    text = open(path).read().split('\n')
    names, zones, i = [], OrderedDict(), 0
    while i < len(text):
        line = text[i].strip()
        if re.match(r'(?i)variables', line):
            names = [a or b for a, b in re.findall(r'"([^"]*)"|([^\s,"]+)', line.split('=', 1)[1])]
        elif re.match(r'(?i)zone', line):
            title = re.search(r'T\s*=\s*"?([^",\s]+)', line).group(1)
            ni, nj, nk = (int(re.search(r'\b%s\s*=\s*(\d+)' % c, line).group(1)) for c in 'IJK')
            cc = re.search(r'\[(\d+)-(\d+)\]\s*=\s*CELLCENTERED', line, re.I)
            sizes = [(max(ni - 1, 1) * max(nj - 1, 1) * max(nk - 1, 1))
                     if cc and int(cc.group(1)) <= v <= int(cc.group(2)) else ni * nj * nk
                     for v in range(1, len(names) + 1)]
            vals, i = [], i + 1
            while len(vals) < sum(sizes) and i < len(text):
                s = text[i].strip()
                if s and not s[0].isalpha():
                    vals += [float(x) for x in s.split()]
                elif re.match(r'(?i)(zone|dataset|text)', s):
                    break
                i += 1
            zone, off = OrderedDict(), 0
            for n, s in zip(names, sizes):
                zone[n] = vals[off:off + s]
                off += s
            zone['dims'] = (ni, nj, nk)
            zones[title] = zone
            continue
        i += 1
    return zones


def cell_centres(zone):
    I, J, K = zone['dims']
    xs, ys, zs = zone['x'], zone['y'], zone['z']
    out = []
    for k in range(max(K - 1, 1)):
        for j in range(max(J - 1, 1)):
            for i in range(max(I - 1, 1)):
                c = [0.0, 0.0, 0.0]
                for dk in (0, 1) if K > 1 else (0,):
                    for dj in (0, 1) if J > 1 else (0,):
                        for di in (0, 1) if I > 1 else (0,):
                            n = (k + dk) * I * J + (j + dj) * I + (i + di)
                            c[0] += xs[n]; c[1] += ys[n]; c[2] += zs[n]
                m = (2 if K > 1 else 1) * (2 if J > 1 else 1) * (2 if I > 1 else 1)
                out.append((c[0] / m, c[1] / m, c[2] / m))
    return out


def max_abs_diff(a, b):
    return max((abs(p - q) for p, q in zip(a, b)), default=0.0)
