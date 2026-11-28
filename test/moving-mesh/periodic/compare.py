#!/usr/bin/env python3
"""compare.py <ring dir> <sector dir> <tol> [--vars T,p_pg] [--nodes] [--differ]

Plan 14 S9d, Part B: a periodic sector (generate.py) against every sector of the ring it stands for, from the two
OUTPUT/field.tec. For each named cell variable the largest |v_sector - v_ring| / |v_ring| over the cells of the 6 sectors
of the ring (default: T); with --nodes also the largest distance of a node of the sector, turned by s x 60 deg about x
(in the sense of the generator's cart), from the node of sector s of the ring, divided by the outer radius (a moving
mesh: the nodes are part of the answer).
Exit 0 when every number is <= tol; with --differ, exit 0 when the largest is > tol (a control that must not be the ring).
Python 3.6, no numpy (upstream's fusstec.read_tec).
"""
import math
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'test', 'numerics-features', 'patch_test'))
sys.path.insert(0, os.path.join(ROOT, 'test', 'numerics-features', 'periodic_sector'))
from fusstec import read_tec  # noqa: E402
from make_sector import NS, NX, NT, NR, NTA, DTH, R1, cart  # noqa: E402

# the sense in which sector s of the ring is turned from sector 0 about x, from the generator itself (plan 14 S9d
# D-S9d-6 mirrors upstream's mesh, z -> -z: a right-handed block, sector s at -s DTH; upstream's at +s DTH)
SENSE = 1.0 if cart((0.0, 1.0, 0.0), 1)[2] > 0.0 else -1.0


def cells(zones, var):
    """var of one sector in the order (k, j, i), j over the NT cells in theta (blocks A, B or the single block)"""
    if len(zones) == 1:
        z = zones[0][var]
        return [z[(k * NT + j) * NX + i] for k in range(NR) for j in range(NT) for i in range(NX)]
    a, b = zones[0][var], zones[1][var]
    return [a[(k * NTA + j) * NX + i] if j < NTA else b[(k * (NT - NTA) + j - NTA) * NX + i]
            for k in range(NR) for j in range(NT) for i in range(NX)]


def turn(p, s):
    c, w = math.cos(s * DTH), SENSE * math.sin(s * DTH)
    return (p[0], c * p[1] - w * p[2], w * p[1] + c * p[2])


args = sys.argv[1:]
differ = '--differ' in args
nodes = '--nodes' in args
names = ['T']
if '--vars' in args:
    names = args[args.index('--vars') + 1].split(',')
ring = list(read_tec(os.path.join(args[0], 'OUTPUT', 'field.tec')).values())
sec = list(read_tec(os.path.join(args[1], 'OUTPUT', 'field.tec')).values())
tol = float(args[2])
nb = len(ring) // NS
if nb != len(sec) or nb * NS != len(ring):
    sys.exit('compare.py: %d ring zones and %d sector zones do not make %d sectors' % (len(ring), len(sec), NS))
worst = {}
for v in names:
    vs = cells(sec, v)
    w = 0.0
    for s in range(NS):
        vr = cells(ring[nb * s:nb * s + nb], v)
        w = max(w, max(abs(a - b) / max(abs(b), 1e-300) for a, b in zip(vs, vr)))
    worst[v] = w
if nodes:
    w = 0.0
    for s in range(NS):
        for z_s, z_r in zip(sec, ring[nb * s:nb * s + nb]):
            for p in zip(z_s['x'], z_s['y'], z_s['z'], z_r['x'], z_r['y'], z_r['z']):
                q = turn(p[0:3], s)
                w = max(w, math.sqrt(sum((q[c] - p[3 + c]) ** 2 for c in range(3))) / R1)
    worst['nodes'] = w
top = max(worst.values())
print('%s vs %s: %s; tolerance %.0e%s' % (os.path.basename(args[1].rstrip('/')), os.path.basename(args[0].rstrip('/')),
      ', '.join('%s %.3e' % kv for kv in worst.items()), tol, ', must exceed it' if differ else ''))
sys.exit(0 if ((top > tol) if differ else (top <= tol)) else 1)
