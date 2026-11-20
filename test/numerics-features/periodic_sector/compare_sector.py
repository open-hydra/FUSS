#!/usr/bin/env python3
"""compare_sector.py <ring dir> <sector dir> <tol> [--differ]

max |T_sector - T_ring| / T over the cells of every sector of the ring run (ring: blocks 2s+1, 2s+2; ring1: block s+1)
against the sector run (blocks A and B, or the single block of sector1: its first NTA cells in theta are A), from
OUTPUT/field.tec. Exit 1 when the difference exceeds tol; with --differ, exit 1 when it does not (the symmetric control
must not be the ring).
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'patch_test'))
sys.path.insert(0, HERE)
from fusstec import read_tec  # noqa: E402
from make_sector import NS, NX, NT, NR, NTA  # noqa: E402


def sector_cells(zones):
    """T of one sector in the order (k, j, i), j over the NT cells in theta"""
    if len(zones) == 1:
        z = zones[0]
        return [z[(k * NT + j) * NX + i] for k in range(NR) for j in range(NT) for i in range(NX)]
    a, b = zones
    out = []
    for k in range(NR):
        for j in range(NT):
            for i in range(NX):
                if j < NTA:
                    out.append(a[(k * NTA + j) * NX + i])
                else:
                    out.append(b[(k * (NT - NTA) + j - NTA) * NX + i])
    return out


ring = [z['T'] for z in read_tec(os.path.join(sys.argv[1], 'OUTPUT', 'field.tec')).values()]
sec = sector_cells([z['T'] for z in read_tec(os.path.join(sys.argv[2], 'OUTPUT', 'field.tec')).values()])
tol = float(sys.argv[3])
differ = '--differ' in sys.argv
nb = len(ring) // NS          # blocks per sector of the ring
worst, spread = 0.0, 0.0
r0 = sector_cells(ring[0:nb])
for s in range(NS):
    rs = sector_cells(ring[nb * s:nb * s + nb])
    worst = max(worst, max(abs(a - b) / abs(b) for a, b in zip(sec, rs)))
    spread = max(spread, max(abs(a - b) / abs(b) for a, b in zip(rs, r0)))
print('%s vs ring: max relative difference %.3e over the %d sectors (ring sectors among themselves %.3e), '
      'tolerance %.0e%s' % (os.path.basename(sys.argv[2].rstrip('/')), worst, NS, spread, tol,
                            ', must exceed it' if differ else ''))
ok = (worst > tol) if differ else (worst <= tol)
sys.exit(0 if ok else 1)
