#!/usr/bin/env python3
"""compare_split.py <one dir> <two dir> [tol]: max |T_two - T_one| / T of the cells of the two-block run (block 2 offset by
8 cells in i) against the one-block run, from OUTPUT/field.tec; exit 1 above tol (default 1e-12)."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'patch_test'))
from fusstec import read_tec  # noqa: E402

NX, CUT = 16, 8
one = list(read_tec(os.path.join(sys.argv[1], 'OUTPUT', 'field.tec')).values())[0]['T']
two = [z['T'] for z in read_tec(os.path.join(sys.argv[2], 'OUTPUT', 'field.tec')).values()]
worst = 0.0
for j in range(NX):
    for i in range(NX):
        t2 = two[0][j * CUT + i] if i < CUT else two[1][j * (NX - CUT) + i - CUT]
        worst = max(worst, abs(t2 - one[j * NX + i]) / abs(one[j * NX + i]))
tol = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-12
print('max relative difference two blocks vs one block: %.3e (tolerance %.0e)' % (worst, tol))
sys.exit(0 if worst <= tol else 1)
