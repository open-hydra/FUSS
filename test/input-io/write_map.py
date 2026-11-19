#!/usr/bin/env python3
# write_map.py <ic.tec> <map file> <rank block 1> <rank block 2> ...: a decomposition map (ATLAS mesh decomposer format) for
# the blocks of a Tecplot ASCII file: new parent i0 i1 j0 j1 k0 k1 ni nj nk cells rank, with the sizes read from the zone headers
import re, sys
ic, out, ranks = sys.argv[1], sys.argv[2], [int(r) for r in sys.argv[3:]]
dims = []
for line in open(ic):
    m = re.search(r'I=(\d+),\s*J=(\d+),\s*K=(\d+)', line)
    if m:
        dims.append(tuple(max(int(x) - 1, 1) for x in m.groups()))
with open(out, 'w') as f:
    f.write('# decomposition map written by write_map.py\n')
    for b, (ni, nj, nk) in enumerate(dims, 1):
        f.write('%d 0 1 %d 1 %d 1 %d %d %d %d %d %d\n' % (b, ni, nj, nk, ni, nj, nk, ni * nj * nk, ranks[b - 1]))
    f.write('# end of block lines\n')
