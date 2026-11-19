#!/usr/bin/env python3
"""ConnectionSplit: the same mesh as one block and as two blocks joined by a block connection (101).

  make_split.py case <one|two|two-reversed> <case dir>

The mesh is the `rand` family of the patch test (16x16 cells, nodes moved at random, one cell thick); the two-block cases
cut it between cell 8 and 9 and share the node bits. All walls 302 at 1000 K, the initial field a smooth non-uniform
function of the cell centre, RK3, 50 iterations. `two-reversed` lists the records of every block face in the opposite
order. With one flux per connection face and the stencil metric, the two-block run must give the fields of the one-block
run (to round-off: the ghost metric of the connection is built from the source nodes) and the two orders of the records
the same bits (run_split.sh).
"""
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'patch_test'))
import patch_case as pc  # noqa: E402

NX = 16
CUT = 8


def field(c):
    return 1000.0 + 150.0 * math.sin(3.0 * c[0]) * math.cos(2.0 * c[1]) + 40.0 * c[1]


def blocks_of(kind):
    full = pc.build('rand')[0][0]
    if kind == 'one':
        b = full
        for f in (1, 2, 3, 4):
            b.faces[f] = ('wall', 302)
        b.faces[5] = b.faces[6] = ('null',)
        return [b]
    b1 = pc.Block('Block1', CUT, NX, 1, lambda i, j, k: full.nodes[k][j][i])
    b2 = pc.Block('Block2', NX - CUT, NX, 1, lambda i, j, k: full.nodes[k][j][i + CUT])
    for b in (b1, b2):
        b.faces.update({3: ('wall', 302), 4: ('wall', 302), 5: ('null',), 6: ('null',)})
    b1.faces[1] = ('wall', 302); b1.faces[2] = ('conn', 2, 1)
    b2.faces[1] = ('conn', 1, 2); b2.faces[2] = ('wall', 302)
    return [b1, b2]


def write_case(kind, out):
    blocks = blocks_of(kind)
    for sub in ('MESH', 'INPUT', 'OUTPUT'):
        os.makedirs(os.path.join(out, sub), exist_ok=True)
    mesh, ic, bc = [' VARIABLES ="x" "y" "z"\n'], [' VARIABLES ="x" "y" "z" T matID\n'], []
    for bn, b in enumerate(blocks, start=1):
        zone = (' ZONE  T = %s, I=%d, J=%d, K=%d, DATAPACKING=BLOCK, VARLOCATION=([1-3]=NODAL), SOLUTIONTIME=0.0, STRANDID = 0\n'
                % (b.name, b.nx + 1, b.ny + 1, b.nz + 1))
        flat = b.nodes_flat()
        mesh += [zone, '\n'.join(pc.e(v) for v in flat) + '\n']
        tcell = [field(b.cell_centre(i, j, k)) for k in range(1, b.nz + 1) for j in range(1, b.ny + 1)
                 for i in range(1, b.nx + 1)]
        ic += [zone.replace(b.name, 'B%d-SP' % bn).replace('VARLOCATION=([1-3]=NODAL)',
                                                           'VARLOCATION=([1-3]=NODAL,[4-5]=CELLCENTERED)'),
               '\n'.join(pc.e(v) for v in flat + tcell + [1.0] * len(tcell)) + '\n']
        for f in range(1, 7):
            kindf = b.faces[f]
            cells = b.face_cells(f)
            if kind == 'two-reversed':
                cells = cells[::-1]
            for (i, j, k) in cells:
                if kindf[0] == 'null':
                    bc.append(pc.iv(bn, i, j, k, f, 0))
                elif kindf[0] == 'conn':
                    ob, of = kindf[1], kindf[2]
                    oi = 1 if of == 1 else blocks[ob - 1].nx
                    bc += [pc.iv(bn, i, j, k, f, 101), pc.iv(ob, oi, j, k, of, 1, 0, 0, 1)]
                else:
                    bc += [pc.iv(bn, i, j, k, f, 302), pc.rv(1000.0)]
    open(os.path.join(out, 'MESH', 'mesh.tec'), 'w').write(''.join(mesh))
    open(os.path.join(out, 'INPUT', 'ic.tec'), 'w').write(''.join(ic))
    open(os.path.join(out, 'INPUT', 'bc.txt'), 'w').write('\n'.join(bc) + '\n')
    open(os.path.join(out, 'INPUT', 'phase.txt'), 'w').write('solid-bulk\nmat 1\n')
    with open(os.path.join(out, 'INPUT', 'properties.dat'), 'w') as fh:
        fh.write('TITLE = "Mass Thermodynamic Properties"\n'
                 'VARIABLES = "Temperature", "Cp", "Density", "Conductivity", "Energy"\nZONE T="mat"\nI=2500, F=POINT\n')
        fh.writelines('%.1f %.6f %.6f %.6f %.6f\n' % (t, pc.CP, pc.RHO, pc.K, pc.RHO * pc.CP * t) for t in range(1, 2501))
    ini = (pc.INI.replace('iter-threshold = 1', 'iter-threshold = 50').replace('time-scheme = euler', 'time-scheme = RK3')
           .replace('sol-diter     = 1', 'sol-diter     = 50'))
    open(os.path.join(out, 'input.ini'), 'w').write(ini)


if __name__ == '__main__':
    if len(sys.argv) == 4 and sys.argv[1] == 'case' and sys.argv[2] in ('one', 'two', 'two-reversed'):
        write_case(sys.argv[2], sys.argv[3])
    else:
        sys.exit(__doc__)
