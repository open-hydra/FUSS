"""Inputs of the chimera_mg case: two overlapping single-material blocks.

Block1 covers [0, 0.56] x [0, 0.1] with 112 x 16 cells, Block2 covers [0.47, 1.0] x [0, 0.1]
with 104 x 12 cells: 0.09 m of genuine overlap and different spacings in BOTH directions
(equal dx would make boundary facelet centres of the two blocks coincide and the conformal
search would pair them as 101 connections). Every dimension is divisible by 4 and the
coarsest level keeps >= 3 cells: the case supports 3 multigrid levels. Writes MESH/mesh.tec,
INPUT/{ic.tec,phase.txt,properties.dat}; bc<m>.txt come from ATLAS BCB (BC-chimera = true,
MG-levels = 3 in input.ini), which also fills SOLUTION/ with its chimera diagnostics.
"""
import os

BLOCKS = [  # name, x0, lx, nx, ny
    ('Block1', 0.00, 0.56, 112, 16),
    ('Block2', 0.47, 0.53, 104, 12),
]
LY, LZ = 0.1, 0.01
T0, TMAX = 1500.0, 3000
K, RHO, CP = 10.0, 8000.0, 450.0


def e(v):
    return f' {v:.15E}'


def zone(name, nx, ny, extra=''):
    return (f' ZONE  T = {name}, I={nx + 1}, J={ny + 1}, K=2, DATAPACKING=BLOCK, '
            f'VARLOCATION=([1-3]=NODAL{extra}), SOLUTIONTIME=0.0, STRANDID = 0\n')


def nodes(x0, lx, nx, ny):
    xs = [x0 + i / nx * lx for j in range(ny + 1) for i in range(nx + 1)] * 2
    ys = [j / ny * LY for j in range(ny + 1) for i in range(nx + 1)] * 2
    zs = [z for z in (-LZ / 2, LZ / 2) for _ in range((nx + 1) * (ny + 1))]
    return xs + ys + zs


here = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(here, 'MESH'), exist_ok=True)
os.makedirs(os.path.join(here, 'INPUT'), exist_ok=True)
os.makedirs(os.path.join(here, 'SOLUTION'), exist_ok=True)

with open(os.path.join(here, 'MESH', 'mesh.tec'), 'w') as f:
    f.write(' VARIABLES ="x" "y" "z"\n')
    for name, x0, lx, nx, ny in BLOCKS:
        f.write(zone(name, nx, ny) + '\n'.join(e(v) for v in nodes(x0, lx, nx, ny)) + '\n')

with open(os.path.join(here, 'INPUT', 'ic.tec'), 'w') as f:
    f.write(' VARIABLES ="x" "y" "z" T matID\n')
    for bn, (name, x0, lx, nx, ny) in enumerate(BLOCKS, start=1):
        f.write(zone(f'B{bn}-SP', nx, ny, ',[4-5]=CELLCENTERED'))
        f.write('\n'.join(e(v) for v in nodes(x0, lx, nx, ny)
                          + [T0] * (nx * ny) + [1.0] * (nx * ny)) + '\n')

open(os.path.join(here, 'INPUT', 'phase.txt'), 'w').write('solid-bulk\nmat1 1\n')

with open(os.path.join(here, 'INPUT', 'properties.dat'), 'w') as f:
    f.write('TITLE = "Mass Thermodynamic Properties"\n'
            'VARIABLES = "Temperature", "Cp", "Density", "Conductivity", "Energy"\n'
            f'ZONE T="mat1"\nI={TMAX}, F=POINT\n')
    for t in range(1, TMAX + 1):
        f.write(f'{t:.1f} {CP:.6f} {RHO:.6f} {K:.6f} {RHO * CP * t:.6f}\n')
