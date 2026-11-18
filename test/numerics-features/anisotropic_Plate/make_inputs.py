"""Inputs of the anisotropic_Plate case: 128 x 16 cells on [0,1] x [0,0.00125].

Every cell has aspect ratio 100 (dx = 7.8e-3, dy = 7.8e-5): the reference case for the
documented full-coarsening multigrid limit, where V/W cycles gain little because coarsening
in the stiff direction cannot represent the strong y coupling (the honest rates are on the
docs V&V page). All six faces are isothermal (302) with the exact values of the linear field
T = 1000 + 200 xh + 300 yh in NORMALIZED coordinates, sampled at the face centres of each of
the 3 multigrid levels (bc.txt, bc2.txt, bc3.txt), so the discrete steady state is the exact
linear field on every level. Derived from the FAS verification generator (fuss-deps/vfy/
make_case.py --nx 128 --ny 16 --ly 0.00125 --levels 3 --simple-bc).
"""
import os

NX, NY = 128, 16
LX, LY, LZ = 1.0, 0.00125, 0.01
LEVELS = 3
A, B, C = 1000.0, 200.0, 300.0
K, RHO, CP = 52.0, 7850.0, 460.0
TMIN, TMAX, T0 = 1, 2500, 1300.0


def e(v):
    return f' {v:.15E}'


def bc_lines(nx, ny):
    lines = []
    for f in (1, 2, 3, 4):
        ncell = ny if f in (1, 2) else nx
        for c in range(1, ncell + 1):
            i, j = {1: (1, c), 2: (nx, c), 3: (c, 1), 4: (c, ny)}[f]
            xh, yh = {1: (0.0, (j - .5) / ny), 2: (1.0, (j - .5) / ny),
                      3: ((i - .5) / nx, 0.0), 4: ((i - .5) / nx, 1.0)}[f]
            lines += [''.join(f'{v:8d}' for v in (1, i, j, 1, f, 302)),
                      f'    {A + B * xh + C * yh:.15E},']
    for f in (5, 6):
        lines += [''.join(f'{v:8d}' for v in (1, i, j, 1, f, 0))
                  for j in range(1, ny + 1) for i in range(1, nx + 1)]
    return lines


here = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(here, 'MESH'), exist_ok=True)
os.makedirs(os.path.join(here, 'INPUT'), exist_ok=True)

nodes = ([i / NX * LX for j in range(NY + 1) for i in range(NX + 1)] * 2
         + [j / NY * LY for j in range(NY + 1) for i in range(NX + 1)] * 2
         + [z for z in (-LZ / 2, LZ / 2) for _ in range((NX + 1) * (NY + 1))])
zone = (f' ZONE  T = Block1, I={NX + 1}, J={NY + 1}, K=2, DATAPACKING=BLOCK, '
        'VARLOCATION=([1-3]=NODAL), SOLUTIONTIME=0.0, STRANDID = 0\n')
with open(os.path.join(here, 'MESH', 'mesh.tec'), 'w') as f:
    f.write(' VARIABLES ="x" "y" "z"\n' + zone + '\n'.join(e(v) for v in nodes) + '\n')

with open(os.path.join(here, 'INPUT', 'ic.tec'), 'w') as f:
    f.write(' VARIABLES ="x" "y" "z" T matID\n' + zone.replace('Block1', 'B1-SP').replace(
        'VARLOCATION=([1-3]=NODAL)', 'VARLOCATION=([1-3]=NODAL,[4-5]=CELLCENTERED)'))
    f.write('\n'.join(e(v) for v in nodes + [T0] * (NX * NY) + [1.0] * (NX * NY)) + '\n')

for m in range(1, LEVELS + 1):
    fn = 'bc.txt' if m == 1 else f'bc{m}.txt'
    open(os.path.join(here, 'INPUT', fn), 'w').write(
        '\n'.join(bc_lines(NX >> (m - 1), NY >> (m - 1))) + '\n')

open(os.path.join(here, 'INPUT', 'phase.txt'), 'w').write('solid-bulk\nmat 1\n')

with open(os.path.join(here, 'INPUT', 'properties.dat'), 'w') as f:
    f.write('TITLE = "Mass Thermodynamic Properties"\n'
            'VARIABLES = "Temperature", "Cp", "Density", "Conductivity", "Energy"\n'
            f'ZONE T="mat"\nI={TMAX - TMIN + 1}, F=POINT\n')
    for t in range(TMIN, TMAX + 1):
        f.write(f'{t:.1f} {CP:.6f} {RHO:.6f} {K:.6f} {RHO * CP * t:.6f}\n')
