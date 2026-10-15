"""Inputs of the multimat_liner case: a 4-layer liner wall in 2 blocks, material jumps inside both.

Slab of 64 x 128 cells on [0, 0.02] x [0, 0.5] (x through the thickness, cells 0.3125 mm x 3.9 mm),
four material layers through x: coat 2.5 mm (k = 5), liner 7.5 mm (k = 300), jacket 5 mm
(k = 15), shell 5 mm (k = 150) W/m/K. Hot side x = 0: convection 303, h = 2000, Tref falling
3000 -> 1500 K along y; cold side x = 0.02: convection 303, h = 10000, Tref = 300 K; y ends
adiabatic. The slab is split along y into two 64 x 64 blocks joined by a block connection (101),
so on 2 ranks each owns a block and all material jumps lie inside the blocks. bc<m>.txt is written
for m = 1..4 (LEVELS), the boundary conditions of each multigrid level sampled at its own face
centres: the blocks coarsen to 8 x 8 and the layer boundaries (8, 32, 48 cells through x) stay on
grid lines of every level. The shipped input.ini runs one level; the FAS figure campaign
(convergence_figure/) runs the case on 4.

The interior faces between two materials carry the series conductivity of their half cells and
the time step takes the face conductivities: the run goes at vnn 1.0. Before October 2026 the face
took the k of the lower cell (300 at the liner -> jacket jump, cells i = 32 | 33, 20 times the
jacket's own) and with the cell-k time step the run died within a few iterations at vnn 0.5.
Derived from the FAS stress-campaign liner (fuss-deps/vfy/liner_case.py), single block there.
"""
import os

NX, NY, NB, LEVELS = 64, 128, 2, 4
NYB = NY // NB
LX, LY, LZ = 0.02, 0.5, 0.01
MATS = [  # name, x_end / LX, rho, cp, k
    ('coat',   0.125, 6000.0, 500.0, 5.0),
    ('liner',  0.500, 8900.0, 390.0, 300.0),
    ('jacket', 0.750, 7900.0, 480.0, 15.0),
    ('shell',  1.000, 2800.0, 900.0, 150.0),
]
H_HOT, TREF_HOT0, TREF_HOT1 = 2000.0, 3000.0, 1500.0
H_COLD, TREF_COLD = 10000.0, 300.0
T0, TMAX = 800.0, 3500

INI = """[FUSS-Parameters]
newrun = true
res-threshold = 1e-8
iter-threshold = 1500000

[FUSS-Numerics]
time-scheme = RK3
vnn = 1.0
time-accurate = false
integration-variables = cons

[FUSS-IO]
sol-diter     = 100000000
sol-overwrite = true
sol-format    = tecplot ascii
ini-diter = 100000000
res-diter = 100
shell-diter = 1000
"""


def e(v):
    return f' {v:.15E}'


def iv(*v):
    return ''.join(f'{x:8d}' for x in v)


def rv(*v):
    return ''.join(f'    {x:.15E},' for x in v)


def matid(xh):
    for m, (_, xe, _, _, _) in enumerate(MATS, start=1):
        if xh < xe:
            return float(m)
    return float(len(MATS))


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    for sub in ('MESH', 'INPUT'):
        os.makedirs(os.path.join(here, sub), exist_ok=True)

    mesh, ic = [' VARIABLES ="x" "y" "z"\n'], [' VARIABLES ="x" "y" "z" T matID\n']
    for b in range(NB):
        j0 = b * NYB
        nodes = ([i / NX * LX for j in range(NYB + 1) for i in range(NX + 1)] * 2
                 + [(j0 + j) / NY * LY for j in range(NYB + 1) for i in range(NX + 1)] * 2
                 + [z for z in (-LZ / 2, LZ / 2) for _ in range((NX + 1) * (NYB + 1))])
        zone = (f' ZONE  T = Block{b + 1}, I={NX + 1}, J={NYB + 1}, K=2, DATAPACKING=BLOCK, '
                'VARLOCATION=([1-3]=NODAL), SOLUTIONTIME=0.0, STRANDID = 0\n')
        mesh += [zone, '\n'.join(e(v) for v in nodes) + '\n']
        cells = [T0] * (NX * NYB) + [matid((i + .5) / NX) for j in range(NYB) for i in range(NX)]
        ic += [zone.replace(f'Block{b + 1}', f'B{b + 1}-SP').replace(
            'VARLOCATION=([1-3]=NODAL)', 'VARLOCATION=([1-3]=NODAL,[4-5]=CELLCENTERED)'),
            '\n'.join(e(v) for v in nodes + cells) + '\n']
    open(os.path.join(here, 'MESH', 'mesh.tec'), 'w').write(''.join(mesh))
    open(os.path.join(here, 'INPUT', 'ic.tec'), 'w').write(''.join(ic))

    for m in range(1, LEVELS + 1):  # bc<m>.txt: level m, its own cell counts and face centres
        nx, nyb, ny = NX >> (m - 1), NYB >> (m - 1), NY >> (m - 1)
        lines = []
        for b in range(1, NB + 1):
            for j in range(1, nyb + 1):  # face 1: hot gas, Tref falls along y
                yh = ((b - 1) * nyb + j - .5) / ny
                lines += [iv(b, 1, j, 1, 1, 303), rv(H_HOT, 0.0, TREF_HOT0 + (TREF_HOT1 - TREF_HOT0) * yh)]
            for j in range(1, nyb + 1):  # face 2: coolant
                lines += [iv(b, nx, j, 1, 2, 303), rv(H_COLD, 0.0, TREF_COLD)]
            for f, jj in ((3, 1), (4, nyb)):  # y ends adiabatic, the block interface a connection
                for i in range(1, nx + 1):
                    if f == 4 and b < NB:
                        lines += [iv(b, i, jj, 1, f, 101), iv(b + 1, i, 1, 1, 3, 1, 0, 0, 1)]
                    elif f == 3 and b > 1:
                        lines += [iv(b, i, jj, 1, f, 101), iv(b - 1, i, nyb, 1, 4, 1, 0, 0, 1)]
                    else:
                        lines += [iv(b, i, jj, 1, f, 0)]
            for f in (5, 6):
                lines += [iv(b, i, j, 1, f, 0) for j in range(1, nyb + 1) for i in range(1, nx + 1)]
        open(os.path.join(here, 'INPUT', 'bc.txt' if m == 1 else f'bc{m}.txt'), 'w').write('\n'.join(lines) + '\n')

    open(os.path.join(here, 'INPUT', 'phase.txt'), 'w').write(
        'solid-bulk\n' + ''.join(f'{n} 1\n' for n, *_ in MATS))
    with open(os.path.join(here, 'INPUT', 'properties.dat'), 'w') as f:
        f.write('TITLE = "Mass Thermodynamic Properties"\n'
                'VARIABLES = "Temperature", "Cp", "Density", "Conductivity", "Energy"\n')
        for name, _, rho, cp, k in MATS:
            f.write(f'ZONE T="{name}"\nI={TMAX}, F=POINT\n')
            f.writelines(f'{t:.1f} {cp:.6f} {rho:.6f} {k:.6f} {rho * cp * t:.6f}\n'
                         for t in range(1, TMAX + 1))
    open(os.path.join(here, 'input.ini'), 'w').write(INI)


if __name__ == '__main__':
    main()
