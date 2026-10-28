"""Writes input.ini, MESH/mesh.tec and INPUT/ of the Skewed3D case (no ATLAS): one 3-D block of rhombohedral cells, the
three edge directions at 115 degrees to each other, split in i into two materials of different conductivity (so the
face-conductivity pass of the time step runs on the cells at the material interface), constant properties, walls 302 at
a constant temperature on the six faces. verify.py checks the step of every cell after one iteration."""
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
NX, NY, NZ = 6, 4, 4                    # cells; the material changes between i = 3 and i = 4
H = 0.05                                # edge length [m]
ANGLE = 115.0                           # angle between any two edges [deg]
K1, K2 = 20.0, 60.0                     # conductivity of the two materials [W/(m K)]
RHO, CP = 7800.0, 500.0
TWALL, T0 = 1000.0, 1300.0
VNN = 0.4


def edges():
    c = math.cos(math.radians(ANGLE))
    s2 = math.sqrt(1 - c * c)
    e1 = (1.0, 0.0, 0.0)
    e2 = (c, s2, 0.0)
    y3 = (c - c * c) / s2
    e3 = (c, y3, math.sqrt(1 - c * c - y3 * y3))
    return [tuple(H * v for v in e) for e in (e1, e2, e3)]


def e(v):
    return f' {v:.15E}'


def main():
    e1, e2, e3 = edges()
    nodes = [[[tuple(i * e1[d] + j * e2[d] + k * e3[d] for d in range(3)) for i in range(NX + 1)]
              for j in range(NY + 1)] for k in range(NZ + 1)]
    flat = [n[d] for d in range(3) for plane in nodes for row in plane for n in row]
    zone = f' ZONE  T = Block1, I={NX + 1}, J={NY + 1}, K={NZ + 1}, DATAPACKING=BLOCK, VARLOCATION=([1-3]=NODAL), SOLUTIONTIME=0.0, STRANDID = 0\n'
    os.makedirs(os.path.join(HERE, 'MESH'), exist_ok=True)
    os.makedirs(os.path.join(HERE, 'INPUT'), exist_ok=True)
    with open(os.path.join(HERE, 'MESH', 'mesh.tec'), 'w') as f:
        f.write(' VARIABLES ="x" "y" "z"\n' + zone + '\n'.join(e(v) for v in flat) + '\n')
    ncell = NX * NY * NZ
    mat = [1.0 if i < NX // 2 else 2.0 for k in range(NZ) for j in range(NY) for i in range(NX)]
    with open(os.path.join(HERE, 'INPUT', 'ic.tec'), 'w') as f:
        f.write(' VARIABLES ="x" "y" "z" T matID\n' + zone.replace('Block1', 'B1-SP').replace(
            'VARLOCATION=([1-3]=NODAL)', 'VARLOCATION=([1-3]=NODAL,[4-5]=CELLCENTERED)'))
        f.write('\n'.join(e(v) for v in flat + [T0] * ncell + mat) + '\n')
    lines = []
    for face in range(1, 7):
        if face in (1, 2):
            cells = [(1 if face == 1 else NX, j, k) for k in range(1, NZ + 1) for j in range(1, NY + 1)]
        elif face in (3, 4):
            cells = [(i, 1 if face == 3 else NY, k) for k in range(1, NZ + 1) for i in range(1, NX + 1)]
        else:
            cells = [(i, j, 1 if face == 5 else NZ) for j in range(1, NY + 1) for i in range(1, NX + 1)]
        for (i, j, k) in cells:
            lines += [''.join(f'{v:8d}' for v in (1, i, j, k, face, 302)), f'    {TWALL:.15E},']
    open(os.path.join(HERE, 'INPUT', 'bc.txt'), 'w').write('\n'.join(lines) + '\n')
    open(os.path.join(HERE, 'INPUT', 'phase.txt'), 'w').write('solid-bulk\nmat1 1\nmat2 1\n')
    with open(os.path.join(HERE, 'INPUT', 'properties.dat'), 'w') as f:
        f.write('TITLE = "Mass Thermodynamic Properties"\nVARIABLES = "Temperature", "Cp", "Density", "Conductivity", "Energy"\n')
        for name, k in (('mat1', K1), ('mat2', K2)):
            f.write(f'ZONE T="{name}"\nI=2500, F=POINT\n')
            f.writelines(f'{t:.1f} {CP:.6f} {RHO:.6f} {k:.6f} {RHO * CP * t:.6f}\n' for t in range(1, 2501))
    with open(os.path.join(HERE, 'input.ini'), 'w') as f:
        f.write(f"""; Time step of 3-D cells skewed in all three pairs of directions: MESH/ and INPUT/ come from make_inputs.py, not ATLAS

[FUSS-Parameters]
newrun = true
res-threshold = 1e-30
iter-threshold = 1

[FUSS-Numerics]
time-scheme = euler
vnn = {VNN}
time-accurate = false
integration-variables = cons

[FUSS-IO]
sol-diter     = 1
sol-overwrite = true
sol-format    = tecplot ascii
res-diter = 1
shell-diter = 1
ini-diter = 100000000
""")


if __name__ == '__main__':
    main()
