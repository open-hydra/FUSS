"""Writes MESH/mesh.tec and INPUT/ of the skewed-plate patch test (no ATLAS: the wall data come from the
exact solution). One block of N x N cells, x sheared by S y; exact T = A + B x + C y. Walls: face 1 (left)
302 with Tw, face 2 (right) 303 with Tref, face 3 (bottom) 301 with qw, face 4 (top) 305 with Tref, each
cell's value set so that the exact solution satisfies the condition."""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
N, S = 20, 1.0                                    # cells per side, shear (45 deg)
A, B, C = 1000.0, 200.0, 300.0                    # T = A + B x + C y [K]
K, RHO, CP = 52.0, 7850.0, 460.0                  # constant properties
H, EPS, SIGMA = 100.0, 0.8, 5.67e-8               # 303 and 305 coefficients
T0 = 1300.0                                       # initial temperature


def exact(x, y):
    return A + B * x + C * y


def face_centre(f, i, j):
    return {1: (S * (j - .5) / N, (j - .5) / N), 2: (1 + S * (j - .5) / N, (j - .5) / N),
            3: ((i - .5) / N, 0.0), 4: ((i - .5) / N + S, 1.0)}[f]


R = (1 + S * S) ** .5
Q_IN = {1: K * (-B + S * C) / R, 2: K * (B - S * C) / R, 3: -K * C, 4: K * C}   # k grad(T).n_out [W/m2]


def tref_305(tw, q):
    """Tref with H (Tref - Tw) + sigma eps (Tref^4 - Tw^4) = q."""
    x = tw
    for _ in range(60):
        x -= (H * (x - tw) + SIGMA * EPS * (x ** 4 - tw ** 4) - q) / (H + 4 * SIGMA * EPS * x ** 3)
    return x


def e(v):
    return f' {v:.15E}'


def main():
    xs = [[[i / N + S * j / N for i in range(N + 1)] for j in range(N + 1)] for _ in range(2)]
    ys = [[[j / N for i in range(N + 1)] for j in range(N + 1)] for _ in range(2)]
    zs = [[[z for i in range(N + 1)] for j in range(N + 1)] for z in (-0.005, 0.005)]
    nodes = [v for arr in (xs, ys, zs) for plane in arr for row in plane for v in row]
    zone = (f' ZONE  T = Block1, I={N + 1}, J={N + 1}, K=2, DATAPACKING=BLOCK, VARLOCATION=([1-3]=NODAL), '
            'SOLUTIONTIME=0.0, STRANDID = 0\n')
    with open(os.path.join(HERE, 'MESH', 'mesh.tec'), 'w') as f:
        f.write(' VARIABLES ="x" "y" "z"\n' + zone + '\n'.join(e(v) for v in nodes) + '\n')
    with open(os.path.join(HERE, 'INPUT', 'ic.tec'), 'w') as f:
        f.write(' VARIABLES ="x" "y" "z" T matID\n' + zone.replace('Block1', 'B1-SP').replace(
            'VARLOCATION=([1-3]=NODAL)', 'VARLOCATION=([1-3]=NODAL,[4-5]=CELLCENTERED)'))
        f.write('\n'.join(e(v) for v in nodes + [T0] * N * N + [1.0] * N * N) + '\n')
    lines = []
    for f in (1, 2, 3, 4):
        for c in range(1, N + 1):
            i, j = {1: (1, c), 2: (N, c), 3: (c, 1), 4: (c, N)}[f]
            tw, q = exact(*face_centre(f, i, j)), Q_IN[f]
            t, vals = {1: (302, [tw]), 2: (303, [H, 0.0, tw + q / H]), 3: (301, [q]),
                       4: (305, [H, EPS, tref_305(tw, q)])}[f]
            lines += [''.join(f'{v:8d}' for v in (1, i, j, 1, f, t)), ''.join(f'    {v:.15E},' for v in vals)]
    for f in (5, 6):
        lines += [''.join(f'{v:8d}' for v in (1, i, j, 1, f, 0)) for j in range(1, N + 1) for i in range(1, N + 1)]
    open(os.path.join(HERE, 'INPUT', 'bc.txt'), 'w').write('\n'.join(lines) + '\n')
    open(os.path.join(HERE, 'INPUT', 'phase.txt'), 'w').write('solid-bulk\nmat 1\n')
    with open(os.path.join(HERE, 'INPUT', 'properties.dat'), 'w') as f:
        f.write('TITLE = "Mass Thermodynamic Properties"\n'
                'VARIABLES = "Temperature", "Cp", "Density", "Conductivity", "Energy"\nZONE T="mat"\nI=2500, F=POINT\n')
        f.writelines(f'{t:.1f} {CP:.6f} {RHO:.6f} {K:.6f} {RHO * CP * t:.6f}\n' for t in range(1, 2501))


if __name__ == '__main__':
    main()
