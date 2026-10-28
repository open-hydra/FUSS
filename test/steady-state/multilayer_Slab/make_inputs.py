"""Writes MESH/mesh.tec and INPUT/ of the multilayer_Slab test (no ATLAS: one block, the wall data are
the two end temperatures). Three layers in series through x inside ONE block, each with its own
material and its own cell width, so that every interface joins two materials AND two cell sizes:
coat k = 5 (8 cells of 0.3125 mm), liner k = 300 (6 cells of 1.25 mm), jacket k = 15 (10 cells of 0.5 mm).
Walls 302 at x = 0 (1000 K) and x = L (400 K); y and z faces adiabatic. The exact steady field is
linear in each layer, with the same heat flux q = (T_L - T_R) / sum(L_i / k_i) through all three."""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
LAYERS = [  # name, thickness [m], cells, k [W/(m K)], rho, cp
    ('coat',   0.0025, 8,  5.0,   6000.0, 500.0),
    ('liner',  0.0075, 6,  300.0, 8900.0, 390.0),
    ('jacket', 0.0050, 10, 15.0,  7900.0, 480.0),
]
T_LEFT, T_RIGHT, T0 = 1000.0, 400.0, 700.0
NY, LY, LZ = 2, 0.002, 0.001
TMAX = 2500

R_TOTAL = sum(t / k for _, t, _, k, _, _ in LAYERS)
Q = (T_LEFT - T_RIGHT) / R_TOTAL          # heat flux along +x [W/m2]


def nodes_x():
    xs = [0.0]
    for _, t, n, _, _, _ in LAYERS:
        x0 = xs[-1]
        xs += [x0 + t * (i + 1) / n for i in range(n)]
    return xs


def exact(x):
    """Exact steady temperature at x: T_LEFT minus q times the thermal resistance from x = 0."""
    r, x0 = 0.0, 0.0
    for _, t, _, k, _, _ in LAYERS:
        if x <= x0 + t:
            return T_LEFT - Q * (r + (x - x0) / k)
        r, x0 = r + t / k, x0 + t
    return T_RIGHT


def e(v):
    return f' {v:.15E}'


def main():
    for sub in ('MESH', 'INPUT'):
        os.makedirs(os.path.join(HERE, sub), exist_ok=True)
    xs = nodes_x()
    nx = len(xs) - 1
    ys = [LY * j / NY for j in range(NY + 1)]
    zs = [-LZ / 2, LZ / 2]
    nodes = ([x for z in zs for y in ys for x in xs] + [y for z in zs for y in ys for x in xs]
             + [z for z in zs for y in ys for x in xs])
    zone = (f' ZONE  T = Block1, I={nx + 1}, J={NY + 1}, K=2, DATAPACKING=BLOCK, VARLOCATION=([1-3]=NODAL), '
            'SOLUTIONTIME=0.0, STRANDID = 0\n')
    with open(os.path.join(HERE, 'MESH', 'mesh.tec'), 'w') as f:
        f.write(' VARIABLES ="x" "y" "z"\n' + zone + '\n'.join(e(v) for v in nodes) + '\n')
    mat = [float(m) for m, (_, _, n, _, _, _) in enumerate(LAYERS, start=1) for _ in range(n)]
    with open(os.path.join(HERE, 'INPUT', 'ic.tec'), 'w') as f:
        f.write(' VARIABLES ="x" "y" "z" T matID\n' + zone.replace('Block1', 'B1-SP').replace(
            'VARLOCATION=([1-3]=NODAL)', 'VARLOCATION=([1-3]=NODAL,[4-5]=CELLCENTERED)'))
        f.write('\n'.join(e(v) for v in nodes + [T0] * nx * NY + mat * NY) + '\n')
    iv = lambda *v: ''.join(f'{x:8d}' for x in v)
    lines = []
    for j in range(1, NY + 1):
        lines += [iv(1, 1, j, 1, 1, 302), f'    {T_LEFT:.15E},']
    for j in range(1, NY + 1):
        lines += [iv(1, nx, j, 1, 2, 302), f'    {T_RIGHT:.15E},']
    for f, jj in ((3, 1), (4, NY)):
        lines += [iv(1, i, jj, 1, f, 0) for i in range(1, nx + 1)]
    for f in (5, 6):
        lines += [iv(1, i, j, 1, f, 0) for j in range(1, NY + 1) for i in range(1, nx + 1)]
    open(os.path.join(HERE, 'INPUT', 'bc.txt'), 'w').write('\n'.join(lines) + '\n')
    open(os.path.join(HERE, 'INPUT', 'phase.txt'), 'w').write(
        'solid-bulk\n' + ''.join(f'{n} 1\n' for n, *_ in LAYERS))
    with open(os.path.join(HERE, 'INPUT', 'properties.dat'), 'w') as f:
        f.write('TITLE = "Mass Thermodynamic Properties"\n'
                'VARIABLES = "Temperature", "Cp", "Density", "Conductivity", "Energy"\n')
        for name, _, _, k, rho, cp in LAYERS:
            f.write(f'ZONE T="{name}"\nI={TMAX}, F=POINT\n')
            f.writelines(f'{t:.1f} {cp:.6f} {rho:.6f} {k:.6f} {rho * cp * t:.6f}\n' for t in range(1, TMAX + 1))


if __name__ == '__main__':
    main()
