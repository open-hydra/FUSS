#!/usr/bin/env python3
"""Periodic faces (201) with ours' moving mesh and type-2 pyrolysis (plan 14 S9d, Part B): the decks.

  generate.py <variant> <kind> <out dir>    one case (the physics gate writes the sectors and the rings with it)
  generate.py --decks [variant,...]          rewrite the committed decks (of the named variants; default every one)

The geometry, the boundary values and the initial field are those of upstream's PeriodicSector
(test/numerics-features/periodic_sector/make_sector.py, imported, not edited): a sector of a helical annular wall about
the x axis, 6 sectors in the ring, 6 x 6 x 4 cells (x, theta, r) per sector, the side faces periodic (201) in a sector,
block connections (101) in the ring. Kinds: sector (blocks A, B), ring (12 blocks), sector1 (one block periodic with
itself), ring1 (6 blocks), sym (sector with symmetry faces 300 instead of the pair: the control that must not give the
ring). This script changes the run (time-accurate, a moving mesh, the inner wall, pyrolysis):

  translate   law = translation along the axis of the pair (x): the copies of the sector move with it
  prescribed  law = prescribed, a displacement along x that depends on x only (amp, kx of the first row)
  receding    law = nodeshift, the inner wall (face 5) an ablating wall 306 at a uniform rate (sdot-override)
  darcy       gas-closure = type2 on a static mesh, the inner wall 303 open to the gas, K transversely isotropic about x
  ablation    type2 and nodeshift, the inner wall 306 with the TACOT B' table, its env pattern periodic in theta, the
              outer face adiabatic
  ablation-tight  ablation with [FUSS-Darcy] tol = 1e-13 (the fallback RED of the join's rate and windows; never committed)
Variants for checks by hand (never committed; run_guard_demos.sh edits the committed decks for its across, wobble and
aniso demos and writes only the openseam base with this script):
  across      translate with the velocity across the axis (stops at set-up)
  wobble      prescribed with the displacement across the axis (stops at the first step)
  aniso       darcy with K transversely isotropic about z (not invariant under the pair's rotation: stops at set-up)
  axis        receding on a sector of a full cylinder (R0 = 0), the end face x = 0 receding: stops at Check_Folded_Cells,
              before the seam's own-image refusal it was meant for (plan 14 S9d Part B, O-AXIS); its axis node is its own
              periodic image
  openseam    receding, sector1, the last theta column of the inner wall a 301: the receding patch is not closed
              across the periodic seam (stops at set-up)
"""
import os
import shutil
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'test', 'numerics-features', 'periodic_sector'))
import make_sector as ms  # noqa: E402

PYRO_SRC = os.path.join(ROOT, 'test', 'pyrolysis', 'ablation', 'smoothing', 'w5_2blk', 'INPUT')
PYRO_FILES = ('pyrolysis.dat', 'virgin.dat', 'char.dat', 'gas.dat', 'properties.dat', 'phase.txt')
FUSS_SH = os.path.join(ROOT, 'test', 'moving-mesh', 'receding_4blk', 'FUSS.sh')
ITERS = 200
ITERS_OF = {'ablation': 20}
# The ablation deck (plan 14 S9d Part B, M-B3-d, measured on the ring): env tables = w5_2blk's env_hi.dat with rhoeUeCH0
# x 0.3 (hi) and x 0.15 (lo), and the outer face (6; upstream's 301 at -3.0e4 W/m^2) adiabatic (inner_wall).
#  - w5_2blk's own tables (0.33 / 0.27): every variant (sector, sym, ring) stopped "mesh tangled" after 31 steps, at
#    ~3.4 mm of recession: on make_sector's randomly perturbed mesh the first line edges are up to 53 degrees off the
#    receding direction, nodeshift moves the wall nodes along them and the first-layer cells shear until one folds (the
#    min-thickness is not reached; a receding deck on this mesh must stay well below that depth);
#  - x 0.1 / x 0.05 with the -3.0e4 outer face: the outer layer cooled from 446 K to 261 K in 16 steps and the run
#    stopped at step 17 ("energy outside table range", block 2 cell (2,2,4); the tables start at 256 K), the wall at
#    780-960 K, the low side's chemistry nearly off (rate 23x lower than the high side's);
#  - these tables, the outer face adiabatic (step 20 of 20, t ~15 s): Tw 1160-1203 K (hi) / 990-1019 K (lo), the outer
#    layer >= 506 K, recession <= 0.72 mm (hi) / 0.25 mm (lo) (21 % of the tangle depth), rate and recession ~3x apart
#    across the seam (m-B2-3/4 need a rate that differs there), min scaled corner Jacobian 0.303 (0.343 at the start).
ENV_SCALE = {'env_hi.dat': 0.3, 'env_lo.dat': 0.15}

DECKS = (   # committed decks: (variant, kind, directory relative to the root)
    ('translate', 'sector', 'test/moving-mesh/periodic/translate'),
    ('prescribed', 'sector', 'test/moving-mesh/periodic/prescribed'),
    ('receding', 'sector', 'test/moving-mesh/periodic/receding'),
    ('receding', 'sector1', 'test/moving-mesh/periodic/receding1'),
    ('darcy', 'sector', 'test/pyrolysis/darcy/periodic_static'),
    ('darcy', 'sector1', 'test/pyrolysis/darcy/periodic_static1'),
    ('ablation', 'sector', 'test/pyrolysis/ablation/periodic_receding'),
)

HEAD = """[FUSS-Parameters]
newrun = true
res-threshold = 1e-30
time-threshold = 1.0e6
iter-threshold = {iters}

[FUSS-Numerics]
time-scheme = RK3
vnn = 0.25
time-accurate = true
integration-variables = cons
"""
TAIL = """
[FUSS-IO]
sol-diter     = {iters}
sol-overwrite = true
sol-format    = tecplot ascii
res-diter = 1
shell-diter = 50
ini-diter = 100000000
"""
MOTION = {
    'translate': "\n[FUSS-MeshMotion]\nlaw = translation\nvel = 2.0e-3 0.0 0.0\ngcl-tolerance = 1.0e-10\n",
    'across': "\n[FUSS-MeshMotion]\nlaw = translation\nvel = 0.0 2.0e-3 0.0\ngcl-tolerance = 1.0e-10\n",
    'prescribed': ("\n[FUSS-MeshMotion]\nlaw = prescribed\namp = 2.0e-4 0.0 0.0\nkx  = 60.0 0.0 0.0\n"
                   "ky  = 0.0 0.0 0.0\nkz  = 0.0 0.0 0.0\nomega = 0.5\ngcl-tolerance = 1.0e-10\n"),
    'wobble': ("\n[FUSS-MeshMotion]\nlaw = prescribed\namp = 0.0 2.0e-4 0.0\nkx  = 0.0 0.0 0.0\n"
               "ky  = 0.0 0.0 0.0\nkz  = 0.0 0.0 0.0\nomega = 0.5\ngcl-tolerance = 1.0e-10\n"),
    'receding': "\n[FUSS-MeshMotion]\nlaw = nodeshift\n\n[FUSS-Ablation]\nsdot-override = 1.0e-4\n",
    'darcy': "",
    'ablation': "\n[FUSS-MeshMotion]\nlaw = nodeshift\n\n[FUSS-Ablation]\nsmoothing-width = 5\n",
}
MOTION['axis'] = MOTION['openseam'] = MOTION['receding']
MOTION['aniso'] = MOTION['darcy']
MOTION['ablation-tight'] = MOTION['ablation']
# [FUSS-Darcy] tol = 1.0e-10 in every type-2 deck (darcy, aniso, ablation; ablation-tight 1e-13), plan 14 S9d Part B
# M-B3-a: a sector-vs-ring comparison of an iteratively solved field is only as sharp as the iteration. At the default
# 1e-8 the one-block deck (sector1: its closed chain sweeps otherwise than its ring) measured the iteration tolerance,
# and the two-block decks agreed to ~1e-13 only because they sweep as their ring does (a compiler change can break
# that); at 1e-10 the iteration sits on its floor (the sweeps of 1e-12), so every leg measures round-off. Sector vs
# ring, T relative (two blocks | one block | ablation): tol 1e-8 7.1e-13 | 2.9e-10 | 1.8e-13; 1e-9 1.2e-12 | 1.5e-11 |
# 1.8e-13; 1e-10 3.2e-13 | 6.7e-12 | 9.0e-12; 3e-11 1.5e-12 | 9.9e-12 | 2.5e-12; 1e-11 and 1e-12 1.2e-12 | 9.5e-12 |
# 2.5e-12.
PYRO = """
[FUSS-Pyrolysis]
gas-closure = type2

[FUSS-Darcy]
solver = tridiag
p-init = 101325.0
p-wall = 101325.0
"""
ENV_R = ("// t p_w rhoeUeCH0 h_r q_rad chem hconv T_edge   (SI; plan 08 section 2.1)\n"
         "0.0 101325.0 0.0 0.0 0.0 0.0 1000.0 1000.0\n1000000.0 101325.0 0.0 0.0 0.0 0.0 1000.0 1000.0\n")
K_AXIS_X = ("K_v 3.2e-11 0 0 0 1.6e-11 0 0 0 1.6e-11", "K_c 4e-11 0 0 0 2e-11 0 0 0 2e-11")
K_AXIS_Z = ("K_v 1.6e-11 0 0 0 1.6e-11 0 0 0 3.2e-11", "K_c 2e-11 0 0 0 2e-11 0 0 0 4e-11")
WITH_LINE = {101, 103, 201, 301, 302, 303, 304, 305, 306, 307}


def read_records(path):
    """bc.txt as a list of [ints of the record line, data line or None]"""
    lines = open(path).read().split('\n')
    out, i = [], 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        rec = [int(v) for v in lines[i].split()]
        data = None
        if rec[5] in WITH_LINE:
            data = lines[i + 1]
            i += 1
        out.append([rec, data])
        i += 1
    return out


def write_records(path, recs):
    out = []
    for rec, data in recs:
        out.append(''.join('%8d' % v for v in rec))
        if data is not None:
            out.append(data)
    open(path, 'w').write('\n'.join(out) + '\n')


def inner_wall(recs, variant, blks):
    """face 5 (the inner radius, ms 'hot' 303): per variant"""
    for r in recs:
        b, i, j, k, f, code = r[0]
        if f != 5:
            continue
        if variant in ('receding', 'axis', 'openseam'):
            r[0][5], r[1] = 306, '    0.0, 300.0, none, env_r.dat'
        elif variant == 'ablation':
            s = ms.Blk.face_base(blks[b - 1], 5, i, j, k)[2] % 1.0      # the fraction of the sector of the face
            r[0][5], r[1] = 306, '    0.5, 300.0, bprime_tacot_earth.dat, %s' % ('env_hi.dat' if s < 1.0 / 3.0 else 'env_lo.dat')
    if variant == 'ablation':      # the outer face (6, upstream's 301 at -3.0e4 W/m^2) adiabatic: see ENV_SCALE
        for r in recs:
            if r[0][4] == 6 and r[0][5] == 301:
                r[1] = '    %s,' % '0.000000000000000E+00'
    if variant == 'openseam':      # sector1: the last theta column of the inner wall a 301 wall (heat flux 0)
        for r in recs:
            if r[0][4] == 5 and r[0][2] == ms.NT:
                r[0][5], r[1] = 301, '    %s,' % '0.000000000000000E+00'
    if variant == 'axis':          # R0 = 0: face 5 lies on the axis (zero area): symmetry; face 1 (x = 0) recedes
        for r in recs:
            if r[0][4] == 5:
                r[0][5], r[1] = 300, None
            elif r[0][4] == 1:
                r[0][5], r[1] = 306, '    0.0, 300.0, none, env_r.dat'
    return recs


def fuss_sh(out):
    depth = os.path.relpath(out, ROOT).count(os.sep) + 1
    text = open(FUSS_SH).read().replace('MASTERDIR=../../..\n', 'MASTERDIR=%s\n' % '/'.join(['..'] * depth), 1)
    open(os.path.join(out, 'FUSS.sh'), 'w').write(text)
    os.chmod(os.path.join(out, 'FUSS.sh'), 0o755)


def write(variant, kind, out):
    # ablation-tight: the ablation deck with [FUSS-Darcy] tol = 1e-13 (plan 14 S9d Part B section 5, G-DARCY's rule: the
    # fallback RED of mutations m-B2-3 / m-B2-4 when the ablation leg's floor hides them); never committed
    tight = variant == 'ablation-tight'
    if tight:
        variant = 'ablation'
    r0, cart0 = ms.R0, ms.cart
    if variant == 'axis':
        ms.R0 = 0.0
    if kind in ('ring', 'ring1'):
        # the ring's last seam (sector NS-1, s = 1) on the doubles of sector 0's face s = 0: cos and sin of 2 pi and of 0
        # differ by ~1e-17 m, which the surface-node map (exact coordinates) and Assert_Interface_Nodes (bitwise) take
        # for a torn block interface (plan 14 S9d M5 probe (d) 1: the ring tears as generated, runs with an exact seam)
        ms.cart = lambda p, sector: cart0(p, sector - ms.NS) if sector + p[2] >= ms.NS else cart0(p, sector)
    try:
        ms.write_case(kind, out, 'stencil')
        blks = ms.blocks_of(kind)
    finally:
        ms.R0, ms.cart = r0, cart0
    bc = os.path.join(out, 'INPUT', 'bc.txt')
    write_records(bc, inner_wall(read_records(bc), variant, blks))
    iters = ITERS_OF.get(variant, ITERS)
    ini = HEAD.format(iters=iters) + MOTION[variant]
    if variant in ('darcy', 'aniso', 'ablation'):
        ini += PYRO + ('tol = 1.0e-13\n' if tight else 'tol = 1.0e-10\n')
        for f in PYRO_FILES:
            shutil.copyfile(os.path.join(PYRO_SRC, f), os.path.join(out, 'INPUT', f))
        p = os.path.join(out, 'INPUT', 'pyrolysis.dat')
        kv, kc = K_AXIS_Z if variant == 'aniso' else K_AXIS_X
        text = [kv if l.startswith('K_v ') else kc if l.startswith('K_c ') else l for l in open(p).read().split('\n')]
        open(p, 'w').write('\n'.join(text))
        if variant == 'ablation':
            shutil.copyfile(os.path.join(PYRO_SRC, 'bprime_tacot_earth.dat'),
                            os.path.join(out, 'INPUT', 'bprime_tacot_earth.dat'))
            for f, k in sorted(ENV_SCALE.items()):
                rows = []
                for l in open(os.path.join(PYRO_SRC, 'env_hi.dat')).read().split('\n'):
                    t = l.split()
                    if len(t) == 8 and not l.startswith('//'):
                        t[2] = '%.6g' % (float(t[2]) * k)    # rhoeUeCH0 only; h_r, chem, ... as written
                        l = ' '.join(t)
                    rows.append(l)
                open(os.path.join(out, 'INPUT', f), 'w').write('\n'.join(rows))
    if variant in ('receding', 'axis', 'openseam'):
        open(os.path.join(out, 'INPUT', 'env_r.dat'), 'w').write(ENV_R)
    ini += TAIL.format(iters=iters)
    open(os.path.join(out, 'input.ini'), 'w').write(ini)
    fuss_sh(out)
    shutil.rmtree(os.path.join(out, 'OUTPUT'), ignore_errors=True)


if __name__ == '__main__':
    a = sys.argv[1:]
    if a[:1] == ['--decks'] and len(a) <= 2:
        only = a[1].split(',') if len(a) == 2 else [v for v, _, _ in DECKS]
        for variant, kind, rel in DECKS:
            if variant not in only:
                continue
            d = os.path.join(ROOT, rel)
            for sub in ('INPUT', 'MESH'):
                shutil.rmtree(os.path.join(d, sub), ignore_errors=True)
            write(variant, kind, d)
            print('wrote', rel)
    elif len(a) == 3 and a[0] in MOTION and a[1] in ('sector', 'ring', 'sector1', 'ring1', 'sym'):
        write(a[0], a[1], a[2])
    else:
        sys.exit(__doc__)
