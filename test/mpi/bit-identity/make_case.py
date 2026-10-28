#!/usr/bin/env python3
"""Inputs of two tests of the parallel layer, from one lattice of cells sheared by x += S y (so that the tangential
terms of the connection and chimera fluxes are not zero).

  make_case.py bit  [--outdir DIR]   MPIBitIdentity: the same output bits on any layout of ranks x threads
  make_case.py edge [--outdir DIR]   ChimeraEdge: linear field on the chimera pair, walls consistent with it

Both write DIR/MESH/mesh.tec and DIR/INPUT/{ic.tec,bc.txt,phase.txt,properties.dat} (3-D, NK cells in z, DX = DY = DZ).

Blocks (cell lattice, block cell (i,j,k) is at lattice cell (X0+i, Y0+j, Z0+k)):
  1-4  (bit only) a 2 x 2 box of 10 x 10 cells around a vertex the four blocks share, block b = 1 + p + 2 q (p along x,
       q along y), 101 connections between neighbours (identity orientation); outer faces 302 (x-min), 301 (x-max),
       303 (y-min), 305 (y-max), symmetry 300 in z;
  5, 6 a chimera pair (blocks 1 and 2 in `edge`): block 5 (16 x 16 cells) and block 6 (15 x 15) overlap by 3 x 3 cells at the top-right corner of
       block 5 (the bottom-left corner of block 6), on the same lattice: every ghost cell of a chimera face coincides with
       one donor cell (volume fraction 1, no ATLAS needed). Faces 2 and 4 of block 5 and faces 1 and 3 of block 6 are
       chimera where they overlap and walls elsewhere: two chimera faces meet at an edge of each block, and the chimera
       faces cover part of a block face. In `bit`, walls 302 / 301 / 303 / 305 with values that depend on the record,
       symmetry in z; in `edge`, only blocks 5 and 6, every outer face 302 with the wall temperature of the linear field at
       the face centre (z faces too).
  The block sizes are chosen so that the LPT split of the blocks over 2 and 4 ranks puts blocks 5 and 6 on different ranks.
"""
import argparse
import math
import os

NK = 3
DX = 0.01                      # m
SHEAR = 0.3
K, RHO, CP = 52.0, 7850.0, 460.0
TMAX = 3000
LIN = (1000.0, 800.0, 1200.0, 2000.0)   # T = A + B x + C y + D z [K, K/m]; edge test

# id, X0, Y0, ni, nj
BOX = [(1 + p + 2 * q, 10 * p, 10 * q, 10, 10) for q in range(2) for p in range(2)]
PAIR = [(5, 25, 0, 16, 16), (6, 38, 13, 15, 15)]


def node(X, Y, Z):
    """Physical position of lattice node (X, Y, Z): the lattice sheared in x."""
    return ((X + SHEAR * Y) * DX, Y * DX, Z * DX)


def centre(blk, i, j, k):
    """Centre of cell (i, j, k) of a block."""
    _, X0, Y0, _, _ = blk
    pts = [node(X0 + i - 1 + a, Y0 + j - 1 + b, k - 1 + c) for a in (0, 1) for b in (0, 1) for c in (0, 1)]
    return tuple(sum(p[d] for p in pts) / 8.0 for d in range(3))


def face_centre(blk, i, j, k, f):
    """Centre of face f of cell (i, j, k)."""
    _, X0, Y0, _, _ = blk
    lo, hi = {1: (0, 0), 2: (1, 1), 3: (0, 0), 4: (1, 1), 5: (0, 0), 6: (1, 1)}[f]
    ax = (f - 1) // 2
    pts = []
    for a in (0, 1):
        for b in (0, 1):
            off = [0, 0, 0]
            others = [d for d in range(3) if d != ax]
            off[others[0]], off[others[1]] = a, b
            off[ax] = hi
            pts.append(node(X0 + i - 1 + off[0], Y0 + j - 1 + off[1], k - 1 + off[2]))
    return tuple(sum(p[d] for p in pts) / 4.0 for d in range(3))


def e(v):
    return f' {v:.15E}'


def t_init(mode, x, y, z):
    if mode == "edge":
        return LIN[0] + LIN[1] * x + LIN[2] * y + LIN[3] * z
    # bit: a smooth field varying in x, y and z
    return (1500.0 + 300.0 * math.sin(40.0 * x) * math.cos(30.0 * y) + 200.0 * math.sin(25.0 * x + 35.0 * y)
            + 150.0 * math.cos(60.0 * z + 20.0 * x))


def write_ic(path, blocks, mode):
    with open(path, "w") as f:
        f.write(' VARIABLES ="x" "y" "z" T matID\n')
        for blk in blocks:
            b, X0, Y0, ni, nj = blk
            xs = [node(X0 + i, Y0 + j, k) for k in range(NK + 1) for j in range(nj + 1) for i in range(ni + 1)]
            f.write(f' ZONE  T = B{b}-SP, I={ni + 1}, J={nj + 1}, K={NK + 1}, DATAPACKING=BLOCK, '
                    'VARLOCATION=([1-3]=NODAL,[4-5]=CELLCENTERED), SOLUTIONTIME=0.0, STRANDID = 0\n')
            vals = [p[d] for d in range(3) for p in xs]
            vals += [t_init(mode, *centre(blk, i, j, k)) for k in range(1, NK + 1) for j in range(1, nj + 1)
                     for i in range(1, ni + 1)]
            vals += [1.0] * (ni * nj * NK)
            f.write('\n'.join(e(v) for v in vals) + '\n')


def write_mesh(path, blocks):
    with open(path, "w") as f:
        f.write(' VARIABLES ="x" "y" "z"\n')
        for blk in blocks:
            b, X0, Y0, ni, nj = blk
            xs = [node(X0 + i, Y0 + j, k) for k in range(NK + 1) for j in range(nj + 1) for i in range(ni + 1)]
            f.write(f' ZONE  T = Block{b}, I={ni + 1}, J={nj + 1}, K={NK + 1}, DATAPACKING=BLOCK, '
                    'VARLOCATION=([1-3]=NODAL), SOLUTIONTIME=0.0, STRANDID = 0\n')
            f.write('\n'.join(e(p[d]) for d in range(3) for p in xs) + '\n')


GUIDE = {1: (-1, 0, 0), 2: (1, 0, 0), 3: (0, -1, 0), 4: (0, 1, 0), 5: (0, 0, -1), 6: (0, 0, 1)}   # outward normal


def donor(blocks, blk, i, j, k, f, g):
    """Donor cell (block, i, j, k) of the g-th ghost cell of boundary cell (i,j,k) of face f of blk, or None: the cell of
    another block of the lattice that contains it (both ghost layers are required by the caller)."""
    b, X0, Y0, _, _ = blk
    n = GUIDE[f]
    X, Y, Z = X0 + i + g * n[0], Y0 + j + g * n[1], k + g * n[2]
    for o in blocks:
        ob, oX0, oY0, oni, onj = o
        if ob == b:
            continue
        oi, oj, ok = X - oX0, Y - oY0, Z
        if 1 <= oi <= oni and 1 <= oj <= onj and 1 <= ok <= NK:
            return (ob, oi, oj, ok)
    return None


def wall_lines(mode, blk, i, j, k, f):
    b = blk[0]
    if mode == "edge":
        x, y, z = face_centre(blk, i, j, k, f)
        return 302, [t_init(mode, x, y, z)]
    if f in (5, 6):
        return 300, []
    s = 1.0 + 0.01 * (i + 2 * j + 3 * k) + 0.1 * b
    return {1: (302, [900.0 * s]), 2: (301, [2.0e4 * s]), 3: (303, [150.0 * s, 1.0e3 * s, 1100.0 * s]),
            4: (305, [120.0 * s, 0.7 + 0.01 * (i + j) / 20.0, 1200.0 * s])}[f]


def connection(blocks_box, blk, i, j, k, f):
    """(block, i, j, k, face) of the cell across face f of boundary cell (i,j,k) of a box block, or None."""
    b, X0, Y0, ni, nj = blk
    p, q = (b - 1) % 2, (b - 1) // 2
    if f == 2 and p == 0:
        return (b + 1, 1, j, k, 1)
    if f == 1 and p == 1:
        return (b - 1, 10, j, k, 2)
    if f == 4 and q == 0:
        return (b + 2, i, 1, k, 3)
    if f == 3 and q == 1:
        return (b - 2, i, 10, k, 4)
    return None


def write_bc(path, blocks, mode):
    out = []
    for blk in blocks:
        b, X0, Y0, ni, nj = blk
        dims = (ni, nj, NK)
        for f in range(1, 7):
            ax = (f - 1) // 2
            t = [d for d in range(3) if d != ax]
            for c2 in range(1, dims[t[1]] + 1):
                for c1 in range(1, dims[t[0]] + 1):
                    ijk = [0, 0, 0]
                    ijk[t[0]], ijk[t[1]] = c1, c2
                    ijk[ax] = 1 if f % 2 == 1 else dims[ax]
                    i, j, k = ijk
                    head = ''.join(f'{v:8d}' for v in (b, i, j, k, f))
                    if blk in BOX:
                        cn = connection(None, blk, i, j, k, f)
                        if cn is not None:
                            out += [head + f'{101:8d}', ''.join(f'{v:8d}' for v in cn) + f'{1:8d}{0:8d}{0:8d}{1:8d}']
                            continue
                    elif f <= 4:
                        d1, d2 = donor(blocks, blk, i, j, k, f, 1), donor(blocks, blk, i, j, k, f, 2)
                        if d1 is not None and d2 is not None:
                            out += [head + f'{102:8d}', f'{1:8d}{1:8d}',
                                    ''.join(f'{v:8d}' for v in d1) + '    ' + e(1.0),
                                    ''.join(f'{v:8d}' for v in d2) + '    ' + e(1.0)]
                            continue
                    typ, vals = wall_lines(mode, blk, i, j, k, f)
                    out.append(head + f'{typ:8d}')
                    if vals:
                        out.append(''.join(f'    {v:.15E},' for v in vals))
    open(path, "w").write('\n'.join(out) + '\n')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["bit", "edge"])
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()
    blocks = (BOX + PAIR) if a.mode == "bit" else [(1,) + PAIR[0][1:], (2,) + PAIR[1][1:]]   # block numbers 1..nb
    d = os.path.join(a.outdir, "INPUT")
    os.makedirs(d, exist_ok=True)
    os.makedirs(os.path.join(a.outdir, "MESH"), exist_ok=True)
    write_mesh(os.path.join(a.outdir, "MESH", "mesh.tec"), blocks)
    write_ic(os.path.join(d, "ic.tec"), blocks, a.mode)
    write_bc(os.path.join(d, "bc.txt"), blocks, a.mode)
    open(os.path.join(d, "phase.txt"), "w").write("solid-bulk\nmat 1\n")
    with open(os.path.join(d, "properties.dat"), "w") as f:
        f.write('TITLE = "Mass Thermodynamic Properties"\n'
                'VARIABLES = "Temperature", "Cp", "Density", "Conductivity", "Energy"\n'
                f'ZONE T="mat"\nI={TMAX}, F=POINT\n')
        f.writelines(f'{t:.1f} {CP:.6f} {RHO:.6f} {K:.6f} {RHO * CP * t:.6f}\n' for t in range(1, TMAX + 1))


if __name__ == "__main__":
    main()
