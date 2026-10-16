#!/usr/bin/env python3
"""
Free-stream preservation on a moving mesh, four blocks  (plan 05 gate 5.1 on the
2 x 2 split of ../freestream; plan 11 stage 2).

    usage: verify.py [case_dir]        (default: the current directory)

Run it from the case (or a copy of it) after the solve; it reads
OUTPUT/field.tec and INPUT/bc.txt there.

Two checks:

  * FREE STREAM. The physics check of ../freestream/verify.py, unchanged, on
    EVERY zone of OUTPUT/field.tec (that script reads zone 1 only): a uniform
    T0 with adiabatic walls, no source and an arbitrary prescribed motion must
    stay uniform, max |T - T0| and max - min below 1e-9 K. The budget is
    ../freestream/verify.py's: about 1e-10 K of round-off random walk over 300
    steps, an order of magnitude of margin, many orders below what a GCL
    violation gives.

  * SEAMS. Every 101 record of INPUT/bc.txt pairs a cell face with its
    neighbour block's; the four corner nodes of the two faces must be
    identical in field.tec. The prescribed law is analytic in the reference
    node position (Morph_Motion_Prescribed.f90), so both blocks move a shared
    node to the same place; a block whose nodes the writer did not refresh
    (a remote block under MPI) shows up here. Compared as written: field.tec
    carries 15 significant digits, so two doubles one ulp apart can print
    alike -- equal text is necessary, not proof of bit equality.
"""
import os
import re
import sys

T0 = 273.15          # the uniform initial temperature (INPUT/ic.tec)
TOL_ABS = 1.0e-9     # K   (see ../freestream/verify.py)


def read_tec_zones(path):
    """Every zone of a BLOCK-format Tecplot file -> (names, [(dims, {name: [values]})])."""
    with open(path) as fh:
        lines = fh.read().split("\n")
    names = [t.strip('"') for t in lines[0].split("=", 1)[-1].split() if t.strip('"')]
    starts = [n for n, line in enumerate(lines) if line.lstrip().startswith("ZONE")]
    if not starts:
        raise RuntimeError("no ZONE in %s" % path)
    zones = []
    for z, first in enumerate(starts):
        m = re.search(r"I=(\d+),\s*J=(\d+),\s*K=(\d+)", lines[first])
        if not m:
            raise RuntimeError("could not parse the dimensions of zone %d of %s" % (z + 1, path))
        I, J, K = (int(x) for x in m.groups())
        last = starts[z + 1] if z + 1 < len(starts) else len(lines)
        vals = [float(line) for line in lines[first + 1:last] if line.strip()]
        n_nodal = I * J * K
        n_cell = (I - 1) * (J - 1) * max(K - 1, 1)
        data, off = {}, 0
        for name in names:
            n = n_nodal if name in ("x", "y", "z") else n_cell
            data[name] = vals[off:off + n]
            off += n
        if off != len(vals):
            raise RuntimeError("zone %d of %s: %d values, %d expected" % (z + 1, path, len(vals), off))
        zones.append(((I, J, K), data))
    return names, zones


def face_nodes(i, j, k, f):
    """The four corner nodes (0-based) of face f of cell (i, j, k) (1-based),
    ordered by (first, second) tangential offset."""
    out = []
    for a in (0, 1):
        for c in (0, 1):
            if f in (1, 2):
                out.append((i - 1 if f == 1 else i, j - 1 + a, k - 1 + c))
            elif f in (3, 4):
                out.append((i - 1 + a, j - 1 if f == 3 else j, k - 1 + c))
            else:
                out.append((i - 1 + a, j - 1 + c, k - 1 if f == 5 else k))
    return out


def seam_check(zones, bc_path):
    """(pairs of nodes compared, mismatches, max |difference|) over every 101 record."""
    with open(bc_path) as fh:
        lines = [line for line in fh.read().split("\n")]
    two = (101, 103, 301, 302, 303, 304, 305, 306, 307)
    n, compared, bad, worst = 0, 0, 0, 0.0
    while n < len(lines) and lines[n].strip():
        b, i, j, k, f, ty = (int(t) for t in lines[n].split())
        if ty == 101:
            t = [int(x) for x in lines[n + 1].split()]
            bs, i_s, js, ks, fs = t[:5]
            if t[5:] != [1, 0, 0, 1]:
                raise RuntimeError("bc.txt record at line %d: orientation %s, only 1 0 0 1 is checked"
                                   % (n + 1, t[5:]))
            (I, J, K), A = zones[b - 1]
            (Is, Js, Ks), B = zones[bs - 1]
            for (p, q) in zip(face_nodes(i, j, k, f), face_nodes(i_s, js, ks, fs)):
                la = p[0] + I * (p[1] + J * p[2])
                lb = q[0] + Is * (q[1] + Js * q[2])
                for c in ("x", "y", "z"):
                    d = abs(A[c][la] - B[c][lb])
                    compared += 1
                    if A[c][la] != B[c][lb]:
                        bad += 1
                        worst = max(worst, d)
        n += 2 if ty in two else 1
    return compared, bad, worst


def main():
    case = sys.argv[1] if len(sys.argv) > 1 else "."
    solution = os.path.join(case, "OUTPUT", "field.tec")
    try:
        _, zones = read_tec_zones(solution)
    except (IOError, OSError):
        print("ERROR: file not found - %s" % solution)
        return 1

    T = [t for _, data in zones for t in data.get("T", [])]
    if not T:
        print("ERROR: no T field in the solution")
        return 1

    max_dev = max(abs(t - T0) for t in T)
    spread = max(T) - min(T)
    compared, bad, worst = seam_check(zones, os.path.join(case, "INPUT", "bc.txt"))

    print()
    print("------------------------------------------------------------")
    print(" FREE-STREAM PRESERVATION ON A MOVING MESH, %d BLOCKS" % len(zones))
    print(" (uniform field + arbitrary mesh motion + no forcing)")
    print("------------------------------------------------------------")
    print("  blocks / cells               %d / %d" % (len(zones), len(T)))
    print("  initial uniform temperature  %.6f K" % T0)
    print("  min / max temperature        %.12f / %.12f K" % (min(T), max(T)))
    print("  max |T - T0|                 %.3e K" % max_dev)
    print("  max - min  (field spread)    %.3e K" % spread)
    print("  tolerance                    %.3e K" % TOL_ABS)
    print("  seam node coordinates        %d compared, %d differ (max %.3e m)" % (compared, bad, worst))
    print("------------------------------------------------------------")

    ok_fs = max_dev < TOL_ABS and spread < TOL_ABS
    ok_seam = compared > 0 and bad == 0
    if ok_fs:
        print("  free stream: PASS  (max |T - T0| < %.0e K)" % TOL_ABS)
    else:
        print("  free stream: FAIL  -- the mesh motion is injecting energy (GCL violation).")
    if ok_seam:
        print("  seams:       PASS  (every node shared by two blocks is written alike)")
    elif compared == 0:
        print("  seams:       FAIL  -- INPUT/bc.txt holds no 101 record: not a split case.")
    else:
        print("  seams:       FAIL  -- blocks disagree on the position of a shared node.")
    ok = ok_fs and ok_seam
    print("  Result: %s" % ("PASS" if ok else "FAIL"))
    print()
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
