#!/usr/bin/env python3
"""
Energy conservation on a moving mesh, four blocks  (plan 05 gate 5.4 on the
2 x 2 split of ../conservation; plan 11 stage 2).

    usage: verify.py [case_dir]        (default: the current directory)

Run it from the case (or a copy of it) after the solve; it reads INPUT/ic.tec,
INPUT/properties.dat, INPUT/bc.txt and OUTPUT/field.tec there.

Two checks:

  * CONSERVATION. The physics check of ../conservation/verify.py on EVERY
    zone of INPUT/ic.tec and OUTPUT/field.tec (that script reads zone 1 only,
    and on this case it would see block 1 exchange heat with its neighbours:
    an imbalance of order 1e-7). Total energy E = sum h_i V_i over all blocks,
    the volumes recomputed here from the node coordinates by the divergence
    theorem, zone by zone with each zone's own dimensions, must be the same at
    the start and at the end: relative imbalance below 1e-10.

    `taper-to-boundary` tapers the motion to zero on the boundary of EACH
    block (Morph_Motion_Prescribed.f90: sin(pi i/dim) per block), so here the
    seams are pinned too and the interior of every block moves on its own.
    The motion therefore differs from the unsplit case's by construction and
    the two runs are not comparable; conservation must hold all the same,
    since every ALE flux is still internal and the 101 interfaces exchange
    conduction only.

  * SEAMS. Every 101 record of INPUT/bc.txt pairs a cell face with its
    neighbour block's; the four corner nodes of the two faces must be
    identical in field.tec. Compared as written (15 significant digits, so
    equal text is necessary, not proof of bit equality).
"""
import os
import re
import sys

TOL_REL = 1.0e-10      # exact conservation is expected; this is a round-off budget


def read_tec_zones(path):
    """Every zone of a BLOCK-format Tecplot file -> (names, [(dims, {name: [values]})]).
    ic.tec writes  VARIABLES ="x" "y" "z" T matID  (trailing names not quoted),
    field.tec  VARIABLES = x y z T ...; both are handled."""
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


def cell_volumes(dims, data):
    """Signed hexahedral volumes by the divergence theorem, from node coords
    (../conservation/verify.py, applied to one zone)."""
    I, J, K = dims
    X, Y, Z = data["x"], data["y"], data["z"]

    def nid(i, j, k):
        return i + j * I + k * I * J

    # faces of the hex, as corner indices into the 8-node list, inward-ordered
    FACES = ((0, 2, 3, 1), (4, 5, 7, 6), (0, 1, 5, 4),
             (2, 6, 7, 3), (0, 4, 6, 2), (1, 3, 7, 5))

    def tri_flux(p, q, r):
        u = [q[c] - p[c] for c in range(3)]
        v = [r[c] - p[c] for c in range(3)]
        n = (u[1] * v[2] - u[2] * v[1],
             u[2] * v[0] - u[0] * v[2],
             u[0] * v[1] - u[1] * v[0])
        cen = [(p[c] + q[c] + r[c]) / 3.0 for c in range(3)]
        return 0.5 * sum(cen[c] * n[c] for c in range(3))

    def quad_flux(a, b, c, d):
        m = [0.25 * (a[t] + b[t] + c[t] + d[t]) for t in range(3)]
        return tri_flux(a, b, m) + tri_flux(b, c, m) + tri_flux(c, d, m) + tri_flux(d, a, m)

    idx = {(0, 0, 0): 0, (0, 0, 1): 1, (0, 1, 0): 2, (0, 1, 1): 3,
           (1, 0, 0): 4, (1, 0, 1): 5, (1, 1, 0): 6, (1, 1, 1): 7}
    vols = []
    for k in range(max(K - 1, 1)):
        for j in range(J - 1):
            for i in range(I - 1):
                pts = [None] * 8
                for a in (0, 1):
                    for b in (0, 1):
                        for c in (0, 1):
                            nd = nid(i + a, j + b, k + c)
                            pts[idx[(a, b, c)]] = (X[nd], Y[nd], Z[nd])
                v = sum(quad_flux(pts[f[0]], pts[f[1]], pts[f[2]], pts[f[3]]) for f in FACES)
                vols.append(-v / 3.0)          # inward-ordered faces -> negate
    return vols


def h_of_T(props):
    """Interpolator for volumetric enthalpy h(T) from the property table."""
    Ts, hs = [], []
    with open(props) as fh:
        for line in fh:
            p = line.split()
            if len(p) == 5:
                try:
                    row = [float(x) for x in p]
                except ValueError:
                    continue
                Ts.append(row[0])
                hs.append(row[4])

    def interp(T):
        if T <= Ts[0]:
            return hs[0]
        if T >= Ts[-1]:
            return hs[-1]
        lo, hi = 0, len(Ts) - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if Ts[mid] <= T:
                lo = mid
            else:
                hi = mid
        w = (T - Ts[lo]) / (Ts[hi] - Ts[lo])
        return hs[lo] + w * (hs[hi] - hs[lo])

    return interp


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
        lines = fh.read().split("\n")
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
                    compared += 1
                    if A[c][la] != B[c][lb]:
                        bad += 1
                        worst = max(worst, abs(A[c][la] - B[c][lb]))
        n += 2 if ty in two else 1
    return compared, bad, worst


def main():
    case = sys.argv[1] if len(sys.argv) > 1 else "."
    ic = os.path.join(case, "INPUT", "ic.tec")
    solution = os.path.join(case, "OUTPUT", "field.tec")
    try:
        _, z0 = read_tec_zones(ic)
        _, z1 = read_tec_zones(solution)
    except (IOError, OSError) as e:
        print("ERROR: %s" % e)
        return 1
    if len(z0) != len(z1) or any(a[0] != b[0] for a, b in zip(z0, z1)):
        print("ERROR: %s and %s do not hold the same blocks" % (ic, solution))
        return 1

    hT = h_of_T(os.path.join(case, "INPUT", "properties.dat"))
    V0, V1, E0, E1, T0s, T1s = 0.0, 0.0, 0.0, 0.0, [], []
    for (dims0, d0), (dims1, d1) in zip(z0, z1):
        v0 = cell_volumes(dims0, d0)
        v1 = cell_volumes(dims1, d1)
        V0 += sum(v0)
        V1 += sum(v1)
        E0 += sum(hT(t) * v for t, v in zip(d0["T"], v0))
        # the solution file carries h directly; use it rather than round-tripping T
        h1 = d1["h"] if "h" in d1 else [hT(t) for t in d1["T"]]
        E1 += sum(hh * v for hh, v in zip(h1, v1))
        T0s += d0["T"]
        T1s += d1["T"]

    dE = E1 - E0
    rel = abs(dE) / abs(E0)
    compared, bad, worst = seam_check(z1, os.path.join(case, "INPUT", "bc.txt"))

    print()
    print("------------------------------------------------------------")
    print(" ENERGY CONSERVATION ON A MOVING MESH, %d BLOCKS" % len(z1))
    print(" (non-uniform field, adiabatic walls, interior mesh motion)")
    print("------------------------------------------------------------")
    print("  blocks / cells             %d / %d" % (len(z1), len(T1s)))
    print("  initial T range            %.2f .. %.2f K" % (min(T0s), max(T0s)))
    print("  final   T range            %.2f .. %.2f K" % (min(T1s), max(T1s)))
    print("  total volume  start/end    %.12e / %.12e m^3" % (V0, V1))
    print("  total energy  start        %.12e J" % E0)
    print("  total energy  end          %.12e J" % E1)
    print("  absolute imbalance         %.6e J" % dE)
    print("  relative imbalance         %.3e" % rel)
    print("  tolerance                  %.3e" % TOL_REL)
    print("  seam node coordinates      %d compared, %d differ (max %.3e m)" % (compared, bad, worst))
    print("------------------------------------------------------------")

    ok_e = rel < TOL_REL
    ok_seam = compared > 0 and bad == 0
    if ok_e:
        print("  conservation: PASS  (relative energy imbalance < %.0e)" % TOL_REL)
    else:
        print("  conservation: FAIL  -- the ALE flux (or a 101 interface) is not conservative.")
    if ok_seam:
        print("  seams:        PASS  (every node shared by two blocks is written alike)")
    elif compared == 0:
        print("  seams:        FAIL  -- INPUT/bc.txt holds no 101 record: not a split case.")
    else:
        print("  seams:        FAIL  -- blocks disagree on the position of a shared node.")
    ok = ok_e and ok_seam
    print("  Result: %s" % ("PASS" if ok else "FAIL"))
    print()
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
