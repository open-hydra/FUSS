#!/usr/bin/env python3
"""
Energy conservation on a moving mesh  (plan 05, gate 5.4).

Adiabatic walls, no source, and interior-only mesh motion (the boundary nodes
are pinned by `taper-to-boundary`). Total energy

    E = sum_cells  h_i * V_i         h = volumetric enthalpy [J/m^3]

must therefore be identical at the start and the end of the run. Conduction
moves energy around; it does not create or destroy it, and neither should the
ALE swept-volume flux.

INDEPENDENT ORACLE
------------------
The cell volumes here are recomputed FROM THE NODE COORDINATES using the
divergence theorem, in this script, rather than read back from the solver. So a
mistake in the solver's volume would show up as an imbalance instead of
cancelling out on both sides of the comparison.

WHAT THIS CATCHES THAT GATE 5.1 CANNOT
--------------------------------------
5.1 uses a uniform field, where every choice of face value agrees. Only a
non-uniform field can distinguish them: the swept volume must carry the UPWIND
enthalpy, chosen on the sign of dV, so that the two cells sharing a face pick
the same value and their contributions cancel. A centred average passes 5.1 and
fails this test.
"""
import re
import sys

TOL_REL = 1.0e-10      # exact conservation is expected; this is a round-off budget

IC = "INPUT/ic.tec"
SOLUTION = "OUTPUT/field.tec"
PROPS = "INPUT/properties.dat"


def read_tec(path):
    with open(path) as fh:
        lines = fh.read().split("\n")
    # ic.tec writes  VARIABLES ="x" "y" "z" T matID   -- the trailing names are
    # NOT quoted, while field.tec quotes all of them. Handle both.
    header = lines[0].split("=", 1)[1] if "=" in lines[0] else lines[0]
    names = [tok.strip('"') for tok in header.split() if tok.strip('"')]
    m = re.search(r"I=(\d+),\s*J=(\d+),\s*K=(\d+)", lines[1])
    I, J, K = (int(x) for x in m.groups())
    vals = []
    for line in lines[2:]:
        s = line.strip()
        if not s:
            continue
        try:
            vals.append(float(s))
        except ValueError:
            pass
    n_nodal = I * J * K
    n_cell = (I - 1) * (J - 1) * max(K - 1, 1)
    data, off = {}, 0
    for name in names:
        n = n_nodal if name in ("x", "y", "z") else n_cell
        data[name] = vals[off:off + n]
        off += n
    return (I, J, K), data


def cell_volumes(dims, data):
    """Signed hexahedral volumes by the divergence theorem, from node coords."""
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

    vols = []
    for k in range(max(K - 1, 1)):
        for j in range(J - 1):
            for i in range(I - 1):
                p = []
                for a in (i, i + 1):
                    for b in (j, j + 1):
                        for c in (k, k + 1):
                            nd = nid(a, b, c)
                            p.append((X[nd], Y[nd], Z[nd]))
                # reorder to the (i,j,k) corner convention used by FACES
                idx = {(0, 0, 0): 0, (0, 0, 1): 1, (0, 1, 0): 2, (0, 1, 1): 3,
                       (1, 0, 0): 4, (1, 0, 1): 5, (1, 1, 0): 6, (1, 1, 1): 7}
                pts = [None] * 8
                t = 0
                for a in (0, 1):
                    for b in (0, 1):
                        for c in (0, 1):
                            pts[idx[(a, b, c)]] = p[t]
                            t += 1
                v = sum(quad_flux(pts[f[0]], pts[f[1]], pts[f[2]], pts[f[3]]) for f in FACES)
                vols.append(-v / 3.0)          # inward-ordered faces -> negate
    return vols


def h_of_T():
    """Interpolator for volumetric enthalpy h(T) from the property table."""
    Ts, hs = [], []
    with open(PROPS) as fh:
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


def main():
    try:
        dims0, d0 = read_tec(IC)
        dims1, d1 = read_tec(SOLUTION)
    except FileNotFoundError as e:
        print(f"ERROR: {e}")
        return 1

    hT = h_of_T()
    V0 = cell_volumes(dims0, d0)
    V1 = cell_volumes(dims1, d1)

    E0 = sum(hT(t) * v for t, v in zip(d0["T"], V0))
    # the solution file carries h directly; use it rather than round-tripping T
    h1 = d1["h"] if "h" in d1 else [hT(t) for t in d1["T"]]
    E1 = sum(hh * v for hh, v in zip(h1, V1))

    dE = E1 - E0
    rel = abs(dE) / abs(E0)

    print()
    print("------------------------------------------------------------")
    print(" ENERGY CONSERVATION ON A MOVING MESH")
    print(" (non-uniform field, adiabatic walls, interior mesh motion)")
    print("------------------------------------------------------------")
    print(f"  cells                      {len(V0)}")
    print(f"  initial T range            {min(d0['T']):.2f} .. {max(d0['T']):.2f} K")
    print(f"  final   T range            {min(d1['T']):.2f} .. {max(d1['T']):.2f} K")
    print(f"  total volume  start/end    {sum(V0):.12e} / {sum(V1):.12e} m^3")
    print(f"  total energy  start        {E0:.12e} J")
    print(f"  total energy  end          {E1:.12e} J")
    print(f"  absolute imbalance         {dE:.6e} J")
    print(f"  relative imbalance         {rel:.3e}")
    print(f"  tolerance                  {TOL_REL:.3e}")
    print("------------------------------------------------------------")

    ok = rel < TOL_REL
    if ok:
        print(f"  Result: PASS  (relative energy imbalance < {TOL_REL:.0e})")
    else:
        print("  Result: FAIL  -- the ALE flux is not conservative.")
        print("          The swept volume must carry the UPWIND enthalpy so the")
        print("          two cells sharing a face cancel; a centred or own-cell")
        print("          value breaks the telescoping sum.")
    print()
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
