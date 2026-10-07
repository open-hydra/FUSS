#!/usr/bin/env python3
"""
Free-stream preservation on a moving mesh  (plan 05, gate 5.1).

A uniform temperature field with adiabatic walls, no source, and an arbitrary
prescribed mesh motion must stay exactly uniform. Anything else means the ALE
convective flux and the cell-volume change disagree, i.e. the discrete
Geometric Conservation Law is violated.

WHY THE TOLERANCE IS WHAT IT IS
-------------------------------
The per-step GCL residual is at round-off, ~1e-15 relative to the cell volume
(measured directly by the MORPH unit tests). Over N steps the worst case is a
random walk of those residuals scaled by T0, so the drift budget is roughly

    N * 1e-15 * T0  ~=  300 * 1e-15 * 273.15  ~=  1e-10 K

The tolerance below is set an order of magnitude above that, and is still many
orders BELOW what an actual GCL violation would produce: getting the swept
volume wrong by even a per-cent of the volume change gives a drift of order
(dV/V)*T0 ~ 1e-3 K per step. So this test has a very wide margin between
"round-off" and "broken" -- it is not a tolerance tuned to pass.
"""
import re
import sys

T0 = 273.15          # the uniform initial temperature, from ICB-Block1
TOL_ABS = 1.0e-9     # K   (see the note above)

SOLUTION = "OUTPUT/field.tec"


def read_tec_cellcentred(path):
    """Return (variable_names, {name: [values]}) for a BLOCK-format Tecplot file."""
    with open(path) as fh:
        lines = fh.read().split("\n")

    names = [t.strip('"') for t in lines[0].split("=", 1)[-1].split() if t.strip('"')]  # VARIABLES names, quoted (ORION <= v1.6) or bare (ORION >= v1.7.0 writes `VARIABLES = x y z T ...`)
    m = re.search(r"I=(\d+),\s*J=(\d+),\s*K=(\d+)", lines[1])
    if not m:
        raise RuntimeError("could not parse zone dimensions")
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

    out, off = {}, 0
    for name in names:
        n = n_nodal if name in ("x", "y", "z") else n_cell
        out[name] = vals[off:off + n]
        off += n
    return names, out


def main():
    try:
        _, data = read_tec_cellcentred(SOLUTION)
    except FileNotFoundError:
        print(f"ERROR: file not found - {SOLUTION}")
        return 1

    T = data.get("T")
    if not T:
        print("ERROR: no T field in the solution")
        return 1

    dev = [abs(t - T0) for t in T]
    max_dev = max(dev)
    spread = max(T) - min(T)

    print()
    print("------------------------------------------------------------")
    print(" FREE-STREAM PRESERVATION ON A MOVING MESH")
    print(" (uniform field + arbitrary mesh motion + no forcing)")
    print("------------------------------------------------------------")
    print(f"  cells                        {len(T)}")
    print(f"  initial uniform temperature  {T0:.6f} K")
    print(f"  min / max temperature        {min(T):.12f} / {max(T):.12f} K")
    print(f"  max |T - T0|                 {max_dev:.3e} K")
    print(f"  max - min  (field spread)    {spread:.3e} K")
    print(f"  tolerance                    {TOL_ABS:.3e} K")
    print("------------------------------------------------------------")

    ok = max_dev < TOL_ABS and spread < TOL_ABS
    if ok:
        print(f"  Result: PASS  (max |T - T0| < {TOL_ABS:.0e} K)")
    else:
        print("  Result: FAIL  -- the mesh motion is injecting energy.")
        print("          A drift that grows with mesh velocity and does not")
        print("          shrink under refinement indicates a GCL violation:")
        print("          the swept volumes and the volume change disagree.")
    print()
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
