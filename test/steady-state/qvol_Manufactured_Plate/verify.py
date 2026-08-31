#!/usr/bin/env python3
"""
2-D steady-state conduction with a NON-UNIFORM volumetric source.

Method of manufactured solutions. For constant conductivity k the steady heat
equation is

    k d2T/dx2 + k d2T/dy2 + q(x,y) = 0

Choosing

    T_exact(x,y) = T0 + A sin(pi x / Lx) sin(pi y / Ly)

gives, by direct substitution,

    q(x,y) = k A [ (pi/Lx)^2 + (pi/Ly)^2 ] sin(pi x / Lx) sin(pi y / Ly)

and T_exact = T0 on all four edges, which matches the fixed-temperature
boundary condition used by this case. INPUT/st.tec carries exactly that q.

WHY THIS CASE EXISTS
--------------------
test/steady-state/qvol_Plate is the only other case with an INPUT/st.tec, and
its source field is UNIFORM (1e6 W/m^3, a single distinct value over all 2400
cells). A uniform field cannot detect an indexing error in the source reader,
because every cell holds the same number. This case uses a source with 2398
distinct values over 2400 cells and compares the whole field, so a misread of
INPUT/st.tec shows up directly.

The tolerance is a discretisation-error budget, not a fitted value: the scheme
is second order on a 60x40 uniform mesh, and the observed error on the correct
code is well inside TOL_LINF.
"""

import os
import re
import sys
import math

import numpy as np

FIELD_TEC   = os.path.join(os.path.dirname(__file__), "OUTPUT", "field.tec")
TARGET_ZONE = "Block1"

# Manufactured-solution parameters -- must match INPUT/st.tec and input.ini
LX, LY = 0.6, 0.4          # plate dimensions      [m]
T0     = 273.15            # edge temperature      [K]  (= Twall in input.ini)
AMP    = 100.0             # peak temperature rise [K]

# Acceptance thresholds on the temperature error over the whole field.
TOL_LINF = 1.0             # [K]  = 1 % of the 100 K rise
TOL_L2   = 0.5             # [K]
# The source read-back is an exact identity, not a discretised quantity: the
# solver should hand back the same numbers it read from INPUT/st.tec. The
# threshold is a round-trip ASCII-formatting allowance only.
TOL_QVOL = 1.0             # [W/m3] against a peak of ~4.6e5


def parse_tec(filepath):
    """Minimal Tecplot BLOCK reader: returns (variables, {zone: {header, raw}})."""
    with open(filepath) as fh:
        lines = fh.readlines()

    variables, zones, header, raw = [], {}, None, []
    for line in lines:
        s = line.strip()
        if re.match(r"(?i)variables", s):
            variables = re.findall(r'"([^"]+)"', s)
            # unquoted trailing names (the writer does not always quote them)
            tail = re.sub(r'(?i)^variables\s*=', '', s)
            for tok in re.sub(r'"[^"]*"', ' ', tail).replace(',', ' ').split():
                if tok not in variables:
                    variables.append(tok)
            continue
        if re.match(r"(?i)zone", s):
            if header is not None:
                zones[header["name"]] = {"header": header, "raw": raw}
            raw = []
            name = re.search(r'(?i)\bT\s*=\s*([^,]+)', s)
            header = {
                "name": name.group(1).strip() if name else "zone%d" % len(zones),
                "I": int(re.search(r'(?i)\bI\s*=\s*(\d+)', s).group(1)),
                "J": int(re.search(r'(?i)\bJ\s*=\s*(\d+)', s).group(1)),
                "K": int(re.search(r'(?i)\bK\s*=\s*(\d+)', s).group(1)),
            }
            continue
        if header is not None and s:
            raw.extend(float(v) for v in s.split())
    if header is not None:
        zones[header["name"]] = {"header": header, "raw": raw}
    return variables, zones


def main():
    if not os.path.isfile(FIELD_TEC):
        print("ERROR: %s not found -- run './FUSS.sh solve' first." % FIELD_TEC)
        sys.exit(1)

    variables, zones = parse_tec(FIELD_TEC)

    zone = None
    for name, z in zones.items():
        if TARGET_ZONE in name or len(zones) == 1:
            zone = z
            break
    if zone is None:
        print("ERROR: zone '%s' not found. Available: %s"
              % (TARGET_ZONE, ", ".join(zones)))
        sys.exit(1)

    hdr = zone["header"]
    I, J, K = hdr["I"], hdr["J"], hdr["K"]
    nn = I * J * K                          # nodal values per variable
    nc = (I - 1) * (J - 1) * max(K - 1, 1)  # cell-centred values per variable

    raw = np.asarray(zone["raw"], dtype=float)
    x, y = raw[0:nn], raw[nn:2 * nn]

    # First cell-centred variable after the three nodal coordinate blocks is T.
    off = 3 * nn
    if len(raw) < off + nc:
        print("ERROR: zone too short: %d values, expected at least %d"
              % (len(raw), off + nc))
        sys.exit(1)
    T = raw[off:off + nc]

    # Cell centres from the nodal coordinates (i fastest, then j, then k).
    x3 = x.reshape(K, J, I)
    y3 = y.reshape(K, J, I)
    xc = 0.25 * (x3[0, :-1, :-1] + x3[0, :-1, 1:] + x3[0, 1:, :-1] + x3[0, 1:, 1:])
    yc = 0.25 * (y3[0, :-1, :-1] + y3[0, :-1, 1:] + y3[0, 1:, :-1] + y3[0, 1:, 1:])

    T2 = T.reshape(J - 1, I - 1)

    # ------------------------------------------------------------------
    # Direct check of the source reader.
    # The solver echoes the qvol it actually used back into field.tec, so we
    # can compare it against the analytic q that INPUT/st.tec was built from.
    # This isolates "did the source get read correctly" from "did the solve
    # converge", and is what fails loudly if the st.tec index is wrong.
    # ------------------------------------------------------------------
    q_err_linf = None
    if "qvol" in variables:
        iq = variables.index("qvol")          # position among ALL variables
        ncc = iq - 3                          # index among cell-centred ones
        qs = off + ncc * nc
        if len(raw) >= qs + nc:
            q_used = raw[qs:qs + nc].reshape(J - 1, I - 1)
            coef = 52.0 * AMP * ((np.pi / LX) ** 2 + (np.pi / LY) ** 2)
            q_exact = coef * np.sin(np.pi * xc / LX) * np.sin(np.pi * yc / LY)
            q_err_linf = float(np.max(np.abs(q_used - q_exact)))

    T_exact = T0 + AMP * np.sin(np.pi * xc / LX) * np.sin(np.pi * yc / LY)
    err = T2 - T_exact

    linf = float(np.max(np.abs(err)))
    l2 = float(np.sqrt(np.mean(err ** 2)))
    imax, jmax = np.unravel_index(np.argmax(np.abs(err)), err.shape)

    print()
    print("------------------------------------------------------------")
    print(" 2-D STEADY STATE, NON-UNIFORM VOLUMETRIC SOURCE")
    print(" Method of manufactured solutions")
    print("------------------------------------------------------------")
    print("  Cells                       %d x %d" % (I - 1, J - 1))
    print("  Peak T (exact / computed)   %8.3f / %8.3f K"
          % (float(np.max(T_exact)), float(np.max(T2))))
    print("  L-inf error                 %8.4f K   (tol %.3f)" % (linf, TOL_LINF))
    print("  L-2   error                 %8.4f K   (tol %.3f)" % (l2, TOL_L2))
    print("  worst cell (i,j)            (%d, %d)  at x=%.4f y=%.4f"
          % (jmax, imax, xc[imax, jmax], yc[imax, jmax]))
    if q_err_linf is not None:
        print("  source read-back L-inf      %8.3f W/m3  (tol %.3f)"
              % (q_err_linf, TOL_QVOL))
    print("------------------------------------------------------------")

    ok = (linf <= TOL_LINF) and (l2 <= TOL_L2)
    if q_err_linf is not None:
        ok = ok and (q_err_linf <= TOL_QVOL)
    print()
    if ok:
        print("  Result: PASS  (L-inf < %.3f K, L-2 < %.3f K)" % (TOL_LINF, TOL_L2))
    else:
        print("  Result: FAIL")
        print("  A misread of INPUT/st.tec shifts the source field and shows up")
        print("  here as a systematic error plus large errors near the i-max edge.")
    print()
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
