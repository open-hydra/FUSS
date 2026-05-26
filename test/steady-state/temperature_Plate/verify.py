import os
import re
import sys

import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

FIELD_TEC    = os.path.join(os.path.dirname(__file__), "OUTPUT", "field.tec")
SVG_OUT      = os.path.join(os.path.dirname(__file__), "2D-steady-state-temperature-plate-field.svg")
TARGET_ZONE  = "Block1"
TARGET_X     = 0.2
TARGET_Y     = 0.2
REFERENCE_T  = 273.15 + 260.5 # K


def parse_tec(filepath):

    with open(filepath, "r") as fh:
        lines = fh.readlines()

    variables = []
    zones = {}
    current_header = None
    raw = []

    for line in lines:
        stripped = line.strip()

        if re.match(r"(?i)variables", stripped):
            variables = re.findall(r'"([^"]+)"', stripped)
            continue

        if re.match(r"(?i)zone", stripped):
            if current_header is not None:
                zones[current_header["name"]] = {"header": current_header, "raw": raw}
            raw = []

            name_m = re.search(r'T\s*=\s*(\w+)',   stripped, re.IGNORECASE)
            i_m    = re.search(r'\bI\s*=\s*(\d+)', stripped, re.IGNORECASE)
            j_m    = re.search(r'\bJ\s*=\s*(\d+)', stripped, re.IGNORECASE)
            k_m    = re.search(r'\bK\s*=\s*(\d+)', stripped, re.IGNORECASE)
            cc_m   = re.search(r'\[(\d+)-(\d+)\]\s*=\s*CELLCENTERED',
                               stripped, re.IGNORECASE)

            current_header = {
                "name"    : name_m.group(1) if name_m else f"Zone{len(zones)+1}",
                "ni"      : int(i_m.group(1)) if i_m else 1,
                "nj"      : int(j_m.group(1)) if j_m else 1,
                "nk"      : int(k_m.group(1)) if k_m else 1,
                "cc_start": int(cc_m.group(1)) if cc_m else 0,
                "cc_end"  : int(cc_m.group(2)) if cc_m else 0,
            }
            continue

        try:
            raw.append(float(stripped))
        except ValueError:
            pass

    if current_header is not None:
        zones[current_header["name"]] = {"header": current_header, "raw": raw}

    return variables, zones


def extract_vars(variables, zone_info):

    hdr = zone_info["header"]
    raw = zone_info["raw"]
    ni, nj, nk       = hdr["ni"], hdr["nj"], hdr["nk"]
    cc_start, cc_end = hdr["cc_start"], hdr["cc_end"]

    n_nodal = ni * nj * nk
    n_cell  = max(ni-1,1) * max(nj-1,1) * max(nk-1,1)

    result = {}
    idx = 0
    for v_idx, vname in enumerate(variables, start=1):
        n = n_cell if (cc_start <= v_idx <= cc_end) else n_nodal
        result[vname] = np.array(raw[idx : idx + n])
        idx += n

    return result


def build_cell_centres(x_nodes, y_nodes, ni_nodes, nj_nodes):
    ni_cells = ni_nodes - 1
    nj_cells = nj_nodes - 1

    x2d = x_nodes[: ni_nodes * nj_nodes].reshape(nj_nodes, ni_nodes)
    y2d = y_nodes[: ni_nodes * nj_nodes].reshape(nj_nodes, ni_nodes)

    x_cc = 0.25 * (x2d[:nj_cells, :ni_cells] + x2d[:nj_cells, 1:]
                 + x2d[1:,        :ni_cells] + x2d[1:,        1:])
    y_cc = 0.25 * (y2d[:nj_cells, :ni_cells] + y2d[:nj_cells, 1:]
                 + y2d[1:,        :ni_cells] + y2d[1:,        1:])

    return x_cc, y_cc


def bilinear_interp(x_cc, y_cc, T_2d, x_q, y_q):
    x_u = np.unique(x_cc[0, :])
    y_u = np.unique(y_cc[:, 0])

    ix = np.searchsorted(x_u, x_q) - 1
    iy = np.searchsorted(y_u, y_q) - 1
    ix = max(0, min(ix, len(x_u) - 2))
    iy = max(0, min(iy, len(y_u) - 2))

    x0, x1 = x_u[ix], x_u[ix + 1]
    y0, y1 = y_u[iy], y_u[iy + 1]

    tx = (x_q - x0) / (x1 - x0)
    ty = (y_q - y0) / (y1 - y0)

    T00 = T_2d[iy,     ix    ]
    T10 = T_2d[iy,     ix + 1]
    T01 = T_2d[iy + 1, ix    ]
    T11 = T_2d[iy + 1, ix + 1]

    T_interp = ((1 - tx) * (1 - ty) * T00
              +      tx  * (1 - ty) * T10
              + (1 - tx) *      ty  * T01
              +      tx  *      ty  * T11)

    cells = [
        {"ic": ix,   "jc": iy,   "x": x0, "y": y0, "T": T00, "w": (1-tx)*(1-ty)},
        {"ic": ix+1, "jc": iy,   "x": x1, "y": y0, "T": T10, "w":    tx *(1-ty)},
        {"ic": ix,   "jc": iy+1, "x": x0, "y": y1, "T": T01, "w": (1-tx)*   ty },
        {"ic": ix+1, "jc": iy+1, "x": x1, "y": y1, "T": T11, "w":    tx *   ty },
    ]

    return float(T_interp), cells


def save_colormap(x_cc, y_cc, T_2d, point_A, svg_path):
    """Save a filled-contour temperature colourmap to an SVG file."""
    fig, ax = plt.subplots(figsize=(7.5, 5.5))

    T_celsius = T_2d - 273.15

    x_1d = x_cc[0, :]
    y_1d = y_cc[:, 0]

    pcm = ax.pcolormesh(
        x_1d, y_1d, T_celsius,
        cmap="RdBu_r",
        shading="nearest",
        vmin = 100.0,
        vmax = 1000.0
    )

    cbar = fig.colorbar(pcm, ax=ax, pad=0.02)
    cbar.set_label("Temperature (°C)", fontsize=11)
    cbar.ax.tick_params(labelsize=9)

    # Boundary outline
    nx, ny = x_cc.shape[1], x_cc.shape[0]
    xmin, xmax = float(x_cc[0, 0] - (x_cc[0, 1]-x_cc[0, 0])/2), \
                 float(x_cc[0,-1] + (x_cc[0,-1]-x_cc[0,-2])/2)
    ymin, ymax = float(y_cc[0, 0] - (y_cc[1, 0]-y_cc[0, 0])/2), \
                 float(y_cc[-1, 0] + (y_cc[-1, 0]-y_cc[-2, 0])/2)
    ax.plot([xmin, xmax, xmax, xmin, xmin],
            [ymin, ymin, ymax, ymax, ymin],
            "k-", linewidth=1.2)

    # Point A marker
    ax.plot(point_A[0], point_A[1], "o",
            color="cyan", markersize=7, markeredgecolor="black",
            markeredgewidth=0.8,
            label=f"P ({point_A[0]}, {point_A[1]}) m",
            zorder=5)
    ax.legend(fontsize=9, loc="upper right")

    ax.set_xlabel("x (m)", fontsize=11)
    ax.set_ylabel("y (m)", fontsize=11)
    ax.set_title(
        "Steady State Temperature Distribution\n",
        fontsize=11,
        fontweight="bold"
    )
    ax.set_aspect("equal")
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))
    ax.tick_params(labelsize=9)

    fig.tight_layout()
    fig.savefig(svg_path, format="svg", bbox_inches="tight")
    plt.close(fig)
    print(f"  Colourmap saved → {svg_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():

    if not os.path.isfile(FIELD_TEC):
        print(f"ERROR: file not found – {FIELD_TEC}")
        sys.exit(1)

    variables, zones = parse_tec(FIELD_TEC)

    if TARGET_ZONE not in zones:
        available = ", ".join(zones.keys())
        print(f"ERROR: zone '{TARGET_ZONE}' not found. Available zones: {available}")
        sys.exit(1)

    zone_vars = extract_vars(variables, zones[TARGET_ZONE])

    for var in ("x", "y", "T"):
        if var not in zone_vars:
            print(f"ERROR: variable '{var}' not found in zone {TARGET_ZONE}.")
            print(f"       Available variables: {', '.join(variables)}")
            sys.exit(1)

    hdr      = zones[TARGET_ZONE]["header"]
    ni, nj   = hdr["ni"], hdr["nj"]
    ni_cells = ni - 1
    nj_cells = nj - 1

    x_cc, y_cc = build_cell_centres(zone_vars["x"], zone_vars["y"], ni, nj)
    T_2d = zone_vars["T"].reshape(nj_cells, ni_cells)

    # Bilinear interpolation at point A
    T_interp, cells = bilinear_interp(x_cc, y_cc, T_2d, TARGET_X, TARGET_Y)

    error_abs = T_interp - REFERENCE_T
    error_pct = abs(error_abs) / REFERENCE_T * 100.0

    sep = "-" * 48

    print()
    print(sep)
    print(" 2-D STEADY STATE THERMAL ANALYSIS ")
    print(" NAFEMS Benchmark T9 (i) ")
    print(sep)
    print(f"  {'Quantity':<20} {'Value (K)':>10}  {'Value (°C)':>10}")
    print(sep)
    print(f"  {'NAFEMS Reference':<20} "
          f"{REFERENCE_T:>10.3f}  {REFERENCE_T - 273.15:>10.3f}")
    print(f"  {'FUSS Solution':<20} "
          f"{T_interp:>10.3f}  {T_interp - 273.15:>10.3f}")
    print(sep)
    print(f"  {'Absolute error':<20} {error_abs:>10.3f}  {'K':>10}")
    print(f"  {'Relative error':<20} {error_pct:>10.4f}  {'%':>10}")
    print(sep)
    print()

    if error_pct < 1.0:
        print("  Result: PASS  (relative error < 1 %)")
    else:
        print("  Result: FAIL  (relative error >= 1 %)")

    print()
    save_colormap(x_cc, y_cc, T_2d, (TARGET_X, TARGET_Y), SVG_OUT)
    print()


if __name__ == "__main__":
    main()