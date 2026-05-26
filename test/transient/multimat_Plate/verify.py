import os
import re
import sys

import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

BASE_DIR = os.path.dirname(__file__)

TARGET_FILES = [
    ("field1.tec",   0.005),
    ("field21.tec",  0.105),
    ("field41.tec",  0.205),
    ("field61.tec",  0.305),
    ("field81.tec",  0.405),
    ("field101.tec", 0.505),
]

SVG_OUT = os.path.join(BASE_DIR, "2D-transient-multimat-plate-comparison.svg")

# Simulation BC values (used for normalisation)
T_HOT  = 101.0   # K – hot wall (y = 0,  face3)
T_COLD = 100.0   # K – cold wall (y = 2, face4)

# Analytical solution parameters (matches analytical.py)
T_0     = 1.0        # normalised hot-wall temperature
T_L     = 0.0        # normalised cold-wall temperature
L       = 2.0        # domain height [m]
ALPHA   = 1.0        # thermal diffusivity [m²/s]
N_TERMS = 100        # Fourier series terms

# Colour cycle (one per time snapshot)
COLORS = ["black", "tab:red", "tab:green", "tab:blue", "tab:orange", "tab:purple"]


# ---------------------------------------------------------------------------
# Analytical solution
# ---------------------------------------------------------------------------

def analytical_T(y_arr, t):
    T = np.zeros_like(y_arr, dtype=float)
    for j, yj in enumerate(y_arr):
        T_ss = T_0 + yj / L * (T_L - T_0)
        T_tr = 0.0
        for r in range(1, N_TERMS + 1):
            A_r = 2.0 / (r * np.pi) * ((-1) ** r * T_L - T_0)
            T_tr += A_r * np.sin(r * np.pi * yj / L) * np.exp(
                -r ** 2 * np.pi ** 2 * ALPHA * t / L ** 2
            )
        T[j] = T_ss + T_tr
    return T


# ---------------------------------------------------------------------------
# Tecplot parser – Block 2 T profile along J at central I
# ---------------------------------------------------------------------------

def parse_block2_profile(filepath, i_central=None):

    with open(filepath, "r") as fh:
        lines = fh.readlines()

    # Locate Block2 zone header
    b2_idx = None
    for idx, line in enumerate(lines):
        if re.match(r"\s*ZONE", line) and "Block2" in line:
            b2_idx = idx
            break

    if b2_idx is None:
        raise RuntimeError(f"Block2 zone not found in {filepath}")

    zone_line = lines[b2_idx]

    # Parse dimensions
    ni = int(re.search(r"\bI\s*=\s*(\d+)", zone_line).group(1))
    nj = int(re.search(r"\bJ\s*=\s*(\d+)", zone_line).group(1))
    nk = int(re.search(r"\bK\s*=\s*(\d+)", zone_line).group(1))

    sol_time = float(
        re.search(r"SOLUTIONTIME\s*=\s*([\d.E+\-]+)", zone_line, re.IGNORECASE).group(1)
    )

    ni_c = ni - 1   # cells in I
    nj_c = nj - 1   # cells in J
    nk_c = nk - 1   # cells in K (= 1 for 2-D plane)

    n_nodal = ni * nj * nk
    n_cell  = ni_c * nj_c * nk_c   # 400 for Block2

    if i_central is None:
        i_central = ni_c // 2   # default: middle column

    # Collect float values after the zone header line
    needed = 3 * n_nodal + n_cell   # x, y, z (nodal) + T (cell-centred)
    raw = []
    line_idx = b2_idx + 1
    while len(raw) < needed and line_idx < len(lines):
        try:
            raw.append(float(lines[line_idx].strip()))
        except ValueError:
            pass
        line_idx += 1

    if len(raw) < needed:
        raise RuntimeError(
            f"Not enough data in Block2 of {filepath}: "
            f"need {needed}, got {len(raw)}"
        )

    raw = np.array(raw)
    y_flat = raw[n_nodal : 2 * n_nodal]   # y is the second nodal variable
    T_flat = raw[3 * n_nodal : 3 * n_nodal + n_cell]

    # In BLOCK format the fastest index is I.
    # Use the k=0 layer of nodal data: first ni*nj values.
    y_nodes_k0 = y_flat[: ni * nj].reshape(nj, ni)

    # Cell-centre y at the chosen I column
    y_cc = 0.5 * (y_nodes_k0[:nj_c, i_central] + y_nodes_k0[1 : nj_c + 1, i_central])

    # T at the central I column: T_flat[j * ni_c + i_central]
    T_raw = np.array([T_flat[j * ni_c + i_central] for j in range(nj_c)])

    # Normalise to [0, 1] matching the analytical convention
    T_norm = (T_raw - T_COLD) / (T_HOT - T_COLD)

    return sol_time, y_cc, T_norm


# ---------------------------------------------------------------------------
# Comparison plot
# ---------------------------------------------------------------------------

def save_comparison_plot(results, svg_path):
    """Overlay analytical (solid) and FUSS (dashed + markers) T profiles.

    Parameters
    ----------
    results : list of (sol_time, y_cc, T_norm_fuss, y_an, T_an)
    """
    fig, ax = plt.subplots(figsize=(7.5, 5.5))

    # Fine y grid for smooth analytical curves
    y_fine = np.linspace(0.0, L, 300)

    for (sol_time, y_cc, T_fuss, _, _), color in zip(results, COLORS):
        T_an_fine = analytical_T(y_fine, sol_time)

        label = f"t = {sol_time:.3f} s"
        ax.plot(y_fine, T_an_fine, color=color, linewidth=1.6,
                label=f"{label} - Analytical")
        ax.plot(y_cc, T_fuss, color=color, linewidth=0.0,
                marker="o", markersize=4, markerfacecolor=color,
                markeredgewidth=1.2, label=f"{label} - FUSS")

    ax.set_xlabel("y  (m)", fontsize=11)
    ax.set_ylabel("Normalised temperature  T*", fontsize=10)
    ax.set_title(
        "Temperature Profile for Conductive Material",
        fontsize=12, fontweight="bold"
    )

    ax.set_xlim([0.0, L])
    ax.set_ylim(-0.005, 1.005)
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.05))
    ax.tick_params(axis="both", which="major", direction="in", length=6, labelsize=9)
    ax.tick_params(axis="both", which="minor", direction="in", length=3)

    # Two-column legend (analytical / FUSS pairs grouped by colour)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, fontsize=7.5, ncol=2,
              loc="upper right", framealpha=0.85)

    fig.tight_layout()
    fig.savefig(svg_path, format="svg", bbox_inches="tight")
    plt.close(fig)
    print(f"  Comparison plot saved → {svg_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    # Check all target files exist
    missing = []
    for fname, _ in TARGET_FILES:
        fp = os.path.join(BASE_DIR, "SOLUTION", fname)
        if not os.path.isfile(fp):
            missing.append(fp)
    if missing:
        for p in missing:
            print(f"ERROR: file not found – {p}")
        sys.exit(1)

    sep = "-" * 72

    print()
    print(sep)
    print(" 2-D TRANSIENT THERMAL ANALYSIS ")
    print(" Multi-Material Plate – Verification")
    print(sep)
    print(f"  {'File':<14} {'t_nominal':>10} {'t_file':>10} "
          f"{'L-inf err':>12} {'L2 err':>12}")
    print(sep)

    results       = []
    linf_errors   = []
    l2_errors     = []

    for fname, t_nominal in TARGET_FILES:
        fp = os.path.join(BASE_DIR, "SOLUTION", fname)
        sol_time, y_cc, T_fuss = parse_block2_profile(fp)

        # Evaluate analytical solution at the same y positions
        T_an = analytical_T(y_cc, sol_time)

        err      = np.abs(T_fuss - T_an)
        linf_err = float(err.max())
        l2_err   = float(np.sqrt(np.mean(err ** 2)))

        linf_errors.append(linf_err)
        l2_errors.append(l2_err)
        results.append((sol_time, y_cc, T_fuss, y_cc, T_an))

        print(f"  {fname:<14} {t_nominal:>10.3f} {sol_time:>10.4f} "
              f"{linf_err:>12.6f} {l2_err:>12.6f}")

    print(sep)
    max_linf = max(linf_errors)
    mean_l2  = float(np.mean(l2_errors))
    print(f"  {'Max L-inf error':<44} {max_linf:>12.6f}")
    print(f"  {'Mean L2 error':<44} {mean_l2:>12.6f}")
    print(sep)
    print()

    if max_linf < 0.02:
        print(f"  Result: PASS  (max L-inf error < 2 %)")
    else:
        print(f"  Result: FAIL  (max L-inf error >= 2 % )")

    print()

    save_comparison_plot(results, SVG_OUT)
    print()


if __name__ == "__main__":
    main()