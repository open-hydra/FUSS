import os
import sys

import matplotlib
matplotlib.use("svg")           # non-interactive backend for SVG export
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

BASE_DIR  = os.path.dirname(__file__)
P1_FILE   = os.path.join(BASE_DIR, "OUTPUT", "p1.txt")
P2_FILE   = os.path.join(BASE_DIR, "OUTPUT", "p2.txt")
TW_FILE   = os.path.join(BASE_DIR, "Tw_time.dat")

SVG_BC    = os.path.join(BASE_DIR, "1D-transient-oscillating-rod-bc.svg")
SVG_PROBE = os.path.join(BASE_DIR, "1D-transient-oscillating-rod-probes.svg")

T_BASELINE   = 2000.0  # K
REFERENCE_DT = 94.69

def load_probe(filepath):
    data = np.loadtxt(filepath)
    return data[:, 0], data[:, 1]


def load_tw(filepath):
    times, temps = [], []
    with open(filepath, "r") as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped or stripped.lower() == "periodic":
                continue
            parts = stripped.split()
            times.append(float(parts[0]))
            temps.append(float(parts[1]))
    return np.array(times), np.array(temps)


def save_bc_plot(t_tw, T_tw, svg_path):
    """Save the oscillating boundary condition plot as an SVG."""
    fig, ax = plt.subplots(figsize=(8, 4))

    ax.plot(t_tw, T_tw, color="tab:red", linewidth=1.8,
            label=r"Right wall BC  $T_w(t) = a\,(\sin\omega t + \theta)$")
    ax.axhline(T_BASELINE, color="tab:blue", linewidth=1.4, linestyle="--",
               label="Left wall BC (fixed, $T = 2000$ K)")

    ax.set_xlabel("Time  (s)", fontsize=11)
    ax.set_ylabel("Temperature  (K)", fontsize=11)
    ax.set_title(
        "1-D Oscillating Transient Thermal Analysis – Boundary Conditions\n"
        r"$\omega/2\pi = 0.0125$ Hz,  $\theta = 0$,  $a = 1000$ K",
        fontsize=10, fontweight="bold"
    )
    ax.legend(fontsize=9)
    ax.tick_params(labelsize=9)
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.0f"))
    ax.set_xlim([0.0, 80.0])

    fig.tight_layout()
    fig.savefig(svg_path, format="svg", bbox_inches="tight")
    plt.close(fig)
    print(f"  BC plot saved → {svg_path}")


def save_probe_plot(t, T_mean, T_p1, T_p2, T_ref, svg_path):
    """Save the probe temperature evolution plot as an SVG."""
    fig, ax = plt.subplots(figsize=(8, 4))

    ax.plot(t, T_p1, color="tab:orange", linewidth=1.0,
            alpha=0.55, label="probe 1")
    ax.plot(t, T_p2, color="tab:green", linewidth=1.0,
            alpha=0.55, label="probe 2")
    ax.plot(t, T_mean, color="tab:blue", linewidth=2.0,
            label="mean probe")
    ax.axhline(T_ref + T_BASELINE, color="tab:red", linewidth=1.4, linestyle="--",
               label=f"NAFEMS Reference")

    ax.set_xlabel("Time  (s)", fontsize=11)
    ax.set_ylabel("Temperature  (K)", fontsize=11)
    ax.set_title(
        "1-D Oscillating Transient Thermal Analysis – Probe Temperature at x ≈ 0.09 m",
        fontsize=10, fontweight="bold"
    )
    ax.legend(fontsize=9)
    ax.tick_params(labelsize=9)
    ax.set_xlim([0.0, 160.0])
    ax.set_ylim([1900.0, 2150.0])
    
    fig.tight_layout()
    fig.savefig(svg_path, format="svg", bbox_inches="tight")
    plt.close(fig)
    print(f"  Probe plot saved → {svg_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():

    for fp in (P1_FILE, P2_FILE, TW_FILE):
        if not os.path.isfile(fp):
            print(f"ERROR: file not found – {fp}")
            sys.exit(1)

    print("  Loading probe data …", end="", flush=True)
    t1, T1 = load_probe(P1_FILE)
    t2, T2 = load_probe(P2_FILE)
    print(" done.")

    # Build common time axis and mean temperature
    if len(t1) == len(t2) and np.allclose(t1, t2):
        t_common = t1
        T_mean   = 0.5 * (T1 + T2)
    else:
        t_common = t1
        T2_interp = np.interp(t_common, t2, T2)
        T_mean = 0.5 * (T1 + T2_interp)

    T_max  = float(np.max(T_mean))
    dT_max = T_max - T_BASELINE           # peak temperature rise above baseline [K = °C]

    error_abs = dT_max - REFERENCE_DT
    error_pct = abs(error_abs) / REFERENCE_DT * 100.0

    sep = "-" * 38

    print()
    print(sep)
    print(" 1-D TRANSIENT THERMAL ANALYSIS ")
    print(" NAFEMS Benchmark T5 ")
    print(sep)
    print(f"  {'Quantity':<20} {'Value':>10}")
    print(sep)
    print(f"  {'NAFEMS Reference':<20} "
          f"{REFERENCE_DT:>10.3f}")
    print(f"  {'FUSS Solution':<20} "
          f"{dT_max:>10.3f}")
    print(sep)
    print(f"  {'Absolute error':<20} {error_abs:>10.3f}")
    print(f"  {'Relative error':<20} {error_pct:>10.4f}  {'%':>1}")
    print(sep)
    print()

    if error_pct < 1.0:
        print("  Result: PASS  (relative error < 1 %)")
    else:
        print("  Result: FAIL  (relative error >= 1 %)")

    print()
    t_tw, T_tw = load_tw(TW_FILE)
    save_bc_plot(t_tw, T_tw, SVG_BC)
    save_probe_plot(t_common, T_mean, T1, T2, REFERENCE_DT, SVG_PROBE)
    print()


if __name__ == "__main__":
    main()
