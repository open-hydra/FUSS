"""
Verification of the variable_cp_Plate test: uniform heating of an adiabatic plate with a temperature-dependent
specific heat, cp(T) = A + B T (make_properties.py), constant density, uniform source qvol.

With every wall adiabatic and a uniform source the temperature stays uniform, and the heat equation reduces to

    rho cp(T) dT/dt = qvol   =>   rho [ A (T - T0) + B/2 (T^2 - T0^2) ] = qvol t

whose positive root is the reference T(t). The heat capacity is what is tested: an update that used
d(rho cp T)/dT = rho (A + 2 B T) instead of rho cp would follow rho [ A (T - T0) + B (T^2 - T0^2) ] = qvol t,
much slower (shown for reference). The volumetric energy stored, rho [A T + B T^2/2] (the Energy column), must
grow by qvol t exactly.

Reads every OUTPUT/field*.tec (sol-overwrite = false: field<iteration>.tec every sol-diter, field.tec at the
end), prints the comparison and writes the SVG plot used in the documentation.
"""
import glob
import os
import re
import sys

import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SVG_OUT = os.path.join(HERE, "2D-transient-variable-cp-plate.svg")

A, B = 300.0, 0.4        # cp = A + B T [J/(kg K)]   (as in make_properties.py)
RHO = 7850.0             # [kg/m3]
QVOL = 5.0e7             # [W/m3]   (input.ini, STB-Block1)
T0 = 300.0               # [K]      (input.ini, ICB-Block1)
TOL = 1e-3               # largest error allowed, relative to the temperature rise


def T_exact(t):
    """Uniform temperature with heat capacity rho cp(T): rho (A T + B T^2/2) grows by qvol t."""
    c = RHO * (A * T0 + 0.5 * B * T0**2) + QVOL * np.asarray(t)
    return (-A + np.sqrt(A**2 + 2.0 * B * c / RHO)) / B


def T_dcpT(t):
    """The same with heat capacity d(rho cp T)/dT: rho (A T + B T^2) grows by qvol t (reference only)."""
    c = RHO * (A * T0 + B * T0**2) + QVOL * np.asarray(t)
    return (-A + np.sqrt(A**2 + 4.0 * B * c / RHO)) / (2.0 * B)


def energy(T):
    """Volumetric enthalpy of make_properties.py: h(T) = rho (A T + B T^2 / 2)."""
    return RHO * (A * np.asarray(T) + 0.5 * B * np.asarray(T)**2)


def read_field(path):
    """Solution time and cell temperatures of an ASCII Tecplot field written by FUSS (one or more zones)."""
    with open(path) as fh:
        text = fh.read()
    names = [a or b for a, b in re.findall(r'"([^"]*)"|([^\s,"]+)',
                                            re.search(r'(?i)variables\s*=\s*(.*)', text).group(1))]
    t, temps = None, []
    for zone in re.split(r'(?im)^\s*zone', text)[1:]:
        header, _, body = zone.partition('\n')
        while not re.match(r'\s*[-+0-9.]', body.split('\n', 1)[0]):   # continuation lines of the header
            line, _, body = body.partition('\n')
            header += ' ' + line
        ni, nj, nk = (int(re.search(rf'\b{d}\s*=\s*(\d+)', header).group(1)) for d in 'IJK')
        m = re.search(r'SOLUTIONTIME\s*=\s*([-+0-9.EeDd]+)', header)
        t = float(m.group(1).replace('D', 'E').replace('d', 'e')) if m else t
        cc = re.search(r'\[(\d+)-(\d+)\]\s*=\s*CELLCENTERED', header)
        first_cc = int(cc.group(1)) if cc else len(names) + 1
        values = np.array(body.split(), dtype=float)
        n_node, n_cell = ni * nj * nk, max(ni - 1, 1) * max(nj - 1, 1) * max(nk - 1, 1)
        offset = 0
        for v, name in enumerate(names, start=1):
            n = n_cell if v >= first_cc else n_node
            if name == 'T':
                temps.append(values[offset:offset + n])
            offset += n
    return t, np.concatenate(temps)


def main():
    files = sorted(glob.glob(os.path.join(HERE, "OUTPUT", "field*.tec")))
    if not files:
        print("ERROR: no OUTPUT/field*.tec: run ./FUSS.sh solve first")
        sys.exit(1)
    data = sorted(read_field(f) for f in files)
    t = np.array([d[0] for d in data])
    T_mean = np.array([d[1].mean() for d in data])
    spread = max(np.ptp(d[1]) for d in data)
    T_ref = T_exact(t)
    rel = np.abs(T_mean - T_ref) / np.maximum(T_ref - T0, 1e-30)
    t_end, T_end = t[-1], T_mean[-1]
    e_err = (energy(T_end) - energy(T0) - QVOL * t_end) / (QVOL * t_end)

    sep = "-" * 62
    print()
    print(sep)
    print(" 2-D TRANSIENT: UNIFORM HEATING, cp(T) = A + B T, ADIABATIC WALLS")
    print(sep)
    print(f"  {'solutions read':<38} {len(files):>10d}")
    print(f"  {'final time [s]':<38} {t_end:>10.3f}")
    print(f"  {'FUSS temperature [K]':<38} {T_end:>10.3f}")
    print(f"  {'exact, heat capacity rho cp [K]':<38} {T_exact(t_end):>10.3f}")
    print(f"  {'with d(rho cp T)/dT instead [K]':<38} {T_dcpT(t_end):>10.3f}")
    print(f"  {'largest spread across the plate [K]':<38} {spread:>10.2e}")
    print(f"  {'largest error / temperature rise':<38} {rel.max():>10.2e}")
    print(f"  {'stored energy error / (qvol t)':<38} {e_err:>10.2e}")
    print(sep)
    ok = rel.max() < TOL and spread < 1e-6 * T_end
    print(f"  Result: {'PASS' if ok else 'FAIL'}  (error < {TOL:g} of the rise, uniform field)")
    print()

    tt = np.linspace(0.0, t_end, 200)
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.plot(tt, T_exact(tt), color="#2a78d6", lw=2, label=r"exact, $\rho c_p\,\partial T/\partial t = \dot q$")
    ax.plot(tt, T_dcpT(tt), color="#eb6834", lw=2, ls="--",
            label=r"heat capacity $d(\rho c_p T)/dT$ (wrong)")
    ax.plot(t, T_mean, "o", color="#1baf7a", ms=5, label="FUSS")
    ax.set_xlabel("time [s]")
    ax.set_ylabel("temperature [K]")
    ax.grid(True, color="#e4e3dc", lw=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(SVG_OUT)
    print(f"  plot: {SVG_OUT}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
