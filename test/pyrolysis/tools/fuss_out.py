"""Readers for FUSS outputs of the literature cases (2026-10-13): wall*.tec (per-face Tw, qw), probe .txt files, the energy budget."""
import glob
import os
import re


def read_wall(path):
    """{zone title: {'t': time, 'Tw': [...], 'qw': [...]}} of one wall*.tec file."""
    zones = {}
    cur = None
    vals = []

    def close():
        if cur is None:
            return
        I, J, K, t, title = cur
        nn = I * J * K
        nc = max(I - 1, 1) * max(J - 1, 1) * max(K - 1, 1)
        v = vals[3 * nn:]
        zones[title] = {"t": t, "Tw": v[:nc], "qw": v[nc:2 * nc]}
    for line in open(path):
        if line.lstrip().startswith("ZONE"):
            close()
            title = re.search(r"T\s*=\s*([^,\s]+)", line).group(1)
            I = int(re.search(r"\bI=(\d+)", line).group(1)); J = int(re.search(r"\bJ=(\d+)", line).group(1))
            K = int(re.search(r"\bK=(\d+)", line).group(1))
            t = float(re.search(r"SOLUTIONTIME=\s*([-+0-9.eEdD]+)", line).group(1).replace("D", "E").replace("d", "e"))
            cur = (I, J, K, t, title)
            vals = []
        elif cur is not None and not line.lstrip().startswith("VARIABLES"):
            vals += [float(s.replace("D", "E")) for s in line.split()]
    close()
    return zones


def wall_series(case, zone="B1F1"):
    out = []
    for f in glob.glob(os.path.join(case, "OUTPUT", "wall*.tec")):
        z = read_wall(f)
        if zone in z:
            out.append((z[zone]["t"], sum(z[zone]["Tw"]) / len(z[zone]["Tw"]), sum(z[zone]["qw"]) / len(z[zone]["qw"])))
    return sorted(out)


def probe(case, name):
    rows = []
    for line in open(os.path.join(case, "OUTPUT", name + ".txt")):
        s = line.split()
        if s:
            rows.append([float(v) for v in s])
    return rows


def budget(case):
    rows = []
    for line in open(os.path.join(case, "OUTPUT", "energy-budget.txt")):
        if line.startswith("#"):
            continue
        s = line.split()
        if s:
            rows.append([float(v) for v in s])
    return rows      # t E_tot M_s Q_out Q_gas_out M_gas_out Q_src [M_g]
