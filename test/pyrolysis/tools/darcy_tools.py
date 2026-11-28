#!/usr/bin/env python3
"""Readers shared by the verify.py scripts of test/pyrolysis/darcy (plan 07): the Tecplot BLOCK field
file, OUTPUT/darcy.txt, OUTPUT/energy-budget.txt, the case's gas constants and permeability, and the
lagged Darcy coefficient Gamma = p W K/(R T mu) FUSS assembles its face fluxes with. Python 3.6, no numpy."""
import json
import math
import os
import re

R_UNIV = 8.314471469      # FUSS_Parameters_m


def read_tec(path):
    """Tecplot BLOCK file written by ORION: names from the header, cell/nodal split by name; returns
    ((I, J, K), SOLUTIONTIME, {name: [values]})."""
    with open(path) as fh:
        lines = fh.read().split("\n")
    names = [t.strip('"') for t in lines[0].split("=", 1)[1].split() if t.strip('"')]
    m = re.search(r"I=\s*(\d+),\s*J=\s*(\d+),\s*K=\s*(\d+)", lines[1])
    I, J, K = (int(x) for x in m.groups())
    ts = re.search(r"SOLUTIONTIME=([-+0-9.Ee]+)", lines[1])
    t = float(ts.group(1)) if ts else None
    vals = []
    for line in lines[2:]:
        s = line.strip()
        if s:
            try:
                vals.append(float(s))
            except ValueError:
                pass
    nn, nc = I * J * K, (I - 1) * (J - 1) * max(K - 1, 1)
    d, off = {}, 0
    for nm in names:
        n = nn if nm in ("x", "y", "z") else nc
        d[nm] = vals[off:off + n]
        off += n
    return (I, J, K), t, d


def cell_centres_x(dims, d):
    """x of every cell centre (i-major order as the cell arrays), from the node coordinates."""
    I, J, K = dims
    X = d["x"]
    nx, ny, nz = I - 1, J - 1, max(K - 1, 1)
    xc = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                corners = [X[(i + di) + (j + dj) * I + (k + dk) * I * J] for di in (0, 1) for dj in (0, 1) for dk in (0, 1)]
                xc.append(sum(corners) / 8.0)
    return xc


def cell_centres_y(dims, d):
    I, J, K = dims
    Y = d["y"]
    nx, ny, nz = I - 1, J - 1, max(K - 1, 1)
    yc = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                corners = [Y[(i + di) + (j + dj) * I + (k + dk) * I * J] for di in (0, 1) for dj in (0, 1) for dk in (0, 1)]
                yc.append(sum(corners) / 8.0)
    return yc


def read_darcy(path):
    """rows of OUTPUT/darcy.txt: [(iter, time, block, sweeps, resid, dp_last), ...]"""
    rows = []
    for line in open(path):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        p = s.split()
        rows.append((int(p[0]), float(p[1]), int(p[2]), int(p[3]), float(p[4]), float(p[5])))
    return rows


def read_budget(path):
    """rows of OUTPUT/energy-budget.txt: [(t, E_tot, M_s, Q_out, Q_gas_out, M_gas_out, Q_src), ...]"""
    rows = []
    for line in open(path):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        rows.append(tuple(float(x) for x in s.split()))
    return rows


def gas_constants(case):
    """(W, mu, h) of the first row of INPUT/gas.dat (the constant synthetic gas of plan 07's cases)."""
    zone = None
    for line in open(os.path.join(case, "INPUT", "gas.dat")):
        if line.startswith("ZONE"):
            zone = line
            continue
        p = line.split()
        if zone is not None and len(p) == 5:
            try:
                T, pr, M, h, mu = [float(x) for x in p]
            except ValueError:
                continue
            return M, mu, h
    raise SystemExit("gas.dat: no data row")


def case_meta(case):
    return json.load(open(os.path.join(case, "case.json")))


def ini_value(case, key):
    """value of `key = value` in input.ini (first match), or None"""
    for line in open(os.path.join(case, "input.ini")):
        s = line.split(";")[0].strip()
        if s.startswith(key) and "=" in s:
            k, v = s.split("=", 1)
            if k.strip() == key:
                return v.strip()
    return None


def gamma(p, T, W, mu, K):
    """Gamma = rho_g K/mu = p W K/(R T mu) [kg/(m Pa s)], scalar K"""
    return p * W * K / (R_UNIV * T * mu)


def exact_profile(x, L, pw, pi, T, W, mu, K):
    """plan 07 section 4.2: p^2 = p_w^2 + (pi R T mu/(W K)) (L^2 - x^2)"""
    a = pi * R_UNIV * T * mu / (W * K)
    return math.sqrt(pw * pw + a * (L * L - x * x))


def selftest():
    return os.environ.get("FUSS_GATE_SELFTEST", "0") == "1"
