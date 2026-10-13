"""lit_case.py -- write 1-D FUSS pyrolysis cases whose material tables come from analytic functions.

Shared by the literature-benchmark generators (make_baum_atreya.py, make_zhang2010.py; 2026-10-13). The file
formats are those of make_tacot_tables.py / make_case.py / make_darcy_case.py; only the material functions are new.
A case is a slab of nx cells along i (2 along j, 1 along k, as every 1-D pyrolysis case).
Python 3.6, no third-party modules.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_case as mc            # noqa: E402
import make_darcy_case as mdc     # noqa: E402
import make_tacot_tables as mtt   # noqa: E402
import make_exact_case as mec     # noqa: E402

SIGMA = 5.67e-8          # FUSS_Parameters_m sigma_SB
NY, DY, DZ = 2, 0.1, 0.1


def write_solid(path, name, Tlo, Thi, cp, h, k, emis=0.9, absorp=0.9):
    rows = [(float(T), [cp(T), h(T), k(T), k(T), k(T), emis, absorp]) for T in range(Tlo, Thi + 1)]
    mtt.write_table(path, mtt.SOLID_VARS, [(name, rows)])


def write_gas(path, name, Tlo, Thi, M, h, mu, plevels):
    zones = []
    for p in plevels:
        rows = [(float(T), [p, M, h(T), mu(T)]) for T in range(Tlo, Thi + 1)]
        zones.append(("%s p=%s" % (name, mtt.FMT % p), rows))
    mtt.write_table(path, mtt.GAS_VARS, zones)


def write_pyro(path, name, phases, eps_g_v, eps_g_c, perm=None, provenance="test/pyrolysis/tools/lit_case.py"):
    const = {"nSolidPhases": len(phases), "phases": phases, "eps_g_v": eps_g_v, "eps_g_c": eps_g_c}
    mtt.write_pyrolysis_dat(path, name, const, provenance, perm=perm)


def reaction(F, A, E, m=1, b=0, Tthr=0.0):
    return {"F": "%.15g" % F, "A": "%.15g" % A, "E": "%.15g" % E, "m": "%.15g" % m, "n": "%.15g" % b,
            "T": "%.15g" % Tthr, "h": "0"}


def phase(rhoI, epsI, reactions):
    return {"rhoI": "%.15g" % rhoI, "epsI": "%.15g" % epsI, "nPyroReac": "%d" % len(reactions), "reactions": reactions}


def tr_305(alpha_q, h, eps, Tinf, Tsur=None):
    """Tref of ONE 305 line, q = h (Tr - T) + eps sigma (Tr^4 - T^4), equal for every T to
    alpha_q + h (Tinf - T) + eps sigma (Tsur^4 - T^4) (Tsur None = Tinf; 0 = no back-radiation term)."""
    Ts = Tinf if Tsur is None else Tsur
    rhs = alpha_q + h * Tinf + eps * SIGMA * Ts ** 4
    lo, hi = 0.0, 1.0e4
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if h * mid + eps * SIGMA * mid ** 4 < rhs:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def write_case(case, name, xs, face1, face2, Tcells, sections, probes=None, meta=None, probe_vars="T",
               masterdir=None, time_files=None):
    """xs: nodes along i; face1/face2: (code, value line) of bc.txt for i = 1 / i = nx; Tcells: initial T per cell;
    sections: ini lines after the header; probes: [(name, i, dtime)] (cell index 1-based, j = 1); masterdir: absolute
    MASTERDIR for FUSS.sh (a case outside the tree), else the relative depth of the case in the tree;
    time_files: {name: [(t, value)]} written beside input.ini (BC time series)."""
    case = os.path.abspath(case)
    nx = len(xs) - 1
    ys = [DY * j for j in range(NY + 1)]
    zs = [0.0, DZ]
    for d in ("INPUT", "MESH"):
        os.makedirs(os.path.join(case, d), exist_ok=True)
    mec.write_ic_T(os.path.join(case, "INPUT", "ic.tec"), xs, ys, zs, Tcells)
    mdc.write_ic(os.path.join(case, "MESH", "mesh.tec"), xs, ys, zs, "Block1", False, 0.0)
    spec = {1: lambda i, j, k: face1, 2: lambda i, j, k: face2,
            3: lambda i, j, k: (300, None), 4: lambda i, j, k: (300, None),
            5: lambda i, j, k: (0, None), 6: lambda i, j, k: (0, None)}
    mdc.write_bc(os.path.join(case, "INPUT", "bc.txt"), nx, NY, spec)
    with open(os.path.join(case, "INPUT", "phase.txt"), "w") as fh:
        fh.write("solid phase\n%s 1\n" % name)
    virgin = mc.read_table_zone1(os.path.join(case, "INPUT", "virgin.dat"), 7)
    _, rho_v, _ = mc.read_pyrolysis_dat(os.path.join(case, "INPUT", "pyrolysis.dat"))
    mc.write_properties_dat(os.path.join(case, "INPUT", "properties.dat"), name, virgin, rho_v)
    lines = list(sections)
    if probes:
        lines += ["", "[FUSS-Probes]"] + ["probe%d = %s" % (n + 1, p[0]) for n, p in enumerate(probes)]
        for pname, i, dtime in probes:
            lines += ["", "[%s]" % pname, "index-position = 1 %d 1 1" % i, "variables = %s" % probe_vars,
                      "dtime = %g" % dtime]
    with open(os.path.join(case, "input.ini"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    for fname, pairs in (time_files or {}).items():
        mc.write_time_file(os.path.join(case, fname), pairs)
    if masterdir:
        with open(mc.DONOR_FUSS_SH) as fh:
            src = fh.read().split("\n")
        with open(os.path.join(case, "FUSS.sh"), "w") as fh:
            fh.write("\n".join("MASTERDIR=%s" % masterdir if s.startswith("MASTERDIR=") else s for s in src))
        os.chmod(os.path.join(case, "FUSS.sh"), 0o755)
    else:
        depth, _ = mc.case_depth(case, "adiabatic0d", 1000.0, None)
        mc.write_fuss_sh(os.path.join(case, "FUSS.sh"), mc.DONOR_FUSS_SH, depth)
    if meta is not None:
        with open(os.path.join(case, "case.json"), "w") as fh:
            json.dump(meta, fh, indent=1, sort_keys=True)
            fh.write("\n")


def ini_common(title_lines, tend, vnn, sol_dt, pyro_lines, darcy_lines=None, generator="test/pyrolysis/tools/lit_case.py"):
    out = ["; " + "=" * 77] + [";  " + s for s in title_lines] + [";", ";  Generated by %s -- edit that, not this." % generator,
                                                                    "; " + "=" * 77]
    out += ["", "[FUSS-Parameters]", "newrun = true", "res-threshold = 1e-30", "time-threshold = %.15g" % tend,
            "iter-threshold = 1000000000"]
    out += ["", "[FUSS-Numerics]", "time-scheme = RK3", "vnn = %g" % vnn, "time-accurate = true",
            "integration-variables = cons", "irs = false", ""] + pyro_lines
    if darcy_lines:
        out += [""] + darcy_lines
    out += ["", "[FUSS-IO]", "sol-dtime     = %.15g" % sol_dt, "sol-overwrite = false", "sol-format    = tecplot ascii",
            "ini-diter     = 1000000000"]
    return out
