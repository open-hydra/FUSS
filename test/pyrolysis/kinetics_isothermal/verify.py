#!/usr/bin/env python3
"""Gate 4.2 of plan 06: isothermal kinetics against the closed-form propagator.

A 4-cell slab held at T0 by BC 302 on both i-faces, started at T0, with the
`--isothermal` test tables (cp, h, k constant; h and the gas h_pg the same
linear function of T; eps_g_v = eps_g_c). Then T stays at T0 exactly and every
reaction follows the closed form of d chi/dt = kappa (1 - chi)^m at fixed T:

    1 - chi(t) = [1 + (m - 1) kappa t]^(-1/(m-1))      (m /= 1, chi(0) = 0)
    1 - chi(t) = exp(-kappa t)                         (m = 1)

with kappa = A T0^n exp(-E/(R T0)) above the reaction's threshold, 0 below.

Checks on every OUTPUT/field<i>.tec (t from SOLUTIONTIME), every cell:
  (a) max |chi_r - chi_exact(t)|  <= 1e-8        (exact propagator; the only
      error is the Newton tolerance on T propagated into kappa)
  (b) max |T - T0|                <= 1e-7 K
  (c) rho_s = rho_v - sum_r eps_I rho_I F_r chi_r   to 1e-10 relative
  (d) the type-1 line integral: with uniform cells and a uniform state,
      mdot_pg(s) = s * mdot_pg(1) along every line, to 1e-12 relative, and
      mdot_pg(1) > 0 while any chi_r < 1.

Usage: verify.py [case_dir]   (default: the directory this file sits in / T0700)
Exit 0 on PASS, 1 on FAIL. Python 3.6, no numpy.
"""
import math, os, re, sys

R_UNIV = 8.314471469
TOL_CHI, TOL_T, TOL_RHO, TOL_MDOT = 1.0e-8, 1.0e-7, 1.0e-10, 1.0e-12


def read_tec(path):
    """Tecplot BLOCK file written by ORION: names from the header, cell/nodal split by name."""
    with open(path) as fh:
        lines = fh.read().split("\n")
    names = [t.strip('"') for t in lines[0].split("=", 1)[1].split() if t.strip('"')]
    m = re.search(r"I=(\d+),\s*J=(\d+),\s*K=(\d+)", lines[1])
    I, J, K = (int(x) for x in m.groups())
    ts = re.search(r"SOLUTIONTIME=([-+0-9.Ee]+)", lines[1])
    t = float(ts.group(1)) if ts else None
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
    return (I - 1, J - 1, max(K - 1, 1)), t, data


def read_pyrolysis_dat(path):
    """pyrolysis.dat (plan 06 section 3.2 grammar): one material, its phases and reactions."""
    kv = {}
    for line in open(path):
        line = line.split("//", 1)[0].strip()
        if not line:
            continue
        key, val = line.split(None, 1)
        kv[key] = val.strip()
    nph = int(kv["nSolidPhases"])
    rho_I = [float(kv["rhoI[%d]" % i]) for i in range(1, nph + 1)]
    eps_I = [float(kv["epsI[%d]" % i]) for i in range(1, nph + 1)]
    reac = []
    for i in range(1, nph + 1):
        for j in range(1, int(kv["nPyroReac[%d]" % i]) + 1):
            g = lambda k: float(kv["%s[%d][%d]" % (k, i, j)])
            reac.append(dict(F=g("F"), A=g("A"), E=g("E"), m=g("m"), n=g("n"), T=g("T"),
                             eps_rho_F=eps_I[i - 1] * rho_I[i - 1] * g("F")))
    rho_v = sum(e * r for e, r in zip(eps_I, rho_I))
    return rho_v, reac


def chi_exact(rc, T0, t):
    if T0 < rc["T"] or rc["A"] <= 0.0:
        return 0.0
    kappa = rc["A"] * (T0 ** rc["n"]) * math.exp(-rc["E"] / (R_UNIV * T0))
    if rc["m"] == 1.0:
        return 1.0 - math.exp(-kappa * t)
    br = 1.0 + (rc["m"] - 1.0) * kappa * t
    if br <= 0.0:
        return 1.0
    return 1.0 - br ** (-1.0 / (rc["m"] - 1.0))


def main():
    case = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "T0700")
    rho_v, reac = read_pyrolysis_dat(os.path.join(case, "INPUT", "pyrolysis.dat"))
    dims0, _, ic = read_tec(os.path.join(case, "INPUT", "ic.tec"))
    T0 = ic["T"][0]
    if max(ic["T"]) - min(ic["T"]) != 0.0:
        print("FAIL: the initial condition is not uniform"); return 1
    files = sorted((f for f in os.listdir(os.path.join(case, "OUTPUT")) if re.match(r"field\d*\.tec$", f)),
                   key=lambda f: int(re.sub(r"\D", "", f) or 0))
    if not files:
        print("FAIL: no OUTPUT/field*.tec"); return 1

    worst_chi = worst_T = worst_rho = worst_mdot = 0.0
    ok = True
    nfiles = 0
    for f in files:
        (ni, nj, nk), t, d = read_tec(os.path.join(case, "OUTPUT", f))
        if t is None:
            print("FAIL: %s has no SOLUTIONTIME" % f); return 1
        nfiles += 1
        ncell = ni * nj * nk
        chi = [d["chi_%d" % (r + 1)] for r in range(len(reac))]
        for c in range(ncell):
            for r, rc in enumerate(reac):
                worst_chi = max(worst_chi, abs(chi[r][c] - chi_exact(rc, T0, t)))
            worst_T = max(worst_T, abs(d["T"][c] - T0))
            rho_s = rho_v - sum(rc["eps_rho_F"] * chi[r][c] for r, rc in enumerate(reac))
            worst_rho = max(worst_rho, abs(d["rho_s"][c] - rho_s) / rho_v)
        # (d) along every (j,k) line, s = i = 1..ni
        any_active = any(chi[r][c] < 1.0 for r in range(len(reac)) for c in range(ncell))
        for k in range(nk):
            for j in range(nj):
                base = j * ni + k * ni * nj
                m1 = d["mdot_pg"][base]
                if any_active and not m1 > 0.0:
                    print("FAIL: %s line j=%d k=%d: mdot_pg(1) = %g is not > 0 while reactions are active" % (f, j + 1, k + 1, m1))
                    ok = False
                for s in range(1, ni + 1):
                    ms = d["mdot_pg"][base + s - 1]
                    ref = s * m1
                    worst_mdot = max(worst_mdot, abs(ms - ref) / abs(ref) if ref != 0.0 else abs(ms))

    print("isothermal kinetics gate, T0 = %g K, %d reactions, %d files (t up to %g s)" % (T0, len(reac), nfiles, t))
    print("  max |chi - chi_exact|      = %.3e   (tol %.0e)" % (worst_chi, TOL_CHI))
    print("  max |T - T0|               = %.3e K (tol %.0e)" % (worst_T, TOL_T))
    print("  max |rho_s - rho_s(chi)|/rho_v = %.3e (tol %.0e)" % (worst_rho, TOL_RHO))
    print("  max |mdot_pg(s) - s mdot_pg(1)| rel = %.3e (tol %.0e)" % (worst_mdot, TOL_MDOT))
    ok = ok and worst_chi <= TOL_CHI and worst_T <= TOL_T and worst_rho <= TOL_RHO and worst_mdot <= TOL_MDOT
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
