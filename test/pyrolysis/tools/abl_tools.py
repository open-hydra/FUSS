"""Readers and checks for the ablating-wall (BC 306) cases of plan 08 (2026-10-13).

read_seb       OUTPUT/<prefix>seb.txt: one dict per row, keyed by the header names
faces_306      the 306 records of INPUT/bc.txt: (b, i, j, k, f) -> lambda, T_inf, B' file, env file
read_env       a 306 env table; env_at(env, t) evaluates it as FUSS does (v_i + w (v_{i+1} - v_i), clamped)
seb_closure    plan 08 section 3.2: (a) |F| / max|term| <= 10 seb-tol; (b) beta at the Picard's closed-form
               fixed point zeta/(e^zeta - 1) where none of PATO's beta = 1 exceptions applies; (c) B'_pg and
               mdot_ch consistent with mdot_pg, rhoeUeCH0 and beta
pato_rows      a PATO output/porousMat file (scalar/Ta_plot, scalar/Ta_surfacePatch, mass) or an Amaryllis file
interp         linear interpolation in t of a [(t, v)] series (None outside it)
fronts         PATO's virgin-98 % / char-2 % front depths (massIOModel.C:131-166) from a FUSS field file
"""
import math
import os
import re

SEB_TERMS = ("cond", "conv", "adv_pg", "adv_ch", "alpha_qrad", "eps_Tinf4", "m_eps_Tw4")
CHEM_ON = 1.0 - 1.0e-12          # FUSS's on-switch (Lib_BC_Fluxes_Wall_Ablating.f90 CHEM_ON)
INT_COLS = ("b", "i", "j", "k", "f", "iter", "nclamp_p", "nclamp_Bg", "nclamp_T", "nclamp_rc")


def read_seb(case, prefix=""):
    path = os.path.join(case, "OUTPUT", prefix + "seb.txt")
    rows = []
    with open(path) as fh:
        head = fh.readline()
        if not head.startswith("#"):
            raise SystemExit("%s: no header line" % path)
        names = head[1:].split()
        if names[-1] == "(SI)":
            names = names[:-1]
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            v = line.split()
            if len(v) != len(names):
                raise SystemExit("%s: a row with %d values, the header names %d" % (path, len(v), len(names)))
            rows.append({n: (int(x) if n in INT_COLS else float(x)) for n, x in zip(names, v)})
    return rows


def faces_306(case):
    """{(b, i, j, k, f): (lambda, T_inf, bprime_file, env_file)} from INPUT/bc.txt."""
    out = {}
    lines = open(os.path.join(case, "INPUT", "bc.txt")).read().splitlines()
    n = 0
    while n < len(lines):
        t = lines[n].split()
        if len(t) == 6 and t[5] == "306":
            f2 = [s.strip() for s in lines[n + 1].split(",")]
            out[tuple(int(x) for x in t[:5])] = (float(f2[0]), float(f2[1]), f2[2], f2[3])
            n += 2
            continue
        n += 1
    return out


def read_env(path):
    rows = []
    for line in open(path):
        s = line.strip()
        if not s or s.startswith("//"):
            continue
        rows.append([float(x) for x in s.split()])
    return rows


def env_at(env, t):
    """The env vector (p_w, rhoeUeCH0, h_r, q_rad, chem, hconv, T_edge) at t, as FUSS's time_table_t%eval."""
    if t <= env[0][0]:
        return env[0][1:]
    if t >= env[-1][0]:
        return env[-1][1:]
    i = 0
    while env[i + 1][0] < t:
        i += 1
    w = (t - env[i][0]) / (env[i + 1][0] - env[i][0])
    return [a + w * (b - a) for a, b in zip(env[i][1:], env[i + 1][1:])]


def seb_closure(case, seb_tol=1.0e-8, prefix=""):
    """Plan 08 section 3.2 on every row of seb.txt. Returns (ok, lines, worst) with worst = the three maxima."""
    rows = read_seb(case, prefix)
    faces = faces_306(case)
    envs = {}
    worst = {"F_rel": 0.0, "beta": 0.0, "Bpg_rel": 0.0, "mch_rel": 0.0}
    n_on = n_beta = 0
    for r in rows:
        key = (r["b"], r["i"], r["j"], r["k"], r["f"])
        lam, _tinf, _bf, ef = faces[key]
        if ef not in envs:
            envs[ef] = read_env(os.path.join(case, "INPUT", ef))
        chem = env_at(envs[ef], r["time"])[4]
        terms = [abs(r[k]) for k in SEB_TERMS]
        if max(terms) > 0:
            worst["F_rel"] = max(worst["F_rel"], abs(r["F"]) / max(terms))
        if chem < CHEM_ON:
            continue
        n_on += 1
        rho, beta, mpg, mch = r["rhoeUeCH0"], r["beta"], r["mdot_pg"], r["mdot_ch"]
        rho_f = max(rho, 1.0e-6)
        bpg = mpg / (rho_f * beta)
        worst["Bpg_rel"] = max(worst["Bpg_rel"], abs(r["Bpg"] - bpg) / max(abs(bpg), 1e-300))
        mch_x = r["Bch"] * rho_f * beta
        if mch != 0.0 or mch_x != 0.0:
            worst["mch_rel"] = max(worst["mch_rel"], abs(mch - mch_x) / max(abs(mch_x), 1e-300))
        bprime = bpg + r["Bch"]
        if lam == 0.0 or bprime == 0.0 or mpg <= 1.0e-6 or rho_f <= 1.0e-4:
            continue
        n_beta += 1
        zeta = 2.0 * lam * (mpg + mch) / rho_f
        worst["beta"] = max(worst["beta"], abs(beta - zeta / math.expm1(zeta)))
    lim = 10.0 * seb_tol
    ok_a = worst["F_rel"] <= lim
    ok_b = worst["beta"] <= lim
    ok_c = worst["Bpg_rel"] <= 1.0e-14 and worst["mch_rel"] <= 1.0e-14
    lines = [
        "(seb-a) max |F| / max|term|           %.3e  (<= %.0e, %d rows)  %s" % (worst["F_rel"], lim, len(rows), "PASS" if ok_a else "FAIL"),
        "(seb-b) max |beta - zeta/(e^zeta-1)|   %.3e  (<= %.0e, %d chem-on rows with blowing)  %s" % (worst["beta"], lim, n_beta, "PASS" if ok_b else "FAIL"),
        "(seb-c) B'pg, mdot_ch identities       %.3e, %.3e  (<= 1e-14 relative, %d chem-on rows)  %s"
        % (worst["Bpg_rel"], worst["mch_rel"], n_on, "PASS" if ok_c else "FAIL"),
    ]
    return ok_a and ok_b and ok_c, lines, worst


def pato_rows(path):
    out = []
    for line in open(path):
        s = line.strip()
        if not s or s.startswith(("//", "#")) or s[0].isalpha():
            continue
        try:
            out.append([float(x) for x in s.split()])
        except ValueError:
            pass
    return out


def interp(series, t):
    if not series or t < series[0][0] - 1e-12 or t > series[-1][0] + 1e-12:
        return None
    lo, hi = 0, len(series) - 1
    while hi - lo > 1:
        m = (lo + hi) // 2
        if series[m][0] <= t:
            lo = m
        else:
            hi = m
    (t0, v0), (t1, v1) = series[lo], series[hi]
    return v0 if t1 == t0 else v0 + (v1 - v0) * (t - t0) / (t1 - t0)


def bracket(series, t):
    """(t_lo, t_hi) of the two samples interp() uses at t, or None."""
    if not series or t < series[0][0] or t > series[-1][0]:
        return None
    lo, hi = 0, len(series) - 1
    while hi - lo > 1:
        m = (lo + hi) // 2
        if series[m][0] <= t:
            lo = m
        else:
            hi = m
    return series[lo][0], series[hi][0]


def fronts(path, L, rho_v, rho_c, var="rho_s"):
    """(time, virgin-98 % depth, char-2 % depth) of one FUSS tecplot-ascii field file, block 1, j = k = 1,
    by PATO's rule (massIOModel.C:131-166): cells from the back face up, linear between cell centres, depth
    from the initial surface L; (0, 0) while the surface cell is still above the 98 % density; the char
    front is the surface cell's centre until that cell is below the 2 % density."""
    toks = open(path).read().split()
    k = toks.index("ZONE")
    names = " ".join(toks[:k]).split("=", 1)[1].split()
    s = k
    while not toks[s].startswith("SOLUTIONTIME="):
        s += 1
    zone = " ".join(toks[k:s + 1])
    I = int(re.search(r"\bI=(\d+)", zone).group(1))
    J = int(re.search(r"\bJ=(\d+)", zone).group(1))
    K = int(re.search(r"\bK=(\d+)", zone).group(1))
    t = float(re.search(r"SOLUTIONTIME=([-+0-9.eEdD]+)", zone).group(1).replace("D", "E"))
    s += 1
    nn = I * J * K
    nc = (I - 1) * max(J - 1, 1) * max(K - 1, 1)
    x = [float(v) for v in toks[s:s + I]]
    iv = names.index(var)
    off = s + 3 * nn + (iv - 3) * nc
    rs = [float(v) for v in toks[off:off + I - 1]]
    xc = [0.5 * (x[i] + x[i + 1]) for i in range(I - 1)]
    r98 = rho_c + 0.98 * (rho_v - rho_c)
    r02 = rho_c + 0.02 * (rho_v - rho_c)
    n = len(rs) - 1
    if rs[n] >= r98:
        return t, 0.0, 0.0
    i = 0
    while rs[i] > r98:
        i += 1
    if rs[i] == r98:
        vf = xc[i]
    else:
        w = (rs[i] - r98) / (rs[i] - rs[i - 1])
        vf = (1 - w) * xc[i] + w * xc[i - 1]
    cf = xc[n]
    if rs[n] < r02:
        i -= 1
        while rs[i] > r02:
            i += 1
        if rs[i] == r02:
            cf = xc[i]
        else:
            w = (rs[i] - r02) / (rs[i] - rs[i - 1])
            cf = (1 - w) * xc[i] + w * xc[i - 1]
    return t, abs(L - vf), abs(L - cf)
