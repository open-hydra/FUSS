#!/usr/bin/env python3
"""Exact zero-dimensional oracles for a uniform slab of pyrolysing material (plan 06 verification, 2026-10-08).

Every face adiabatic (BC 301, q = 0), optional uniform volumetric source q [W/m3] from INPUT/st.tec: the state
(T, chi_1..chi_nsub) is spatially uniform and obeys an ODE system built from the case's own 1 K tables
(INPUT/{virgin,char,gas,pyrolysis}.dat, read with the solver's piecewise-linear rule), integrated here with
classical RK4 steps of at most dt that each stay on ONE 1 K table segment and keep every reaction's on/off state
(its threshold) fixed: the right-hand side jumps at every integer kelvin (the segment slopes) and at each threshold,
so a step that straddles a jump would be first order; a step that would leave its segment or cross a threshold is
shortened (bisection) to end on it, and the next step takes the next segment from the direction of T. Pure Python 3.6,
no numpy. Self-difference at dt vs dt/2 (the gates' dt 2e-4 s adiabatic, 1e-4 s heated): <= 5.1e-11 K and 8.2e-13 in
chi on the adiabatic cases, <= 1.7e-9 K and 5.2e-13 on the heated ones (they do not fall further at dt/4: the floor
of the event landing), five or more orders below FUSS's error at the committed step.

  chi_r' = kappa_r(T) (1 - chi_r)^m_r,  kappa_r = A exp(-E/(R T)) [T^n] above the threshold T_r, else 0
  pi_tot = sum_r eps_i rho_i F_r chi_r'                                          [kg/m3/s]
  tau    = sum_r w_r (1 - chi_r),  w_r = eps_i rho_i F_r / (rho_v - rho_c)       (virgin fraction)
  hbar   = (rho_v h_v(T) - rho_c h_c(T)) / (rho_v - rho_c)                       (plan 06 section 3.4)
  E1_T   = tau rho_v h_v'(T) + (1 - tau) rho_c h_c'(T)                           (1 K segment slopes)

  type1 (what FUSS integrates, 06 section 2):   E1_T dT/dt = q - pi_tot (h_pg(T) - hbar)
  type2 (PATO's energy, plan 07; gas storage eps(tau) [rho_g h_g - p], rho_g = p W/(R T), the gas leaving at
         h_g(T), outflow mass pi_tot - d(eps rho_g)/dt):
                       (E1_T + eps rho_g c_pg) dT/dt = q - pi_tot (h_pg - hbar) + p d eps/dt
  type1-tablehbar: type1 with hbar from PATO's 14-row TACOT h columns (linear in T) -- the data difference.

API:  case = Case(case_dir)            tables, q, p_gas, T0, the FUSS field series (t, T, chi[], non-uniformity)
      T, chi = case.integrate(variant, t_end, dt=2e-4, t_eval=[...])   -> lists at t_eval
CLI:  ode0d.py <case_dir> [--variant type1|type2|type1-tablehbar] [--dt 2e-4]  prints FUSS - oracle per field file
"""
import math, os, re, sys

R_UNIV = 8.314471469
PATO_TACOT = "/data10/passarani/.conda/envs/pato/src/volume_pato/pato-3.1/data/Materials/Composites/TACOT"


def read_tec(path):
    """Tecplot BLOCK file written by ORION (or by make_case.py): names from the header, cell/nodal split by name."""
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
    return t, d


def _read_solid(path):
    T, cp, h = [], [], []
    for line in open(path):
        p = line.split()
        if len(p) == 8:
            try:
                row = [float(x) for x in p]
            except ValueError:
                continue
            T.append(row[0]); cp.append(row[1]); h.append(row[2])
    return T, cp, h


def _read_gas(path, p_gas):
    T, W, h, zone = [], [], [], None
    for line in open(path):
        if line.startswith("ZONE"):
            zone = line; continue
        p = line.split()
        if len(p) == 5 and zone is not None:
            try:
                row = [float(x) for x in p]
            except ValueError:
                continue
            if abs(row[1] - p_gas) <= 1e-9 * p_gas:
                T.append(row[0]); W.append(row[2]); h.append(row[3])
    if not T:
        raise SystemExit("gas.dat has no zone at p = %g (the oracle needs the exact p-level)" % p_gas)
    return T, W, h


def _read_pyro(path):
    kv = {}
    for line in open(path):
        line = line.split("//", 1)[0].strip()
        if line:
            k, v = line.split(None, 1)
            kv[k] = v.strip()
    nph = int(kv["nSolidPhases"])
    rho_I = [float(kv["rhoI[%d]" % i]) for i in range(1, nph + 1)]
    eps_I = [float(kv["epsI[%d]" % i]) for i in range(1, nph + 1)]
    reac = []
    for i in range(1, nph + 1):
        for j in range(1, int(kv["nPyroReac[%d]" % i]) + 1):
            g = lambda key: float(kv["%s[%d][%d]" % (key, i, j)])
            reac.append(dict(F=g("F"), A=g("A"), E=g("E"), m=g("m"), n=g("n"), T=g("T"),
                             erf=eps_I[i - 1] * rho_I[i - 1] * g("F")))
    rho_v = sum(e * r for e, r in zip(eps_I, rho_I))
    rho_c = rho_v - sum(r["erf"] for r in reac)
    return rho_v, rho_c, float(kv["eps_g_v"]), float(kv["eps_g_c"]), reac


def read_pato_solid(path):
    """PATO virgin/char: rows `p T cp h ki kj kk eps refl`; the first p block."""
    rows = []
    for line in open(path):
        p = line.split()
        if len(p) == 9:
            try:
                rows.append([float(x) for x in p])
            except ValueError:
                continue
    p0 = rows[0][0]
    rows = [r for r in rows if r[0] == p0]
    return [r[1] for r in rows], [r[3] for r in rows]


def _interp(xs, ys, x):
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    lo, hi = 0, len(xs) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if xs[mid] <= x:
            lo = mid
        else:
            hi = mid
    return ys[lo] + (x - xs[lo]) * (ys[hi] - ys[lo]) / (xs[hi] - xs[lo])


class Case(object):
    def __init__(self, case_dir, pato_tables=PATO_TACOT):
        inp = os.path.join(case_dir, "INPUT")
        self.case_dir = case_dir
        self.Tv, self.cpv, self.hv = _read_solid(os.path.join(inp, "virgin.dat"))
        self.Tc, self.cpc, self.hc = _read_solid(os.path.join(inp, "char.dat"))
        self.rho_v, self.rho_c, self.eps_v, self.eps_c, self.reac = _read_pyro(os.path.join(inp, "pyrolysis.dat"))
        self.p_gas = 101325.0
        for line in open(os.path.join(case_dir, "input.ini")):
            if line.strip().startswith("p-gas"):
                self.p_gas = float(line.split("=")[1].split(";")[0])
        self.Tg, self.Wg, self.hg = _read_gas(os.path.join(inp, "gas.dat"), self.p_gas)
        self.q = 0.0
        st = os.path.join(inp, "st.tec")
        if os.path.exists(st):
            _, dst = read_tec(st)
            qv = dst["qvol"]
            if max(qv) - min(qv) != 0.0:
                raise SystemExit("st.tec is not uniform: the 0-D oracle does not apply")
            self.q = qv[0]
        _, dic = read_tec(os.path.join(inp, "ic.tec"))
        self.T0 = dic["T"][0]
        self.nsub = len(self.reac)
        self.w = [r["erf"] / (self.rho_v - self.rho_c) for r in self.reac]
        self.pTv = self.phv = self.pTc = self.phc = None
        if pato_tables and os.path.isdir(pato_tables):
            self.pTv, self.phv = read_pato_solid(os.path.join(pato_tables, "virgin"))
            self.pTc, self.phc = read_pato_solid(os.path.join(pato_tables, "char"))

    _seg = None      # integer kelvin of the frozen segment [seg, seg + 1] of a step, or None (the solver's rule int(T))
    _kon = None      # frozen on/off state per reaction index during a step, or None (T >= threshold)

    # 1 K tables: value and segment slope at T (the solver's rule; the frozen segment inside an integration step)
    def _lin(self, Tt, v, T):
        i = (self._seg if self._seg is not None else int(T)) - int(Tt[0])
        return v[i] + (T - Tt[i]) * (v[i + 1] - v[i])

    def _slope(self, Tt, v, T):
        i = (self._seg if self._seg is not None else int(T)) - int(Tt[0])
        return v[i + 1] - v[i]

    def kappa(self, T, r, k=None):
        on = (T >= r["T"]) if (self._kon is None or k is None) else self._kon[k]
        if not on or r["A"] <= 0.0:
            return 0.0
        k = r["A"] * math.exp(-r["E"] / (R_UNIV * T))
        if r["n"] != 0.0:
            k *= T ** r["n"]
        return k

    def rhs(self, variant, y):
        T, chi = y[0], y[1:]
        chid = [self.kappa(T, r, k) * max(1.0 - chi[k], 0.0) ** r["m"] for k, r in enumerate(self.reac)]
        pi = sum(r["erf"] * chid[k] for k, r in enumerate(self.reac))
        tau = min(max(sum(self.w[k] * (1.0 - chi[k]) for k in range(self.nsub)), 0.0), 1.0)
        if variant == "type1-tablehbar":
            hbar = (self.rho_v * _interp(self.pTv, self.phv, T) - self.rho_c * _interp(self.pTc, self.phc, T)) / (self.rho_v - self.rho_c)
        else:
            hbar = (self.rho_v * self._lin(self.Tv, self.hv, T) - self.rho_c * self._lin(self.Tc, self.hc, T)) / (self.rho_v - self.rho_c)
        E1T = tau * self.rho_v * self._slope(self.Tv, self.hv, T) + (1.0 - tau) * self.rho_c * self._slope(self.Tc, self.hc, T)
        hpg = self._lin(self.Tg, self.hg, T)
        if variant == "type2":
            W = self._lin(self.Tg, self.Wg, T)
            rho_g = self.p_gas * W / (R_UNIV * T)
            cpg = self._slope(self.Tg, self.hg, T)
            eps = self.eps_c + (self.eps_v - self.eps_c) * tau
            epsd = (self.eps_v - self.eps_c) * (-sum(self.w[k] * chid[k] for k in range(self.nsub)))
            dT = (self.q - pi * (hpg - hbar) + self.p_gas * epsd) / (E1T + eps * rho_g * cpg)
        else:
            dT = (self.q - pi * (hpg - hbar)) / E1T
        return [dT] + chid

    def _rk4(self, variant, y, h):
        k1 = self.rhs(variant, y)
        k2 = self.rhs(variant, [y[i] + 0.5 * h * k1[i] for i in range(len(y))])
        k3 = self.rhs(variant, [y[i] + 0.5 * h * k2[i] for i in range(len(y))])
        k4 = self.rhs(variant, [y[i] + h * k3[i] for i in range(len(y))])
        return [y[i] + h / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) for i in range(len(y))]

    def integrate(self, variant, t_end, dt=2.0e-4, t_eval=()):
        """RK4 from (T0, chi = 0), steps of at most dt on one table segment and one on/off state of every reaction
        (module docstring), landing exactly on every t_eval instant; returns (T, chi) lists there. t_end is kept
        for the callers; the last t_eval instant ends the integration. self.nevents counts the shortened steps."""
        y = [self.T0] + [0.0] * self.nsub
        t = 0.0
        want = sorted(t_eval)
        out_T, out_chi, iw = [], [], 0
        self.nevents = 0
        try:
            while True:
                while iw < len(want) and want[iw] <= t + 1.0e-9:
                    out_T.append(y[0]); out_chi.append(list(y[1:])); iw += 1
                if iw >= len(want):
                    break
                h = min(dt, want[iw] - t)
                self._seg, self._kon = None, None
                up = self.rhs(variant, y)[0] >= 0.0
                T0 = y[0]
                self._seg = int(math.floor(T0)) if up else int(math.ceil(T0)) - 1
                self._kon = [T0 >= r["T"] if up else T0 > r["T"] for r in self.reac]
                # the nearest event in the direction of motion: the segment end or a reaction threshold
                ends = [self._seg + 1.0] if up else [float(self._seg)]
                ends += [r["T"] for r in self.reac if (up and T0 < r["T"]) or (not up and T0 > r["T"])]
                e = min(ends) if up else max(ends)
                yn = self._rk4(variant, y, h)
                if (up and yn[0] > e) or (not up and yn[0] < e):
                    a, b = 0.0, h
                    for _ in range(60):
                        m = 0.5 * (a + b)
                        Tm = self._rk4(variant, y, m)[0]
                        if (up and Tm < e) or (not up and Tm > e):
                            a = m
                        else:
                            b = m
                    h = b
                    yn = self._rk4(variant, y, h)
                    yn[0] = e                    # on the event: the next step picks the next segment from its direction
                    self.nevents += 1
                y = yn
                t += h
        finally:
            self._seg, self._kon = None, None
        return out_T, out_chi

    def fuss_series(self):
        """(t, T, [chi_r], non-uniformity) per OUTPUT/field<i>.tec, sorted in time"""
        outdir = os.path.join(self.case_dir, "OUTPUT")
        files = [f for f in os.listdir(outdir) if re.match(r"field\d*\.tec$", f)]
        series = []
        for f in files:
            t, d = read_tec(os.path.join(outdir, f))
            cols = [d["T"]] + [d["chi_%d" % (k + 1)] for k in range(self.nsub)]
            unif = max(max(c) - min(c) for c in cols)
            series.append((t, d["T"][0], [d["chi_%d" % (k + 1)][0] for k in range(self.nsub)], unif))
        series.sort()
        return series

    def fuss_fields(self):
        """per OUTPUT/field<i>.tec, sorted in time: dict t, T and chi (cell 1), the spatial spreads unif_T and unif_chi
        (per reaction), and dp, the spread of p_pg over the cells (0 without a p_pg field)"""
        outdir = os.path.join(self.case_dir, "OUTPUT")
        out = []
        for f in os.listdir(outdir):
            if not re.match(r"field\d*\.tec$", f):
                continue
            t, d = read_tec(os.path.join(outdir, f))
            chis = [d["chi_%d" % (k + 1)] for k in range(self.nsub)]
            p = d.get("p_pg") or [0.0]
            out.append({"t": t, "T": d["T"][0], "chi": [c[0] for c in chis], "unif_T": max(d["T"]) - min(d["T"]),
                        "unif_chi": [max(c) - min(c) for c in chis], "dp": max(p) - min(p)})
        out.sort(key=lambda r: r["t"])
        return out

    def heat_capacity(self, variant, T, chi):
        """(eps, C): the gas porosity and dE/dT [J/m3/K] of the uniform state, with the gas storage under type2"""
        tau = min(max(sum(self.w[k] * (1.0 - chi[k]) for k in range(self.nsub)), 0.0), 1.0)
        C = tau * self.rho_v * self._slope(self.Tv, self.hv, T) + (1.0 - tau) * self.rho_c * self._slope(self.Tc, self.hc, T)
        eps = self.eps_c + (self.eps_v - self.eps_c) * tau
        if variant == "type2":
            rho_g = self.p_gas * self._lin(self.Tg, self.Wg, T) / (R_UNIV * T)
            C += eps * rho_g * self._slope(self.Tg, self.hg, T)
        return eps, C


EPS_MACHINE = 2.0 ** -52
TOL_UNIF_TYPE2 = 1.0e-5     # K (and chi): type-2 non-uniformity, a floor from measured draws (check_ladder)


def input_value(case_dir, key, default=None):
    """the value of `key = value` in the case's input.ini (comments after ';' dropped)"""
    for line in open(os.path.join(case_dir, "input.ini")):
        s = line.split(";")[0].strip()
        if "=" in s and s.split("=", 1)[0].strip() == key:
            return s.split("=", 1)[1].strip()
    return default


def n_steps(case_dir):
    """the run's step count: the lines of OUTPUT/residual-history.dat (one per step, no header)"""
    with open(os.path.join(case_dir, "OUTPUT", "residual-history.dat")) as fh:
        return sum(1 for line in fh if line.strip())


def check_ladder(name, cases, dt, pato_tables=None):
    """The gate of the 0-D cases (adiabatic0d, heated0d and their type-2 twins): FUSS against the oracle on a ladder of
    runs of one case at vnn, vnn/2, vnn/4 (coarsest first; refinement = 1/vnn, the step scales with vnn).
      error  e_T = max_t |T - T_or|, e_chi = max_t,r |chi_r - chi_or,r| per run. Type 1: the convergence rule
             (tools/conv_rule.py) with p = 2 (the kinetics/energy splitting inside the RK3 stages). Type 2: the error
             must fall on every rung; the observed order on the finest pair is printed, not gated (no expected order
             is stated for the type-2 step: the pressure step is first order in time, the splitting second order).
      non-uniformity of the uniform state, per run. Type 1: every cell sees the same arithmetic but the 301 faces'
             ghost values, one rounding per RK3 stage per step and no amplification (conduction contracts a
             difference, the endothermic source damps it): T spread <= N_steps x 3 x eps x max T, chi spread <=
             N_steps x 3 x eps. Type 2: T and chi spread <= TOL_UNIF_TYPE2, a floor set from measured draws
             (1.8e-7 K adiabatic, 6.7e-7 K heated), not a physics budget. Printed beside it, not gated: the candidate
             argument "E carries -eps_g p, so cells of equal energy whose pressures differ by dp differ in T by
             eps_g dp / C (C = dE/dT)", max over the files with dp the run's own cell spread of p_pg; it ignores the
             flux terms, and the measured T spread exceeds it on both type-2 cases (so it is not adopted).
    The oracle is integrated once per variant on the union of the runs' field-file times (each run writes where its
    own step crosses sol-dtime, so the times differ from rung to rung). Returns True on PASS."""
    import conv_rule as cr
    closures = [input_value(c, "gas-closure", "type1") for c in cases]
    if len(set(closures)) != 1:
        print("%s: FAIL  the runs mix gas closures %s" % (name, closures))
        return False
    variant = "type2" if closures[0] == "type2" else "type1"
    objs = [Case(c, pato_tables=pato_tables) for c in cases]
    o = objs[0]
    for c in objs[1:]:
        if (c.T0, c.q, c.nsub, c.p_gas) != (o.T0, o.q, o.nsub, o.p_gas):
            print("%s: FAIL  %s is not a rung of %s (T0, q, reactions or p_gas differ)" % (name, c.case_dir, o.case_dir))
            return False
    vnn = [float(input_value(c, "vnn")) for c in cases]
    recs = [c.fuss_fields() for c in objs]
    times = sorted(set(r["t"] for rr in recs for r in rr))
    want = ["type1", "type2"] + (["type1-tablehbar"] if o.pTv is not None else [])
    orc = {}
    for v in want:
        T, chi = o.integrate(v, times[-1], dt, times)
        orc[v] = dict(zip(times, zip(T, chi)))
    tag = variant[-1]
    print("%s (%s, oracle %s): T0 = %g K, q = %g W/m3, ladder vnn %s" % (name, closures[0], variant, o.T0, o.q,
                                                                         " ".join("%g" % v for v in vnn)))
    ok = True
    eTs, ecs = [], []
    for c, rr, v in zip(cases, recs, vnn):
        eT = max(abs(r["T"] - orc[variant][r["t"]][0]) for r in rr)
        ec = max(abs(r["chi"][k] - orc[variant][r["t"]][1][k]) for r in rr for k in range(o.nsub))
        eTs.append(eT); ecs.append(ec)
        print("  vnn %-9g %3d files to %.3f s, %6d steps: max |T - T_or%s| %.4e K  max |chi - chi_or%s| %.3e  T_end FUSS %.3f"
              % (v, len(rr), rr[-1]["t"], n_steps(c), tag, eT, tag, ec, rr[-1]["T"]))
    if len(cases) < 2:
        print("  FAIL  the convergence rule needs a ladder of 2 or more runs (vnn, vnn/2, vnn/4)")
        ok = False
    else:
        refine = [1.0 / v for v in vnn]
        for label, errs, fmt in (("max |T - T_or%s| [K]" % tag, eTs, "%.4e"), ("max |chi - chi_or%s|" % tag, ecs, "%.3e")):
            if variant == "type1":
                good, text = cr.check(label, errs, refine, 2, fmt=fmt)
            else:
                good, text = cr.check_falls(label, errs, refine, fmt=fmt)
            ok &= good
            print("  %s  %s" % (text, "PASS" if good else "FAIL"))
    sel = cr.selftest()
    for c, rr, v in zip(cases, recs, vnn):
        uT = max(r["unif_T"] for r in rr)
        uc = [max(r["unif_chi"][k] for r in rr) for k in range(o.nsub)]
        if variant == "type1":
            N = n_steps(c)
            bT = N * 3 * EPS_MACHINE * max(r["T"] for r in rr)
            bc = [N * 3 * EPS_MACHINE] * o.nsub
            how = "N_steps %d x 3 x eps x max T" % N
        else:
            cand = 0.0
            for r in rr:
                eps, C = o.heat_capacity("type2", r["T"], r["chi"])
                cand = max(cand, eps * r["dp"] / C)
            bT, bc = TOL_UNIF_TYPE2, [TOL_UNIF_TYPE2] * o.nsub
            how = "floor from measured draws; candidate eps_g dp / C = %.2e K, dp up to %.3e Pa, %s" % (
                cand, max(r["dp"] for r in rr), "it would FAIL" if uT > cand else "it would pass")
        if sel:
            bT, bc = -1.0, [-1.0] * o.nsub
        good = uT <= bT and all(u <= b for u, b in zip(uc, bc))
        ok &= good
        print("  vnn %-9g non-uniformity T %.2e K (bound %.2e: %s); chi %s (bounds %s)  %s"
              % (v, uT, bT, how, " ".join("%.1e" % u for u in uc), " ".join("%.1e" % b for b in bc),
                 "PASS" if good else "FAIL"))
    t1 = recs[0]
    d2 = max((orc["type2"][r["t"]][0] - orc["type1"][r["t"]][0] for r in t1), key=abs)
    line = "  info (vnn %g): type2 - type1 (exact) max %.4f K" % (vnn[0], d2)
    if "type1-tablehbar" in orc:
        line += "; PATO-table hbar minus type1 max %.3f K" % max(
            (orc["type1-tablehbar"][r["t"]][0] - orc["type1"][r["t"]][0] for r in t1), key=abs)
    print(line)
    print("%s: %s" % (name, "PASS" if ok else "FAIL"))
    return ok


def main():
    args = sys.argv[1:]
    variant, dt, cases = "type1", 2.0e-4, []
    while args:
        a = args.pop(0)
        if a == "--variant":
            variant = args.pop(0)
        elif a == "--dt":
            dt = float(args.pop(0))
        else:
            cases.append(a)
    for cdir in cases:
        c = Case(cdir)
        s = c.fuss_series()
        ts = [x[0] for x in s]
        T, chi = c.integrate(variant, ts[-1], dt, ts)
        eT = [s[i][1] - T[i] for i in range(len(s))]
        ech = [max(abs(s[i][2][k] - chi[i][k]) for k in range(c.nsub)) for i in range(len(s))]
        print("== %s: variant %s, q = %g W/m3, T0 = %g K, %d files to %.3f s, non-uniformity %.1e" % (cdir, variant, c.q, c.T0, len(s), ts[-1], max(x[3] for x in s)))
        print("   max |T_FUSS - oracle| = %.4f K (t = %.2f)   end %.4f K   max |chi - oracle| = %.2e" % (max(abs(e) for e in eT), ts[eT.index(max(eT, key=abs))], eT[-1], max(ech)))


if __name__ == "__main__":
    main()
