#!/usr/bin/env python3
"""Exact zero-dimensional oracles for a uniform slab of pyrolysing material (plan 06 verification, 2026-10-08).

Every face adiabatic (BC 301, q = 0), optional uniform volumetric source q [W/m3] from INPUT/st.tec: the state
(T, chi_1..chi_nsub) is spatially uniform and obeys an ODE system built from the case's own 1 K tables
(INPUT/{virgin,char,gas,pyrolysis}.dat, read with the solver's piecewise-linear rule), integrated here with a
fixed-step classical RK4 (pure Python 3.6, no numpy; step 2e-4 s is 1e3 x below the stability limit of the
non-stiff part and gives < 1e-4 K against scipy's RK45 at rtol 1e-11, see the record of 2026-10-08).

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

    # 1 K tables: value and segment slope at T (the solver's rule)
    def _lin(self, Tt, v, T):
        i = int(T) - int(Tt[0])
        return v[i] + (T - Tt[i]) * (v[i + 1] - v[i])

    def _slope(self, Tt, v, T):
        i = int(T) - int(Tt[0])
        return v[i + 1] - v[i]

    def kappa(self, T, r):
        if T < r["T"] or r["A"] <= 0.0:
            return 0.0
        k = r["A"] * math.exp(-r["E"] / (R_UNIV * T))
        if r["n"] != 0.0:
            k *= T ** r["n"]
        return k

    def rhs(self, variant, y):
        T, chi = y[0], y[1:]
        chid = [self.kappa(T, r) * max(1.0 - chi[k], 0.0) ** r["m"] for k, r in enumerate(self.reac)]
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

    def integrate(self, variant, t_end, dt=2.0e-4, t_eval=()):
        """fixed-step RK4 from (T0, chi = 0), landing exactly on every t_eval instant; returns (T, chi) lists there"""
        y = [self.T0] + [0.0] * self.nsub
        t = 0.0
        want = sorted(t_eval)
        out_T, out_chi, iw = [], [], 0
        max_steps = int(t_end / dt) + len(want) + 10
        for step in range(max_steps):
            while iw < len(want) and want[iw] <= t + 1.0e-9:
                out_T.append(y[0]); out_chi.append(list(y[1:])); iw += 1
            if iw >= len(want):
                break
            h = min(dt, want[iw] - t)
            k1 = self.rhs(variant, y)
            k2 = self.rhs(variant, [y[i] + 0.5 * h * k1[i] for i in range(len(y))])
            k3 = self.rhs(variant, [y[i] + 0.5 * h * k2[i] for i in range(len(y))])
            k4 = self.rhs(variant, [y[i] + h * k3[i] for i in range(len(y))])
            y = [y[i] + h / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) for i in range(len(y))]
            t += h
        if iw < len(want):
            raise RuntimeError("integrate: %d of %d instants not reached (t = %g)" % (len(want) - iw, len(want), t))
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
