"""Exact solutions for the type-2 (Darcy) closure of FUSS (exact-solution benchmarks, 2026-10-10).

Three problems, each with a closed form (a series or a similarity solution with one transcendental
equation), chosen for the parts of the type-2 model that the plan-07 gates check only at steady
state or not at all:

1. pressure_transient_linear: transient gas pressure in a porous slab, small amplitude.
   eps W/(R T) dp/dt = d/dx( p W K/(R T mu) dp/dx ) at uniform constant T; for |p - p0| << p0 this is
   dp/dt = D d2p/dx2 with D = K p0/(eps mu). Slab 0 < x < L, sealed at x = 0 (dp/dx = 0), p = p_w
   at x = L for t > 0, p = p0 at t = 0. Exact (separation of variables):
       (p - p_w)/(p0 - p_w) = sum_n 4 (-1)^n/((2n+1) pi) cos((2n+1) pi x/(2L)) exp(-(2n+1)^2 pi^2 D t/(4 L^2)).
   The linearisation error is O(dp/p0) relative to dp.

2. transpiration_steady: steady conduction with gas transpiration through a porous slab.
   Both faces open (p_in at x = 0, p_out at x = L), Robin heat exchange h_i (Tref_i - T_wall) at both
   faces, constant k and gas cp: k T'' - mdot cp_g T' = 0 with the uniform mass flux mdot (+x):
       T(x) = A + B exp(Pe x/L),  Pe = mdot cp_g L/k,
   A, B from -k T'(0) = h1 (Tref1 - T(0)) and k T'(L) = h2 (Tref2 - T(L)). mdot follows from steady
   Darcy with rho_g = p W/(R T): d(p^2)/dx = -2 mdot R mu T(x)/(W K), so
       mdot = (p_in^2 - p_out^2) W K / (2 R mu int_0^L T dx),
   a fixed point in mdot (T depends on mdot through Pe), solved here to 1e-15.

3. sharp_front_similarity: decomposition front with transpiration through the char (Neumann-type).
   Semi-infinite virgin slab at T0, surface held at Ts > Tp from t = 0, the solid decomposes
   instantly at the isotherm T = Tp (front x_f = 2 lam sqrt(alpha_c t)); the gas produced at the front,
   mdot = (rho_v - rho_c) dx_f/dt, crosses the char to the surface (no storage, impermeable virgin).
   Char:   rho_c cp_c T_t = k_c T_xx + mdot cp_g T_x   ->  T = Ts + (Tp - Ts) [erf(eta+beta) - erf(beta)]
                                                                        / [erf(lam+beta) - erf(beta)]
   Virgin: rho_v cp_v T_t = k_v T_xx                    ->  T = T0 + (Tp - T0) erfc(eta r) / erfc(lam r)
   with eta = x/(2 sqrt(alpha_c t)), beta = lam (rho_v - rho_c) cp_g/(rho_c cp_c), r = sqrt(alpha_c/alpha_v).
   Front (Rankine-Hugoniot of the total energy): k_c T_x(x_f-) - k_v T_x(x_f+) = - V Q_endo,
   V = dx_f/dt, Q_endo = rho_c h_c(Tp) + (rho_v - rho_c) h_g(Tp) - rho_v h_v(Tp) (>0 endothermic),
   i.e.  k_c T'_c(lam) - k_v T'_v(lam) = - 2 lam alpha_c Q_endo  (derivatives in eta); one equation
   for lam, solved by bisection. The pressure in the char (quasi-steady Darcy, rho_g = p W/(R T)):
   p(x)^2 = p_w^2 + (2 mdot R mu/(W K_c)) int_0^x T dx'.

Pure Python 3.6 (math only), no numpy, so every verify.py can import it.
"""
import math

R_UNIV = 8.314471469     # J/mol/K, FUSS_Parameters_m


# ----------------------------------------------------------------------------- 1. transient pressure
def pressure_transient_linear(x, t, L, D, p0, pw, nterms=4000):
    """p(x, t) of problem 1 (linear limit)."""
    if t <= 0.0:
        return pw if x >= L else p0
    s = 0.0
    for n in range(nterms):
        k = 2 * n + 1
        decay = math.exp(-k * k * math.pi * math.pi * D * t / (4.0 * L * L))
        if decay < 1.0e-18:
            break
        s += 4.0 * (-1) ** n / (k * math.pi) * math.cos(k * math.pi * x / (2.0 * L)) * decay
    return pw + (p0 - pw) * s


# ----------------------------------------------------------------------------- 2. steady transpiration
def transpiration_profile(L, k, cpg, mdot, h1, Tref1, h2, Tref2):
    """A, B, Pe of T(x) = A + B exp(Pe x/L) for a given mdot (Robin conditions at both faces)."""
    Pe = mdot * cpg * L / k
    if abs(Pe) < 1.0e-12:
        # pure conduction limit: T = A + B x/L
        # -k B/L = h1 (Tref1 - A);  k B/L = h2 (Tref2 - A - B)
        a11, a12, b1 = h1, -k / L, h1 * Tref1
        a21, a22, b2 = h2, h2 + k / L, h2 * Tref2
        det = a11 * a22 - a12 * a21
        A = (b1 * a22 - a12 * b2) / det
        B = (a11 * b2 - a21 * b1) / det
        return A, B, 0.0, True
    E = math.exp(Pe)
    # -k B Pe/L = h1 (Tref1 - A - B)          ->  h1 A + (h1 - k Pe/L) B = h1 Tref1
    #  k B Pe/L E = h2 (Tref2 - A - B E)      ->  h2 A + (h2 E + k Pe E/L) B = h2 Tref2
    a11, a12, b1 = h1, h1 - k * Pe / L, h1 * Tref1
    a21, a22, b2 = h2, (h2 + k * Pe / L) * E, h2 * Tref2
    det = a11 * a22 - a12 * a21
    A = (b1 * a22 - a12 * b2) / det
    B = (a11 * b2 - a21 * b1) / det
    return A, B, Pe, False


def transpiration_T(x, L, A, B, Pe, linear):
    return A + B * x / L if linear else A + B * math.exp(Pe * x / L)


def transpiration_int_T(L, A, B, Pe, linear):
    """int_0^L T dx."""
    if linear:
        return A * L + 0.5 * B * L
    return A * L + B * L / Pe * (math.exp(Pe) - 1.0)


def transpiration_steady(L, k, cpg, W, mu, K, p_in, p_out, h1, Tref1, h2, Tref2):
    """mdot (kg/m2/s, +x), A, B, Pe of the steady transpiration problem (fixed point in mdot)."""
    mdot = 0.0
    for _ in range(200):
        A, B, Pe, lin = transpiration_profile(L, k, cpg, mdot, h1, Tref1, h2, Tref2)
        I = transpiration_int_T(L, A, B, Pe, lin)
        new = (p_in * p_in - p_out * p_out) * W * K / (2.0 * R_UNIV * mu * I)
        if abs(new - mdot) <= 1.0e-15 * abs(new):
            mdot = new
            break
        mdot = new
    A, B, Pe, lin = transpiration_profile(L, k, cpg, mdot, h1, Tref1, h2, Tref2)
    return mdot, A, B, Pe, lin


def transpiration_p(x, L, mdot, W, mu, K, p_in, A, B, Pe, linear):
    """p(x) from d(p^2)/dx = -2 mdot R mu T/(W K)."""
    if linear:
        I = A * x + 0.5 * B * x * x / L
    else:
        I = A * x + B * L / Pe * (math.exp(Pe * x / L) - 1.0)
    return math.sqrt(p_in * p_in - 2.0 * mdot * R_UNIV * mu * I / (W * K))


# ----------------------------------------------------------------------------- 3. sharp-front similarity
def _front_residual(lam, ac, av, kc, kv, rho_v, rho_c, cpc, cpg, Ts, Tp, T0, Q):
    beta = lam * (rho_v - rho_c) * cpg / (rho_c * cpc)
    r = math.sqrt(ac / av)
    den_c = math.erf(lam + beta) - math.erf(beta)
    dTc = (Tp - Ts) * 2.0 / math.sqrt(math.pi) * math.exp(-(lam + beta) ** 2) / den_c
    dTv = -(Tp - T0) * 2.0 / math.sqrt(math.pi) * r * math.exp(-(lam * r) ** 2) / math.erfc(lam * r)
    return kc * dTc - kv * dTv + 2.0 * lam * ac * Q


def sharp_front_lambda(ac, av, kc, kv, rho_v, rho_c, cpc, cpg, Ts, Tp, T0, Q):
    """lam of x_f = 2 lam sqrt(alpha_c t), by bisection on the front condition."""
    f = lambda lam: _front_residual(lam, ac, av, kc, kv, rho_v, rho_c, cpc, cpg, Ts, Tp, T0, Q)
    lo, hi = 1.0e-9, 1.0
    while f(hi) < 0.0:
        hi *= 2.0
        if hi > 100.0:
            raise ValueError("sharp_front_lambda: no bracket")
    if f(lo) > 0.0:
        raise ValueError("sharp_front_lambda: no bracket at the low end")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0.0:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1.0e-15 * hi:
            break
    return 0.5 * (lo + hi)


def sharp_front_T(x, t, lam, ac, av, rho_v, rho_c, cpc, cpg, Ts, Tp, T0):
    """T(x, t) of problem 3."""
    if t <= 0.0:
        return Ts if x <= 0.0 else T0
    eta = x / (2.0 * math.sqrt(ac * t))
    beta = lam * (rho_v - rho_c) * cpg / (rho_c * cpc)
    r = math.sqrt(ac / av)
    if eta <= lam:
        return Ts + (Tp - Ts) * (math.erf(eta + beta) - math.erf(beta)) / (math.erf(lam + beta) - math.erf(beta))
    return T0 + (Tp - T0) * math.erfc(eta * r) / math.erfc(lam * r)


def sharp_front_position(t, lam, ac):
    return 2.0 * lam * math.sqrt(ac * t)


def sharp_front_mdot(t, lam, ac, rho_v, rho_c):
    """Gas mass flux through the char toward the surface, kg/m2/s (magnitude)."""
    return (rho_v - rho_c) * lam * math.sqrt(ac / t)


def sharp_front_p(x, t, lam, ac, av, rho_v, rho_c, cpc, cpg, Ts, Tp, T0, W, mu, Kc, pw, nint=2000):
    """Char pressure (quasi-steady Darcy); in the virgin (impermeable) the front value is returned."""
    xf = sharp_front_position(t, lam, ac)
    xe = min(x, xf)
    md = sharp_front_mdot(t, lam, ac, rho_v, rho_c)
    # Simpson on [0, xe]
    n = max(2, nint - nint % 2)
    h = xe / n
    s = 0.0
    for i in range(n + 1):
        w = 1.0 if i in (0, n) else (4.0 if i % 2 else 2.0)
        s += w * sharp_front_T(i * h, t, lam, ac, av, rho_v, rho_c, cpc, cpg, Ts, Tp, T0)
    I = s * h / 3.0
    return math.sqrt(pw * pw + 2.0 * md * R_UNIV * mu * I / (W * Kc))


if __name__ == "__main__":
    # self-checks: series at t -> large, conduction limit of transpiration, PDE residuals of 3
    L, D = 0.05, 1.0e-6
    assert abs(pressure_transient_linear(0.0, 1.0e9, L, D, 101325.0, 101326.0) - 101326.0) < 1e-9
    # early time (D t/L^2 = 1e-4): the front has not reached x = 0, p(0) = p0 to the series' accuracy
    assert abs(pressure_transient_linear(0.0, 0.25, L, D, 101325.0, 101326.0) - 101325.0) < 1e-9
    # half-way check against the one-term asymptote at late time (D t/L^2 = 1): higher modes are e^-22 smaller
    late = pressure_transient_linear(0.0, 2500.0, L, D, 101325.0, 101326.0)
    one = 101326.0 - 4.0 / math.pi * math.exp(-math.pi ** 2 / 4.0)
    assert abs(late - one) < 1e-9, (late, one)
    md, A, B, Pe, lin = transpiration_steady(0.05, 0.4, 1000.0, 0.022, 3e-5, 4.7e-13, 2 * 101325.0, 101325.0,
                                             50.0, 300.0, 50.0, 1000.0)
    k, L = 0.4, 0.05
    T = lambda x: transpiration_T(x, L, A, B, Pe, lin)
    dT = lambda x: B * Pe / L * math.exp(Pe * x / L)          # analytic T'(x)
    print("transpiration: mdot %.6e  Pe %.4f  T(0) %.4f  T(L) %.4f" % (md, Pe, T(0), T(L)))
    r0 = -k * dT(0.0) - 50 * (300 - T(0))
    rL = k * dT(L) - 50 * (1000 - T(L))
    print("  Robin residuals at the faces %.2e %.2e W/m2" % (r0, rL))
    assert abs(r0) < 1e-9 and abs(rL) < 1e-9
    I = transpiration_int_T(L, A, B, Pe, lin)
    assert abs((2 * 101325.0) ** 2 - 101325.0 ** 2 - 2 * md * R_UNIV * 3e-5 * I / (0.022 * 4.7e-13)) < 1e-6 * 101325.0 ** 2
    ac, av, kc, kv = 0.4 / (220 * 1000.0), 0.4 / (280 * 1000.0), 0.4, 0.4
    lam = sharp_front_lambda(ac, av, kc, kv, 280.0, 220.0, 1000.0, 1000.0, 1200.0, 600.0, 300.0, 0.0)
    print("sharp front: lam %.10f  x_f(10 s) %.6e m" % (lam, sharp_front_position(10.0, lam, ac)))
    # PDE residual of the char solution: rho_c cp_c T_t - k_c T_xx - mdot cp_g T_x
    t, x, e = 10.0, 0.5 * sharp_front_position(10.0, lam, ac), 1e-6
    f = lambda x, t: sharp_front_T(x, t, lam, ac, av, 280.0, 220.0, 1000.0, 1000.0, 1200.0, 600.0, 300.0)
    Tt = (f(x, t + e) - f(x, t - e)) / (2 * e)
    Tx = (f(x + e, t) - f(x - e, t)) / (2 * e)
    Txx = (f(x + e, t) - 2 * f(x, t) + f(x - e, t)) / (e * e)
    md = sharp_front_mdot(t, lam, ac, 280.0, 220.0)
    print("  char PDE residual / (rho c T_t): %.2e" % ((220 * 1000 * Tt - 0.4 * Txx - md * 1000 * Tx) / (220 * 1000 * Tt)))
