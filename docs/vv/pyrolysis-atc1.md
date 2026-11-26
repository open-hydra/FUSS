# Ablation Test Case 1

Ablation Test Case 1 (ATC1) couples conduction, pyrolysis kinetics and the type-1 gas-enthalpy
advection on a real material (TACOT) over the full decomposition range. Because no analytic solution
exists for this coupled problem, the case is gated on FUSS's own Richardson extrapolation (mesh
convergence of its own model), and compared, for information only, with the independently computed
FIAT solution bundled with PATO 3.1.

Setup and run: `test/pyrolysis/ATC1/README.md`.

---

## Problem

A 0.05 m TACOT slab, 100 cells graded 0.1 toward the heated face (face 2); face 2 follows
`Tw_time.dat` (298 → 1644 K in 0.1 s, held to 60 s); face 1 is adiabatic (BC 301); $T(0) = 298$ K (the
ATC1 specification, [ATC1] v1.1; PATO's own tutorial starts at 300 K). RK3, `cons`, `vnn = 0.25`,
type-1 gas closure, real TACOT tables (the same reaction set as
[Pyrolysis Kinetics and 0-D Decomposition](pyrolysis-kinetics-0d.md)). FUSS brackets each of FIAT's
six thermocouple depths (1, 2, 4, 8, 16, 50 mm) with two index-position probes whose cell centres
straddle the exact depth; the reported value is the linear interpolation between them, then linear in
time.

## Reference

**Gate (anchor).** FUSS's own Richardson extrapolation on 100/200/400 graded cells, at FIAT's instants
$t=1\ldots60$ s and the six depths, where the three meshes are asymptotic (observed order $p$ between
0.5 and 3; outside that window, or on a non-monotone sequence at the 0.01 K level, the 400-cell value
is used instead and its distance to the limit is bounded by $f_{200}-f_{400}\le0.1$ K at those points).
The gate checks that the 100-cell discretisation error, $f_{100}-f_\star$, has not moved by more than
`band_K` (1.0 K) since it was measured: for every point,

$$
\bigl\lvert (T_{\rm FUSS} - f_\star) - (f_{100}-f_\star)\bigr\rvert \le \text{band\_K}.
$$

**Informational.** FIAT, PATO 3.1's own reference solution for ATC1 (`ref/FIAT_T.dat`, columns
$t, T(0), T(1), T(2), T(4), T(8), T(16), T(50\text{ mm})$), and PATO 3.1 itself run on the unmodified
tutorial with the initial and first ramp value at 298 K (`ref/pato-298K`).

## Results

A sample of the anchor (`expected.json`, computed with the 302 ramp read at the RK stage times and
the type-1 gas enthalpy evaluated at the face temperature): $f_\star$, the 100-cell discretisation
error $f_{100}-f_\star$, and the observed order $p$, at the five times FIAT reports and every depth:

| Depth | $t=1$ s | $t=5$ s | $t=10$ s | $t=30$ s | $t=60$ s |
|---|---|---|---|---|---|
| 1 mm | $f_\star=966.992$, $f_{100}-f_\star=-1.868$, $p=1.54$ | $1298.075$, $-0.327$, $p=1.15$ | $1388.035$, $-0.160$, $p=1.0$ | $1486.395$, $-0.048$, $p=0.87$ | $1528.133$, $-0.010$, $p=-0.58$ |
| 2 mm | — | $1019.720$, $-0.440$, $p=3.49$ | $1170.718$, $-0.193$, $p=3.2$ | $1344.064$, $-0.052$ | $1421.196$, $-0.018$ |
| 4 mm | $314.393$, $0.255$, $p=2.17$ | $631.763$, $-0.243$, $p=1.79$ | $831.355$, $-0.103$, $p=1.69$ | $1099.606$, $-0.018$, $p=3.0$ | $1230.052$, $0.024$, $p=2.25$ |
| 8 mm | $298.002$, $0.000$ | $333.455$, $0.152$, $p=2.6$ | $445.483$, $0.071$, $p=4.5$ | $749.621$, $0.003$ | $920.091$, $0.093$, $p=2.86$ |
| 16 mm | — | $298.026$, $0.003$ | $300.938$, $0.051$, $p=1.89$ | $388.427$, $0.119$, $p=1.73$ | $555.785$, $0.038$, $p=1.42$ |
| 50 mm | — | — | — | $298.0$, $0.0$ | $298.305$, $0.012$, $p=2.17$ |

Over the full 46-point window the largest 100-cell discretisation error is 1.868 K, at (1 mm, 1 s);
the per-point order $p$ (shown above where defined) ranges from $-0.58$ to $6.46$ across the window —
a locally fitted exponent at each (depth, $t$), not evidence of one asymptotic rate. $f_\star$ moved by
only 0.007 K at (1 mm, 1 s) under a refinement of the scheme that reaches it (the 302-ramp timing and
the type-1 gas-enthalpy evaluation point): $f_\star$ is a property of the model, not of the scheme.

On the current build the gate's largest move is 0.443 K, at (16 mm, 60 s), well inside the 1.0 K band.

Reported, not gated: the metric against FIAT directly, max $\lvert T_{\rm FUSS}-T_{\rm FIAT}\rvert$
over the six probes and the full window, is 16.03 K (at 2 mm, 4 s); FIAT runs 4-17 K cooler than the
converged FUSS model during each front passage, within $-0.1$ to $+0.7$ K at 1 mm after 20 s; the
pyrolysis front depths ($\rho_s$ 98 %/2 % thresholds between $\rho_c=220$ and $\rho_v=280$ kg/m³,
15.21 mm and 2.91 mm below the surface at the last field write) and the surface $\dot m_{pg}$ there,
$5.387\times10^{-3}$ kg/(m²s).

## Pass criteria

For every one of the 46 gated (depth, $t$) points: $\lvert(T_{\rm FUSS}-f_\star)-(f_{100}-f_\star)\rvert\le1.0$ K,
using the committed $f_\star$ and $f_{100}-f_\star$ of `expected.json`. A change to the model (type 2,
the ablating wall) or to the mesh/probe layout invalidates the anchor and needs a new Richardson
extrapolation.
