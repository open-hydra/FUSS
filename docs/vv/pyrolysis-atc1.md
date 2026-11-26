# Ablation Test Case 1

Ablation Test Case 1 (ATC1) couples conduction, pyrolysis kinetics and the type-1 gas-enthalpy
advection on a real material (TACOT) over the full decomposition range. No analytic solution exists for
this coupled problem: FUSS is compared with PATO 3.1 run on the same problem, and the tolerance is the
disagreement between PATO and FIAT, the second independent code whose ATC1 solution ships with PATO.

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

**Oracle.** PATO 3.1 run on the unmodified ATC1 tutorial with the initial and first ramp value at 298 K
(`ref/pato-298K`: `Ta_plot` at the five thermocouple depths, `Ta_surfacePatch` for the back face).

**Tolerance.** FIAT, PATO 3.1's own reference solution for ATC1 (`ref/FIAT_T.dat`, columns
$t, T(0), T(1), T(2), T(4), T(8), T(16), T(50\text{ mm})$): the largest $\lvert T_{\rm PATO}-T_{\rm FIAT}\rvert$
over the six probes and FIAT's 1 s grid $t=0\ldots60$ s, PATO interpolated linearly in time onto that
grid. It is computed by the check from the committed files each time it runs: 13.18 K, at (1 mm, 1 s).

## Results

| Quantity | Value |
|---|---|
| max over the six probes and $t=1\ldots60$ s (FIAT's grid) of $\lvert T_{\rm FUSS}-T_{\rm PATO}\rvert$ | 4.01 K, at (8 mm, 60 s) |
| the tolerance, max $\lvert T_{\rm PATO}-T_{\rm FIAT}\rvert$ over the same probes and grid | 13.18 K, at (1 mm, 1 s) |

Reported, not gated: the metric against FIAT directly, max $\lvert T_{\rm FUSS}-T_{\rm FIAT}\rvert$
over the six probes and the full window, is 16.03 K (at 2 mm, 4 s); the pyrolysis front depths
($\rho_s$ 98 %/2 % thresholds between $\rho_c=220$ and $\rho_v=280$ kg/m³, 15.21 mm and 2.91 mm below
the surface at the last field write) and the surface $\dot m_{pg}$ there, $5.387\times10^{-3}$ kg/(m²s).

## Pass criteria

max over the six probes and FIAT's 1 s grid $t=1\ldots60$ s of $\lvert T_{\rm FUSS}-T_{\rm PATO}\rvert\le$
the PATO-FIAT spread over the same probes and grid (computed from the committed reference files; no
stored FUSS value enters). FUSS's value at a thermocouple depth is the linear interpolation between the
two bracketing probes, then linear in time; PATO's is linear in time.
