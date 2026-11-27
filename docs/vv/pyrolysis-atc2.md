# Ablation Test Case 2, Chemistry Off

Ablation Test Case 2 with the surface chemistry switched off reduces the surface boundary condition
to a closed-form Robin wall, so that the only source of disagreement between FUSS and an independent
code is the volumetric pyrolysis and conduction model. The oracle is a local PATO 3.1 run of the same
tutorial, edited to the same closed-form balance. Both the type-1 and the type-2 (Darcy) gas closures
are checked against it.

Setup and run: `test/pyrolysis/ATC2-chemistryOff/README.md` (type 1); the type-2 twin
`test/pyrolysis/darcy/ATC2-chemistryOff-type2/` is covered by `test/pyrolysis/darcy/README.md`.

---

## Problem

The mesh and TACOT tables of [ATC1](pyrolysis-atc1.md): a 0.05 m slab, 100 cells graded 0.1 toward the
heated face. $T(0) = 300$ K; face 2 is the closed-form Robin wall BC 303, `hconv = 300` W/(m²K),
`qrad = 0`, `Tref = 2800` K from $t=0$; 60 s. FUSS brackets each of PATO's seven probe depths (1, 2, 4,
8, 12, 16, 24 mm) with two index-position probes, interpolated linearly in $x$, then in time onto PATO's
0.1 s grid by the cubic through the four nearest probe rows. A probe row is written at the first time
step at or after each multiple of 0.1 s; with a linear interpolant in time, that offset alone moved the
type-2 maximum below by 0.05 K between time steps in the ratio 4 : 2 : 1, against 3 × 10⁻⁴ K of
the solution itself; with the cubic the three agree to 2 × 10⁻³ K.

The type-2 case solves the gas pressure ([FUSS-Darcy]) with TACOT's own permeability
($K_v=1.6\times10^{-11}$, $K_c=2\times10^{-11}$ m²); the heated 303 face is open to the gas, the back
face impermeable.

## Reference

A local PATO 3.1 run of `tutorials/1D/AblationTestCase_2.x_chemistryOff`, edited
(`edit_pato_case.py`) so that every row of `constant/porousMat/BoundaryConditions` reads
`hconv = 300, Tedge = 2800, chemistryOn = 0` (ramp and cool-down removed, `endTime = 60`) and the
emissivity column of the TACOT virgin/char tables is set to 0: both codes then evaluate the same
closed-form surface energy balance $h_{\rm conv}(T_{\rm edge}-T_w)$, with no re-radiation. PATO's
`output/porousMat/scalar/Ta_plot` (type 2 also `p_plot`) is read every 0.1 s at the seven depths.

A control run with the shipped (non-zero) emissivity gave a surface temperature of 1563 K at 60 s,
against 2332 K with the edit (2190 K at 1 mm), confirming the zero-emissivity edit took; PATO ran
12184 steps in 40 s wall time.

**Tolerance.** Only one external code solves this variant, so the temperature tolerance is borrowed from
[Ablation Test Case 2.1](ablating-wall.md#ablation-test-case-21-against-pato) (same material, same heating):
the largest $\lvert T_{\rm PATO}-T_{\rm Amaryllis}\rvert$ of that case over its seven depths, computed by the
check from the committed 2.1 reference files each time it runs: 25.35 K (1 mm, 0.3 s).

## Results

### Type 1

| Quantity | Value |
|---|---|
| max over the seven depths and $t=0.1$–60 s of $\lvert T_{\rm FUSS}-T_{\rm PATO}\rvert$ | 14.83 K, at 1 mm / 1.1 s |
| tolerance: the ATC 2.1 temperature spread | 25.35 K |

### Type 2 (Darcy)

| Check | Result |
|---|---|
| (T) max $\lvert T_{\rm FUSS}-T_{\rm PATO}\rvert$ over the seven depths (the open face evaluates the gas enthalpy leaving the domain at the wall temperature) | 10.57 K, at 1 mm / 0.9 s (tolerance: the ATC 2.1 temperature spread, 25.35 K) |
| (acc) same maximum vs. the type-1 run's own maximum, from the type-1 run made in the same check | $10.57\le14.83$ K (type 1 at 1 mm / 1.1 s): type 2 does not make the type-1 answer worse |
| (dir) $\lvert T_{\rm FUSS}-T_{\rm PATO}\rvert$ at the point of that type-1 maximum (1 mm, 1.1 s) | 10.41 K, below the type-1 value there (14.83 K) |
| (p) max $\lvert\Delta p_{\rm FUSS}-\Delta p_{\rm PATO}\rvert / \max\Delta p_{\rm PATO}$ | 0.5 % (11.2 Pa at 24 mm, 0.7 s; tolerance 19.66 %); FUSS peak $\Delta p = 2050.4$ Pa, PATO peak $2049.0$ Pa |
| (sign) $T_{\rm type2}-T_{\rm type1}$ wherever $T_{\rm PATO}\ge700$ K | from $-6.54$ K (4 mm) to $-3.12$ K (1 mm): type 2 runs cooler than type 1 at every hot probe, as the exact 0-D gas-storage shift predicts (table below) |

Type-2 minus type-1 shift at every depth, over the samples where $T_{\rm PATO}\ge700$ K (K, min .. max):

| Depth | 1 mm | 2 mm | 4 mm | 8 mm | 12 mm |
|---|---|---|---|---|---|
| $T_{\rm type2}-T_{\rm type1}$ | $-5.24$ .. $-3.12$ | $-5.90$ .. $-4.17$ | $-6.54$ .. $-4.32$ | $-6.36$ .. $-4.55$ | $-5.00$ .. $-4.60$ |

The largest (least negative) value over every hot sample is $-3.115$ K, at 1 mm: still below the
$0.05$ K slack.

## Pass criteria

- Type 1: $\max\lvert T_{\rm FUSS}-T_{\rm PATO}\rvert$ over the seven depths and $t=0.1$–60 s $\le$ the
  ATC 2.1 temperature spread (computed from that case's reference files; no stored FUSS value enters).
- Type 2, all of: (T) the same maximum $\le$ the ATC 2.1 temperature spread; (acc) the same maximum $\le$
  the type-1 run's own maximum, from the type-1 run made in the same check (removing the type-1
  approximation must not make the answer worse); (dir) at that type-1 maximum's point the type-2
  deviation is strictly smaller; (p) the pressure-rise error $\le$ the sum of the relative temperature,
  viscosity and molar-mass sensitivities of $\Delta p$ at the ATC 2.1 temperature spread (0.42 %/K at
  850 K from the TACOT data, 10.65 %) and ATC 2.1's $\dot m_g$ spread over PATO's $\dot m_g$ peak (9.02 %):
  19.66 % of PATO's peak $\Delta p$; (sign) $T_{\rm type2}-T_{\rm type1}\le0.05$ K
  (the slack of the $-\epsilon\,\mathrm dp/(\rho c)$ gas pressure-work term) wherever $T_{\rm PATO}\ge700$ K,
  i.e. type 2 is never hotter than type 1 at a hot probe beyond that round-off-sized slack.
