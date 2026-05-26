# 2-D Steady State Thermal Analysis – Convective Plate

Steady-state heat conduction in a rectangular solid plate with a prescribed temperature at one edge, convection to ambient at two edges, and an insulated edge. This is a two-dimensional conduction problem with mixed boundary conditions whose reference solution is provided by the NAFEMS P16 thermal benchmark suite. Verification is performed by comparison of the computed wall temperature at point A (x = 0.6 m, y = 0.2 m) against the NAFEMS reference value of 18.25 °C.

**Reference**: NAFEMS Publication P16, "Benchmark Tests for Thermal Analysis", Test 9 (I) and 9 (ii), YR3087, Vol. 2, 1986.

---

## Problem setup

A rectangular plate (0.6 m × 1.0 m) is subjected to a prescribed temperature of 100 °C (373.15 K) along its bottom edge. The right and top edges lose heat by convection to an ambient temperature of 0 °C (273.15 K) with a uniform heat transfer coefficient of 750 W/m²°C. The left edge is thermally insulated.

The governing equation for two-dimensional steady-state conduction with no internal heat generation is:

$$\frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2} = 0$$

The convection boundary condition applied at the right (x = L) and top (y = H) faces is:

$$-k\,\frac{\partial T}{\partial n} = h\,(T - T_\infty)$$

where $h = 750$ W/m²°C and $T_\infty = 0$ °C.

**Boundary conditions**

| Face | Location | Condition | Value |
|---|---|---|---|
| Face 1 | x = 0 (left) | Adiabatic | q = 0 |
| Face 2 | x = 0.6 m (right) | Convection | h = 750 W/m²°C, T_ref = 0 °C |
| Face 3 | y = 0 (bottom) | Prescribed temperature | T = 100 °C (373.15 K) |
| Face 4 | y = 1.0 m (top) | Convection | h = 750 W/m²°C, T_ref = 0 °C |

**Material properties**

| Property | Symbol | Value | Unit |
|---|---|---|---|
| Thermal conductivity | k | 52.0 | W/mK |
| Density | ρ | 7850 | kg/m³ |
| Specific heat | cp | 460 | J/kgK |

## Numerical setup

| Parameter | Value |
|---|---|
| Time scheme | RK3 |
| VNN | 2.0 |
| Steady-state convergence | time-accurate = false |
| Integration variables | primitive |
| Implicit residual smoothing | enabled (β = 0.5) |
| Residual threshold | 1 × 10⁻⁸ |

## Grid structure

The mesh is a 2D structured grid (60 × 100 cells) spanning the physical domain (0.6 m × 1.0 m) with uniform spacing in both directions (Δx = 0.01 m, Δy = 0.01 m).

Boundary conditions (FUSS block-face notation):

- **Face 1** (x = 0): Adiabatic (`type = wall`, `q = 0`)
- **Face 2** (x = 0.6 m): Convection (`type = wall`, `hconv = 750`, `Tref = 273.15`, `qrad = 0`)
- **Face 3** (y = 0): Prescribed wall temperature (`type = wall`, `T = 373.15`)
- **Face 4** (y = 1.0 m): Convection (`type = wall`, `hconv = 750`, `Tref = 273.15`, `qrad = 0`)
- **Faces 5–6**: `null` (degenerate 2-D z-direction)

## Results and verification

The target point A is located at (x = 0.6 m, y = 0.2 m) on the convective face (B1F2). Since the mesh stores temperature at cell centres, the value at y = 0.2 m is estimated as the arithmetic mean of the two flanking cells (centres at y = 0.195 m and y = 0.205 m).

| Location | NAFEMS | FUSS (interpolated) | Error |
|---|---|---|---|
| A (x = 0.6 m, y = 0.2 m) | 18.25 °C | 18.26 °C | 0.004 % |

The interpolated wall temperature agrees with the NAFEMS reference to within 0.01 °C (0.004 % relative error), demonstrating excellent accuracy for this two-dimensional conduction problem with mixed boundary conditions.
