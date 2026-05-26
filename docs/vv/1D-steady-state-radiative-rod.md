# 1-D Steady State Thermal Analysis – Radiative Rod

Steady-state heat conduction in a rectangular solid rod with a prescribed temperature at one end and a radiation boundary condition at the other. This configuration is a standard one-dimensional conduction–radiation problem whose solution requires solving a nonlinear algebraic equation for the free-end temperature. Verification is performed by comparison with the NAFEMS P16 benchmark reference value.

**Reference**: NAFEMS Publication P16, "Benchmark Tests for Thermal Analysis", Test 9 (iii), YR3087, Vol. 2, 1986.

---

## Problem setup

A rectangular plate (0.1 m × 0.01 m) is subject to a prescribed temperature of 1000 K at its left-hand face and radiation to an ambient temperature of 300 K at its right-hand face. The problem reduces to a one-dimensional steady-state balance between conduction along the rod and radiative heat loss at the free end.

The governing equation for steady-state conduction with no internal heat generation is:

$$\frac{d}{dx}\left(k\frac{dT}{dx}\right) = 0$$

which yields a linear temperature profile along the rod. The radiation boundary condition at $x = L$ couples the conductive heat flux to the Stefan–Boltzmann radiative flux:

$$-k\left.\frac{dT}{dx}\right|_{x=L} = \varepsilon\,\sigma\left(T_L^4 - T_\text{ref}^4\right)$$

This nonlinear equation for the unknown temperature $T_L$ at the radiating face must be solved iteratively. Substituting the linear profile $dT/dx = (T_L - T_1)/L$ gives:

$$\frac{k\,(T_1 - T_L)}{L} = \varepsilon\,\sigma\left(T_L^4 - T_\text{ref}^4\right)$$

**Boundary conditions**

| Face | Condition | Value |
|---|---|---|
| Face 1 (x = 0) | Prescribed temperature | T = 1000 K |
| Face 2 (x = L) | Radiation to ambient | ε = 0.98, T_ref = 300 K |
| Faces 3–6 | Adiabatic (symmetry) | — |

**Material properties**

| Property | Symbol | Value | Unit |
|---|---|---|---|
| Thermal conductivity | k | 55.563 | W/mK |
| Density | ρ | 7850 | kg/m³ |
| Specific heat | cp | 460 | J/kgK |
| Stefan–Boltzmann constant | σ | 5.67 × 10⁻⁸ | W/m²K⁴ |
| Emissivity | ε | 0.98 | — |

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

The mesh is a 2D structured grid (100 × 2 cells) spanning the physical domain (0.1 m × 0.01 m) with uniform spacing in both directions. The two-cell depth in the transverse direction makes the problem effectively one-dimensional, consistent with the NAFEMS test specification.

Boundary conditions (FUSS block-face notation):

- **Face 1** (x = 0): Prescribed wall temperature (`type = wall`, `T = 1000`)
- **Face 2** (x = L): Radiation (`type = wall`, `eps = 0.98`, `Tref = 300`)
- **Faces 3–6**: `null` (degenerate 2-D y- and z-directions; effectively adiabatic)

## Results and verification

The FUSS solution is compared against the NAFEMS benchmark target temperature at the radiation face (x = 0.1 m).

| Location | NAFEMS | FUSS | Error |
|---|---|---|---|
| x = 0.1 m (Face 2) | 927 K | 926.356 K | 0.069 % |

The computed wall temperature at the radiation boundary agrees with the NAFEMS reference to within 0.07 %, well below the 1 % acceptance criterion.
