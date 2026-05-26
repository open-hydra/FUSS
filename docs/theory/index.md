# Theoretical Guide

This section provides the theoretical foundations for the physical models and numerical methods implemented in FUSS. The material covers the governing equation of heat conduction in solids, the finite-volume spatial discretisation, the explicit time-integration scheme, stability conditions, and the thermophysical property framework.

<div class="grid cards" markdown>

-   :material-math-integral:{ .lg .middle } __Governing Equations__

    ---

    Transient heat equation in integral form, steady-state limit, volumetric heat generation, and boundary condition types

    [:octicons-arrow-right-24: Governing equations](governing-equations.md)

-   :material-grid:{ .lg .middle } __Spatial Discretisation__

    ---

    Cell-centred finite volume method, diffusive flux evaluation, harmonic-mean conductivity at material interfaces, multi-block connectivity

    [:octicons-arrow-right-24: Spatial discretisation](numerics.md)

-   :material-timer-outline:{ .lg .middle } __Time Integration__

    ---

    SSP Runge–Kutta, VNN diffusive stability condition, steady-state and time-accurate modes, implicit residual smoothing

    [:octicons-arrow-right-24: Time integration](time-integration.md)

-   :material-thermometer:{ .lg .middle } __Material Properties__

    ---

    Solid phase definition, constant thermophysical properties, multi-material block assignment, temperature cap

    [:octicons-arrow-right-24: Material properties](thermo.md)

</div>

---

## Overview

FUSS advances the transient heat equation in conservative form on structured multi-block grids using a cell-centred finite volume method. The numerical pipeline can be summarised as follows:

| Stage | Method | Page |
|-------|--------|------|
| **Governing equation** | Transient (or steady-state) heat diffusion in solids | [Governing Equations](governing-equations.md) |
| **Spatial discretisation** | Cell-centred FVM with second-order central differences | [Spatial Discretisation](numerics.md) |
| **Time marching** | 3-stage SSP RK3 with VNN-limited time step; IRS, multigrid | [Time Integration](time-integration.md) |
| **Material model** | Temperature-dependent $k$, $\rho$, $c_p$ per solid phase | [Material Properties](thermo.md) |

The single prognostic variable is **temperature** $T$. There are no velocity, pressure, or species fields. The solver is therefore substantially simpler than a general CFD code, while retaining the full multi-block structured grid infrastructure and the same explicit RK3 time-stepping kernel.

---
