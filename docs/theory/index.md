# Theoretical Guide

This section describes the physical model and the numerical methods implemented in FUSS: the heat-conduction equation and its boundary conditions, the finite-volume spatial discretisation, the explicit time integration and its stability limit, and the tabulated material properties.

<div class="grid cards" markdown>

-   :material-math-integral:{ .lg .middle } __Governing Equations__

    ---

    Transient heat equation in integral form, enthalpy and integration variables, steady-state limit, boundary-condition types and sign convention

    [:octicons-arrow-right-24: Governing equations](governing-equations.md)

-   :material-grid:{ .lg .middle } __Spatial Discretisation__

    ---

    Cell-centred finite volumes, face gradient and face metric, block connections, wall boundary conditions and ghost cells, grid metrics

    [:octicons-arrow-right-24: Spatial discretisation](numerics.md)

-   :material-timer-outline:{ .lg .middle } __Time Integration__

    ---

    Explicit Euler, RK2 and SSP RK3, VNN time-step limit, steady-state and time-accurate modes, implicit residual smoothing, grid sequencing

    [:octicons-arrow-right-24: Time integration](time-integration.md)

-   :material-thermometer:{ .lg .middle } __Material Properties__

    ---

    Property table per material, linear interpolation, enthalpy column, admissible temperature range, materials per block

    [:octicons-arrow-right-24: Material properties](thermo.md)

</div>

---

## Overview

FUSS advances the transient heat equation in conservative form on structured multi-block grids using a cell-centred finite volume method. The numerical pipeline can be summarised as follows:

| Stage | Method | Page |
|-------|--------|------|
| **Governing equation** | Transient (or steady-state) heat diffusion in solids | [Governing Equations](governing-equations.md) |
| **Spatial discretisation** | Cell-centred FVM, central differences with metric terms; wall flux $q_\text{in} = G\,(T_w - T) + C$ | [Spatial Discretisation](numerics.md) |
| **Time marching** | Explicit Euler (default), RK2 or SSP RK3 with VNN-limited local or global time step; IRS, grid sequencing | [Time Integration](time-integration.md) |
| **Material model** | Tabulated, temperature-dependent $k$, $\rho$, $c_p$ and enthalpy per material | [Material Properties](thermo.md) |

The single unknown is the **temperature** $T$; each cell also carries a material index and a volumetric heat source. There are no velocity, pressure, or species fields.
