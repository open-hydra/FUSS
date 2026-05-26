---
title: Overview
---

# Overview

FUSS (**F**ourier **U**n-**S**teady **S**olver) is an open-source solver for steady-state and time-accurate heat conduction in solid bodies on multi-block structured grids, written in modern Fortran. It targets a wide range of thermal problems — from one-dimensional conduction benchmarks to two- and three-dimensional solids with mixed boundary conditions, volumetric heat generation, and multi-material domains.

---

## Hydra CFD Suite

FUSS is the **solid heat-conduction solver** of the **Hydra** ecosystem — an integrated suite of tools for multi-physics simulation of complex systems.

| Component | Role | Status |
|-----------|------|--------|
| [**ATLAS**](https://github.com/open-hydra/ATLAS) | Pre-processor: mesh prep, initial & boundary conditions, material property data | Separate package |
| [**MOSE**](https://github.com/open-hydra/MOSE) | Solver: compressible Euler/Navier–Stokes with finite-rate chemistry | Sister solver |
| **FUSS** | Solver: transient and steady-state heat conduction in solids | This package |

!!! info "Using FUSS without ATLAS"
    The input files required by FUSS (initial conditions, boundary condition table, material property tables) are **typically produced by ATLAS**. If ATLAS is not available, all input files can be prepared manually; see the [User Guide](user/using.md) for the expected formats.

---

## FUSS Capabilities

FUSS solves the unsteady heat conduction equation for solid continua in conservative form:

$$
\rho(T)\, c_p(T)\, \frac{\partial T}{\partial t}
= \nabla\!\cdot\!\bigl(k(T)\,\nabla T\bigr) + \dot{q}_\text{vol}
$$

where $T$ is the temperature field, $\rho$ the density, $c_p$ the specific heat capacity, $k$ the thermal conductivity, and $\dot{q}_\text{vol}$ a volumetric heat source. Both the transient form and the steady-state limit ($\partial T/\partial t = 0$) are supported.

---

## Physical Models

### Solid material

FUSS models each computational block as a solid with three thermophysical properties — thermal conductivity, density, and specific heat capacity — all of which can depend on temperature. Properties are evaluated by linear interpolation of a precomputed table indexed by integer Kelvin (see [Material Properties](theory/thermo.md)).

| Property | Model |
|----------|-------|
| Thermal conductivity $k(T)$ | Tabulated, linear interpolation |
| Density $\rho(T)$ | Tabulated, linear interpolation |
| Specific heat $c_p(T)$ | Tabulated, linear interpolation |
| Volumetric enthalpy $h(T)$ | Cumulative trapezoid of $\rho c_p$ |
| Multi-material interfaces | Harmonic-mean conductivity (series-resistance correct) |

### Boundary Conditions

FUSS supports a catalogue of boundary behaviours selected via the parameter combination in named `[<name>]` sections (all using `type = wall`). The pre-processor translates each parameter set into an internal numeric code consumed by the solver. Mixed boundary conditions can be applied independently to each block face.

| Behaviour | Selected by | Description |
|-----------|-------------|-------------|
| Block connection | `face* = connection` | Block-to-block interface (continuity of $T$ and $k\nabla T\!\cdot\!\hat{\mathbf n}$) |
| Inactive face | `face* = null` | Degenerate face (e.g. faces 5/6 on 2-D plane cases) |
| Prescribed flux | `q` (or `q-time`) | Specified normal heat flux $q_w$; `q = 0` gives an adiabatic wall |
| Prescribed temperature | `T` (or `T-time`) | Dirichlet BC; constant or time-varying table |
| Convection (+ optional $q_\text{rad}$) | `hconv`, `Tref`, `qrad` | Newton's law of cooling: $q = h(T_\text{ref} - T_\text{wall}) + q_\text{rad}$ |
| Radiation | `eps`, `Tref` | Stefan–Boltzmann: $q = \varepsilon\sigma(T_\text{ref}^4 - T_\text{wall}^4)$ |
| Convection + radiation | `hconv`, `eps`, `Tref` | Sum of convective and radiative fluxes |

### Volumetric Heat Generation

Each block can include a uniform internal heat source $\dot{q}_\text{vol}$ (W/m³) declared in its initial-condition section. The source term is added directly to the cell residual.

---

## Numerical Methods

### Spatial discretisation

| | Details |
|-|---------|
| Framework | Cell-centred finite volume on structured multi-block hexahedral grids |
| Diffusive flux | Second-order central differences |
| Face conductivity | Arithmetic mean within a material; harmonic mean across material interfaces |
| Boundary fluxes | Ghost-cell formulation for Dirichlet, Neumann, convective, and radiative BCs |
| Metric tensor | Computational-to-physical mapping evaluated per face |

### Time integration

| | Details |
|-|---------|
| Explicit scheme | 3-stage strong-stability-preserving Runge–Kutta (SSP RK3) |
| Stability | Von Neumann number (VNN) condition on the diffusion operator |
| Steady-state mode | Local time stepping (pseudo-time marching) |
| Time-accurate mode | Global minimum $\Delta t$, monotone time advancement |
| Acceleration | Implicit Residual Smoothing (IRS); geometric multigrid (2 levels) |

### Parallel computing

| Mode | Details |
|------|---------|
| Shared memory | OpenMP thread-level parallelism within a block |
| Distributed memory | MPI domain decomposition across blocks (ghost-cell exchange) |
| Hybrid | OpenMP + MPI combined runs on HPC clusters |

---

## Code Dependencies

### Required libraries

| Library | Role | Source |
|---------|------|--------|
| [ORION](https://github.com/MarcoGrossi92/ORION) | Multi-format I/O — Tecplot, VTK | Bundled submodule |
| [FiNeR](https://github.com/szaghi/FiNeR) | INI configuration file parser | Bundled submodule |

### Optional libraries

| Library | Role |
|---------|------|
| OpenMP | Shared-memory thread parallelism |
| MPI | Distributed-memory parallelism |
| TecIO | Binary Tecplot output |

### Build toolchain

| Tool | Minimum version |
|------|----------------|
| CMake | 3.23 |
| Fortran compiler | GNU gfortran 11+ or Intel ifx/ifort |
| C / C++ compiler | GCC or ICC (for ORION and the optional TecIO) |
| Python | 3.9+ (for verification scripts) |

---

## Documentation Guide

| Section | What you'll find |
|---------|-----------------|
| [**Getting Started**](getting-started/index.md) | Installation, prerequisites, and first run |
| [**User Guide**](user/index.md) | Running simulations, configuring input files, boundary conditions, output |
| [**Input File Reference**](user/input.md) | `input.ini` structure, all sections, auto-generated parameter registry |
| [**Theory Guide**](theory/index.md) | Governing equations, numerical methods, material property model |
| [**Verification & Validation**](vv/index.md) | NAFEMS benchmarks and analytical test cases |
| [**Developer Guide**](development/index.md) | Repository architecture, testing framework, contribution guidelines |
| [**About**](about/index.md) | License, acknowledgements, and contributors |

---

## License

FUSS is free and open-source software released under the **[GNU General Public License v3.0](about/license.md)** (GPL-3.0).

| Permission | |
|------------|-|
| :white_check_mark: Use freely | For any purpose, including commercial |
| :white_check_mark: Modify | Change the source code as needed |
| :white_check_mark: Distribute | Share original or modified versions |
| :white_check_mark: Patent grant | Contributors grant patent rights |
| :warning: Share-alike | Derivative works must use GPL-3.0 |
| :warning: Disclose source | Source code must be provided when distributing |

Full license text: [`LICENSE`](about/license.md)
