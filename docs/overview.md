---
title: Overview
---

# Overview

FUSS (**F**ourier **U**n-**S**teady **S**olver) is an open-source finite-volume solver for steady and time-accurate heat conduction in solids on multi-block structured grids, written in Fortran. It is built as a library (`libFUSSL`) and a stand-alone executable (`bin/FUSS`).

---

## Hydra CFD Suite

FUSS is the solid heat-conduction solver of [**hydra**](https://github.com/open-hydra/hydra), a multi-physics framework that builds its solvers as libraries and couples them. Inside hydra, FUSS exchanges wall data with the fluid solver through boundary conditions 103 and 104.

| Component | Role |
|-----------|------|
| [**ATLAS**](https://github.com/open-hydra/ATLAS) | Pre-processor: writes the files in `INPUT/` (property tables, initial conditions, boundary conditions, sources) |
| [**MOSE**](https://github.com/open-hydra/MOSE) | Compressible Euler / Navier–Stokes solver; coupled to FUSS inside hydra |
| **FUSS** | Heat conduction in solids (this package) |

!!! info "Using FUSS without ATLAS"
    FUSS reads the `[FUSS-*]` sections of `input.ini` and the files in `INPUT/`. The other sections of the bundled `input.ini` files (`[GRIB-*]`, `[GPB-*]`, `[ICB-*]`, `[BCB-*]`, `[STB-*]`) are read by the pre-processing tools, not by FUSS. Every test case ships its `INPUT/` files, and they can be written by hand or by script in the formats of the [User Guide](user/using.md).

---

## Governing equation

$$
\rho(T)\, c_p(T)\, \frac{\partial T}{\partial t}
= \nabla\!\cdot\!\bigl(k(T)\,\nabla T\bigr) + \dot{q}_\text{vol}
$$

$T$ is the temperature, $\rho$ the density, $c_p$ the specific heat, $k$ the conductivity and $\dot{q}_\text{vol}$ a volumetric source (W/m³). With `integration-variables = cons` (default) FUSS advances the volumetric enthalpy $h(T) = \int \rho c_p\, dT$ and recovers $T$ from it; with `prim` it advances $T$. See [Governing Equations](theory/governing-equations.md).

---

## Physical models

### Materials

Each material is a table in `INPUT/properties.dat` with the columns Temperature, Cp, Density, Conductivity and Energy ($h$), in 1 K steps. Properties are linearly interpolated between the two rows that bracket $T$; $T(h)$ is found by searching the Energy column. Materials are assigned per cell through the `matID` field of the initial condition. See [Material Properties](theory/thermo.md).

### Boundary conditions

Every boundary cell face has a type code in `INPUT/bc.txt`. ATLAS writes the codes from the face keywords and the wall sections of `input.ini`.

| Code | Behaviour | `input.ini` (ATLAS) |
|------|-----------|---------------------|
| 0 | Zero flux (mirror ghost) | `face* = null` |
| 300 | Symmetry: zero flux (mirror ghost) | `face* = symmetry` |
| 101 | Block connection | `face* = connection` |
| 102 | Chimera (overset) interface | `face* = chimera` |
| 103, 104 | Coupling with the hydra fluid solver | hydra only |
| 301 | Heat flux into the solid $q_w$, constant or time table | `q` or `q-time-file` |
| 302 | Wall temperature $T_w$, constant or time table | `T` or `T-time-file` |
| 303 | $q = h_\text{conv}(T_\text{ref} - T_w) + q_\text{rad}$ | `hconv`, `qrad`, `Tref` |
| 304 | $q = \varepsilon\sigma(T_\text{ref}^4 - T_w^4)$ | `eps`, `Tref` |
| 305 | $q = h_\text{conv}(T_\text{ref} - T_w) + \varepsilon\sigma(T_\text{ref}^4 - T_w^4)$ | `hconv`, `eps`, `Tref` |

The walls 301–305 share one routine: $T_w$ and the heat flux into the solid follow from the wall model and the half-cell conduction balance. See [Boundary Conditions](user/boundary-conditions.md).

### Volumetric heat source

$\dot{q}_\text{vol}$ is a cell field read from `INPUT/st.tec` when that file exists (zero otherwise).

---

## Numerical methods

### Space

| | |
|-|-|
| Grid | Structured multi-block hexahedral grids (1-D and 2-D cases have one cell in the unused directions), two ghost layers per block |
| Diffusive flux | Cell-centred finite volume; face gradient from the difference across the face plus tangential central differences |
| Face metric | Inverse of the mean of the two cells' inverse metrics: exact for linear fields on stretched orthogonal meshes |
| Face conductivity | $k$ at the mean of the two cell temperatures; at block connections the interface temperature follows from heat-flux continuity, each side with its own conductivity |
| Walls | Ghost cells $2T_w - T$; boundary flux from the wall model |

See [Spatial Discretisation](theory/numerics.md).

### Time

| | |
|-|-|
| Schemes | `euler` (default), `RK2` (Heun), `RK3` (SSP) |
| Time step | Diffusive limit set by `vnn`; local step for steady runs, global minimum for `time-accurate = true` |
| Implicit residual smoothing | `irs-beta`: two Jacobi sweeps per direction |
| Multigrid | Grid sequencing: coarsest level first, prolongation to the next finer level, no cycles |

See [Time Integration](theory/time-integration.md).

### Parallel computing

| Mode | Details |
|------|---------|
| OpenMP | Cell loops of each block shared among threads |
| MPI | Whole blocks per rank (largest first to the least-loaded rank); exchange of connection ghosts and chimera donors; output gathered on rank 0 |
| Hybrid | MPI ranks with OpenMP threads |

One process and N ranks write identical files. See [Parallel Execution](user/parallel.md).

---

## Dependencies

| Library | Role | Source |
|---------|------|--------|
| [ORION](https://github.com/MarcoGrossi92/ORION) | Tecplot and VTK reading/writing; bundles TecIO | Submodule `lib/ORION` |
| [FiNeR](https://github.com/szaghi/FiNeR) | `input.ini` parser | Submodule `lib/third_party/FiNeR`; its dependencies (PENF, StringiFor, FACE, FLAP, BeFoR64) are cloned by `install.sh` |
| OpenMP | Threads | Optional, `--use-openmp` |
| MPI | Distributed memory | Optional, `--use-mpi` |
| TecIO | Tecplot binary (`.szplt`) | Optional, `--use-tecio`; built by ORION |

Build tools: CMake 3.23 or newer, a Fortran compiler (GNU or Intel), and C and C++ compilers. The `verify.py` scripts of the test cases use Python 3 with NumPy and Matplotlib. See [Installation](getting-started/installation.md).

---

## Documentation guide

| Section | Content |
|---------|---------|
| [Getting Started](getting-started/index.md) | Build and first run |
| [User Guide](user/index.md) | Running cases, `input.ini`, initial and boundary conditions, parallel runs, output |
| [Input Parameters](user/registry.md) | Every `input.ini` key with default and allowed values |
| [Theoretical Guide](theory/index.md) | Equations, discretisation, time integration, material tables |
| [Verification & Validation](vv/index.md) | Benchmark and exact-solution cases, convergence and parallel studies |
| [Developer Guide](development/index.md) | Code structure, tests, contributing |
| [About](about/index.md) | License and acknowledgements |

---

## License

FUSS is released under the [GNU General Public License v3.0](about/license.md).
