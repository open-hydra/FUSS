<p align="center">
  <h1 align="center">FUSS</h1>
  <p align="center"><b>Fourier Un-Steady Solver</b></p>
</p>
<p align="center">
  <img src="https://zenodo.org/badge/1239585241.svg" alt="DOI: https://doi.org/10.5281/zenodo.20395702">
  <a href="https://open-hydra.github.io/FUSS/"><img src="https://img.shields.io/badge/docs-online-brightgreen.svg" alt="Documentation"></a>
  <img src="https://img.shields.io/badge/language-Fortran-734f96.svg" alt="Language: Fortran">
  <a href="https://github.com/open-hydra/FUSS/blob/main/LICENSE"><img src="https://img.shields.io/badge/license-GPLv3-blue.svg" alt="License: GPLv3"></a>
</p>

---

FUSS is an open-source solid heat-conduction solver for the energy equation on multi-block structured grids. Written in modern Fortran, it targets a wide range of steady-state and transient thermal problems in solid materials — from simple heated plates to multi-material assemblies with convective, radiative, and coupled boundary conditions.

## Features

- **Solid heat conduction** — steady-state and time-accurate computations on multi-block structured grids.
- **Multi-material support** — multiple solid regions with temperature-dependent material properties (density, specific heat, thermal conductivity) read from tabulated data.
- **Pyrolysing porous materials** — charring solids with Arrhenius kinetics; the pyrolysis gas by a quasi-steady one-dimensional outflow (type 1) or a Darcy pressure solve (type 2).
- **Explicit time integration** — forward Euler, RK2 and RK3, with optional Implicit Residual Smoothing (IRS).
- **Convergence acceleration** — local time stepping, grid sequencing (multigrid levels) and FAS V/W cycles for steady-state solutions.
- **Rich boundary conditions** — fixed temperature, prescribed heat flux, convective heat transfer with a prescribed flux, surface radiation with or without convection, symmetry, multi-block connection, periodic faces, chimera (overset) blocks, an ablating wall against B′ tables, and coupling to a fluid solver in hydra. Supports volumetric heat sources.
- **Parallel execution** — shared-memory parallelism via OpenMP; MPI support for distributed-memory runs.
- **Flexible I/O** — solution output in Tecplot (ASCII and binary) and VTK formats. Point probes for time-history recording. Restart capability.

## Quick Start

### Prerequisites

| Requirement | Details |
|---|---|
| **CMake** | ≥ 3.23 |
| **Fortran compiler** | GNU (`gfortran`) or Intel/oneAPI (`ifort` / `ifx`) |
| **C/C++ compiler** | Required for the ORION I/O library |

### Build

```bash
git clone --recurse-submodules https://github.com/open-hydra/FUSS.git
cd FUSS

# Build with GNU compilers and OpenMP
./install.sh build --compilers=gnu --use-openmp

# — or with Intel compilers and MPI —
./install.sh build --compilers=intel --use-openmp --use-mpi --use-tecio
```

The executable is placed in `bin/FUSS`.

See the [Installation Guide](https://open-hydra.github.io/FUSS/getting-started/installation/) for all build options, CMake presets, and troubleshooting.

### Run a Verification Case

```bash
cd test/steady-state/temperature_Plate
./FUSS.sh solve
```

See the [Quick Start](https://open-hydra.github.io/FUSS/getting-started/quick-start/) for a full walkthrough.

### Run the Test Suite

```bash
ctest --test-dir build -L fast                        # 25 tests (26 with TecIO), about a minute
ctest --test-dir build -LE "needs-mpi|verification"   # the pull-request check (.github/workflows/ctest.yml)
ctest --test-dir build -L verification -LE gate       # observed order and FAS = single grid, about a minute on 6 cores
ctest --test-dir build -L gate                        # the gate scripts (-L verification selects them too)
ctest --test-dir build                                # all: verification cases, smoke runs, MPI checks
```


## Dependencies

FUSS is built on top of companion libraries: FiNeR v2.3.0 is the Git submodule `lib/FiNeR`, ORION a release (2.1 or a later 2.x) that the configuration checks out in `lib/ORION`:

| Library | Role |
|---|---|
| [ORION](https://github.com/MarcoGrossi92/ORION) | Multi-format I/O (Tecplot, VTK, Plot3D) |
| [FiNeR](https://github.com/szaghi/FiNeR) | INI configuration file parser |
| [MORPH](https://github.com/open-hydra/MORPH) | Mesh motion: motion laws, cell geometry, swept volumes |

FiNeR's own libraries (PENF, StringiFor, FACE, FLAP, BeFoR64) are cloned by the configuration. ORION and FiNeR are checked against the releases FUSS requires (see [Versioning](https://open-hydra.github.io/FUSS/development/versioning/)); MORPH is not; run `git submodule update --init` after pulling a change of the FiNeR pin. Optional external libraries: **OpenMP**, **MPI**, **TecIO**.

## Project Structure

```
FUSS/
├── src/
│   ├── app/           # Main application
│   └── lib/           # Solver library sources
│       ├── base/      #   Core types and data structures
│       ├── config/    #   Configuration and INI parsing
│       ├── diagnostic/#   Diagnostic utilities
│       ├── driver/    #   Solver procedures (setup/solve/postprocess)
│       ├── io/        #   I/O routines (solution, probes, wall, BCs)
│       ├── numerics/  #   Numerical methods
│       │   ├── fluxes/    #   Diffusive fluxes and BC fluxes
│       │   ├── multigrid/ #   Grid sequencing, FAS cycles (restriction, prolongation)
│       │   ├── space/     #   Ghost cells and metrics
│       │   └── time/      #   Time step, explicit RK, IRS
│       ├── parallel/  #   MPI ghost-cell exchange
│       └── physics/   #   Material property evaluations
├── lib/               # FiNeR (submodule) and the ORION release
├── test/              # Test cases
│   ├── steady-state/                 # Steady-state V&V problems (NAFEMS benchmarks)
│   │   ├── temperature_Plate/        #   Fixed-temperature boundary
│   │   ├── convective_Plate/         #   Convective boundary condition
│   │   ├── radiative_Rod/            #   Radiation boundary condition
│   │   ├── convective_qrad_Rod/      #   Convective wall with a prescribed flux
│   │   └── qvol_Plate/               #   Volumetric heat source
│   ├── transient/                    # Time-accurate V&V problems
│   │   ├── oscillating_Rod/          #   1-D oscillating temperature (NAFEMS Test 5)
│   │   ├── multimat_Plate/           #   2-D multi-material transient plate
│   │   └── variable_cp_Plate/        #   Temperature-dependent heat capacity
│   ├── numerics-features/            # Numerical feature checks
│   │   ├── chimera/                  #   Chimera (overset) blocks
│   │   ├── periodic_sector/          #   Periodic faces (BC 201)
│   │   ├── residual_convergence/     #   Convergence strategy comparison
│   │   └── time_convergence/         #   Parallel scaling performance
│   ├── pyrolysis/                    # Pyrolysing solids, the Darcy gas closure, the ablating wall
│   ├── mpi/                          # MPI bit-identity and thread-funnel checks
│   ├── CMakeLists.txt                # CTest definitions
│   └── run_case.sh                   # Runs one case for CTest
├── cmake/             # CMake modules
├── install.sh         # Build helper script
└── CMakeLists.txt
```

## Documentation

Full documentation is available at **[open-hydra.github.io/FUSS](https://open-hydra.github.io/FUSS/)**, covering:

- Installation & quick start
- User guide & input file reference
- Verification & validation results
- Theory guide (governing equations, numerical methods, boundary conditions)

## License

FUSS is free and open-source software released under the [GNU General Public License v3.0](LICENSE).
