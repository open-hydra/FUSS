<p align="center">
  <h1 align="center">FUSS</h1>
  <p align="center"><b>Fourier Un-Steady Solver</b></p>
</p>

<p align="center">
  <a href="https://github.com/open-hydra/FUSS/blob/main/LICENSE"><img src="https://img.shields.io/badge/license-GPLv3-blue.svg" alt="License: GPLv3"></a>
  <a href="https://open-hydra.github.io/FUSS/"><img src="https://img.shields.io/badge/docs-online-brightgreen.svg" alt="Documentation"></a>
  <img src="https://img.shields.io/badge/language-Fortran-734f96.svg" alt="Language: Fortran">
</p>

---

FUSS is an open-source solid heat-conduction solver for the energy equation on multi-block structured grids. Written in modern Fortran, it targets a wide range of steady-state and transient thermal problems in solid materials — from simple heated plates to multi-material assemblies with convective, radiative, and coupled boundary conditions.

## Features

- **Solid heat conduction** — steady-state and time-accurate computations on multi-block structured grids.
- **Multi-material support** — multiple solid regions with temperature-dependent material properties (density, specific heat, thermal conductivity) read from tabulated data.
- **Explicit time integration** — Runge–Kutta schemes with optional Implicit Residual Smoothing (IRS) for stability enhancement.
- **Convergence acceleration** — multigrid cycling and local time stepping for efficient steady-state solutions.
- **Rich boundary conditions** — fixed temperature, prescribed heat flux, convective heat transfer, surface radiation, symmetry, and multi-block connection. Supports volumetric heat sources.
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


## Dependencies

FUSS is built on top of companion libraries, included as Git submodules:

| Library | Role |
|---|---|
| [ORION](https://github.com/MarcoGrossi92/ORION) | Multi-format I/O (Tecplot, VTK, Plot3D) |
| [FiNeR](https://github.com/szaghi/FiNeR) | INI configuration file parser |

Optional external libraries: **OpenMP**, **MPI**, **TecIO**.

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
│       │   ├── fluxes/#     Diffusive flux computation and BCs
│       │   ├── space/ #     Spatial discretization and metrics
│       │   └── time/  #     Time integration (explicit RK, IRS, multigrid)
│       ├── parallel/  #   MPI ghost-cell exchange
│       └── physics/   #   Material property evaluations
├── lib/               # Git submodule dependencies
├── test/              # Test cases
│   ├── steady-state/                 # Steady-state V&V problems (NAFEMS benchmarks)
│   │   ├── temperature_Plate/        #   Fixed-temperature boundary
│   │   ├── convective_Plate/         #   Convective boundary condition
│   │   ├── radiative_Rod/            #   Radiation boundary condition
│   │   └── qvol_Plate/               #   Volumetric heat source
│   ├── transient/                    # Time-accurate V&V problems
│   │   ├── oscillating_Rod/          #   1-D oscillating temperature (NAFEMS Test 5)
│   │   └── multimat_Plate/           #   2-D multi-material transient plate
│   └── numerics-features/            # Numerical feature verification
│       ├── residual_convergence/     #   Convergence strategy comparison
│       └── time_convergence/         #   Parallel scaling performance
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
