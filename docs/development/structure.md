# Code Structure

This page documents the repository layout, the internal architecture
of the FUSS library, and the solver execution pipeline.  All diagrams
use [Mermaid](https://mermaid.js.org/) and render directly in the
documentation.

---

## Repository Layout

```
FUSS/
├── CMakeLists.txt          # Top-level CMake build
├── CMakePresets.json        # Developer presets (compilers, paths)
├── install.sh               # Build / compile / update helper
├── mkdocs.yml               # Documentation site configuration
│
├── src/
│   ├── app/                 # Executables
│   │   ├── main.f90         # FUSS solver entry point
│   │   └── docgen.f90       # Input-parameter documentation generator
│   └── lib/                 # FUSS library (libFUSSL)
│       ├── base/            # Fundamental types and global parameters
│       ├── config/          # Input parsing, registry, setup assignment
│       ├── diagnostic/      # Residual monitoring
│       ├── driver/          # High-level solver orchestration
│       ├── io/              # File I/O (solution, BCs, probes, walls)
│       ├── numerics/        # Numerical methods
│       │   ├── fluxes/      #   Diffusive flux computation
│       │   │   └── bc/      #     Boundary-condition flux routines
│       │   ├── multigrid/   #   Restriction / prolongation
│       │   ├── space/       #   Metrics, ghost-cell filling
│       │   └── time/        #   Time stepping
│       │       └── explicit/#     Runge–Kutta, IRS, state update
│       ├── parallel/        # MPI ghost-cell exchange
│       └── physics/         # Solid material properties
│
├── lib/                     # External dependencies (git submodules)
│   ├── ORION/               # Structured-grid I/O library
│   └── third_party/
│       └── FiNeR/           # INI file parser
│
├── test/                    # Test suite
│   ├── steady-state/        # Steady-state V&V problems (NAFEMS benchmarks)
│   ├── transient/           # Time-accurate V&V problems
│   └── numerics-features/   # Numerical features verification
│       ├── residual_convergence/   # Convergence strategy comparison
│       └── time_convergence/       # Parallel scaling performance
│
├── docs/                    # MkDocs documentation source
├── cmake/                   # CMake modules (flags, OpenMP, MPI, etc.)
├── bin/                     # Built executables (FUSS, DocGen)
└── build/                   # Build artefacts
```

---

## Dependency Graph

The diagram below shows how the FUSS library depends on its
submodules and third-party components.

```mermaid
graph TD
    FUSS["<b>FUSS</b><br/>Solid heat-conduction solver"]
    ORION["<b>ORION</b><br/>Structured-grid I/O"]
    FiNeR["<b>FiNeR</b><br/>INI parser"]

    FUSS --> ORION
    FUSS --> FiNeR

    style FUSS fill:#37474f,stroke:#cfd8dc,color:#fff
    style ORION fill:#1565c0,stroke:#90caf9,color:#fff
    style FiNeR fill:#6a1b9a,stroke:#ce93d8,color:#fff
```

Optional compile-time dependencies (enabled via CMake flags):

```mermaid
graph LR
    FUSS["FUSS"]
    OMP["OpenMP"]
    MPI["MPI"]
    TecIO["TecIO"]

    FUSS -.->|USE_OPENMP| OMP
    FUSS -.->|USE_MPI| MPI
    FUSS -.->|USE_TECIO| TecIO

    style FUSS fill:#37474f,stroke:#cfd8dc,color:#fff
    style OMP fill:#00695c,stroke:#80cbc4,color:#fff
    style MPI fill:#00695c,stroke:#80cbc4,color:#fff
    style TecIO fill:#00695c,stroke:#80cbc4,color:#fff
```

---

## Library Architecture

The FUSS library (`libFUSSL`) is organised in six layers.  Lower
layers have no knowledge of higher layers.

```mermaid
graph TB
    subgraph driver ["<b>driver/</b> — Solver orchestration"]
        Procedures["Procedures_m<br/><i>FUSS_type</i>"]
        WSetup["Wrap_Setup"]
        WSolve["Wrap_Solve"]
        WPost["Wrap_Postprocess"]
        Alloc["Mod_Allocate_Data"]
    end

    subgraph config ["<b>config/</b> — Input & configuration"]
        ReadIni["Read_Ini"]
        Registry["Registry"]
        BackendINI["Backend_INI"]
        AssignSetup["Assign_Setup"]
        ConfigTypes["Config_Types_m"]
    end

    subgraph numerics ["<b>numerics/</b> — Numerical methods"]
        Fluxes["Mod_Fluxes<br/>Diffusive flux"]
        BCFluxes["Mod_BC_Fluxes"]
        Explicit["Mod_Explicit"]
        RK["Lib_RK"]
        IRS["Lib_IRS"]
        Newstate["Mod_Newstate"]
        MG["Mod_Multigrid"]
        Metrics["Mod_Metrics"]
        DT["Mod_dt"]
        Ghost["Lib_Ghost"]
    end

    subgraph physics ["<b>physics/</b> — Material properties"]
        Solid["Lib_Solid<br/><i>κ · ρ · cp · h</i>"]
    end

    subgraph io ["<b>io/</b> — File I/O"]
        IOSol["IO_Solution"]
        IOBC["IO_BC"]
        IOProbes["IO_Probes"]
        IOWall["IO_Wall"]
        LoadTable["Load_Table"]
    end

    subgraph base ["<b>base/</b> — Fundamental types"]
        Global["Global_m"]
        BaseTypes["Base_Types_m"]
        Advanced["Advanced_Types_m"]
        Params["Parameters_m"]
    end

    Procedures --> WSetup & WSolve & WPost
    WSetup --> config & io & base
    WSolve --> Explicit
    Explicit --> Fluxes & BCFluxes & Newstate & DT
    Fluxes --> Solid & Metrics
    Newstate --> RK & IRS
    WSolve --> MG
    config --> base
    io --> base

    style driver fill:#263238,stroke:#90a4ae,color:#eceff1
    style config fill:#1a237e,stroke:#7986cb,color:#e8eaf6
    style numerics fill:#004d40,stroke:#80cbc4,color:#e0f2f1
    style physics fill:#4a148c,stroke:#ce93d8,color:#f3e5f5
    style io fill:#bf360c,stroke:#ff8a65,color:#fbe9e7
    style base fill:#3e2723,stroke:#a1887f,color:#efebe9
```

---

## Module Hierarchy

Each source directory contains Fortran modules following a consistent
naming convention:

| Prefix | Role | Example |
|--------|------|---------|
| `Mod_*` | Module defining types, data, and procedure pointers | `Mod_Fluxes`, `Mod_Metrics` |
| `Lib_*` | Library of pure computational routines | `Lib_RK`, `Lib_Diffusive`, `Lib_Solid` |
| `*_m` | Fundamental type/parameter modules | `Global_m`, `Config_Types_m` |
| `Wrap_*` | High-level driver wrappers | `Wrap_Setup`, `Wrap_Solve` |
| `IO_*` | File I/O routines | `IO_Solution`, `IO_BC` |
| `Read_*` | Input file parsing | `Read_Ini` |
| `Load_*` | Database loading | `Load_Table` |

All public symbols are prefixed with `FUSS_` to avoid namespace
collisions when FUSS is linked as a library (e.g. inside HYDRA
coupling).

---

## Solver Pipeline

The main program (`src/app/main.f90`) creates a `FUSS_type` object
and calls three phases: **setup**, **solve** (in a loop), and
**postprocess**.

```mermaid
sequenceDiagram
    participant Main as main.f90
    participant Setup as Wrap_Setup
    participant Solve as Wrap_Solve
    participant Post as Wrap_Postprocess

    Main->>Setup: FUSS%setup(simulation)
    activate Setup
    Setup->>Setup: Load_Table (material properties)
    Setup->>Setup: Read_Inifile
    Setup->>Setup: Assign_Setup (function pointers)
    Setup->>Setup: Read_IC (initial conditions)
    Setup->>Setup: Setup_Data_Structure
    Setup->>Setup: Setup_Multigrid (if MGL > 1)
    Setup->>Setup: Setup_BC
    Setup->>Setup: Setup_Metrics (all grid levels)
    Setup->>Setup: Setup_Probes
    Setup->>Setup: Initialize_Wall_File
    Setup->>Setup: Fill_Ghost_Cell / Fill_matIDg
    deactivate Setup

    loop while TODO ≤ 2
        Main->>Solve: FUSS%solve(simulation)
        activate Solve
        Solve->>Solve: Explicit_Step (see below)
        Solve->>Solve: Restriction (if multigrid)
        deactivate Solve
        Main->>Post: FUSS%postprocess(simulation)
    end

    Main->>Post: FUSS%postprocess(simulation)
    Note right of Post: Final output + timing
```

### Explicit Step

Each call to `Explicit_Step` performs one complete time step with
Runge–Kutta sub-stages:

```mermaid
flowchart TB
    Start([Explicit_Step entry]) --> DT[Compute Δt<br/>VNN criterion]

    subgraph RK ["Runge–Kutta loop (n_RK stages)"]
        direction TB
        Ghost["Fill ghost cells<br/>(MPI exchange if parallel)"] --> Diff["Diffusive fluxes<br/>(Mod_Fluxes)"]
        Diff --> BCF["Boundary fluxes<br/>(Mod_BC_Fluxes)"]
        BCF --> Qvol["Volumetric source<br/>term (qvol)"]
        Qvol --> Update["RK_Newstate<br/>(IRS + state update)"]
    end

    DT --> RK
    RK --> Diag["Compute residual<br/>+ diagnostics"]
    Diag --> End([Return])

    style Start fill:#37474f,stroke:#cfd8dc,color:#fff
    style End fill:#37474f,stroke:#cfd8dc,color:#fff
    style RK fill:#004d40,stroke:#80cbc4,color:#e0f2f1
```

---

## Diffusive Flux Pipeline

The evaluation of the conductive heat flux at a single cell interface
follows this sequence:

```mermaid
flowchart LR
    A["Cell stencil<br/>(T at 10 points)"] --> B["Face temperature<br/>T = ½(T₁ + T₂)"]
    B --> C["Thermal conductivity<br/>κ = co_Kappa(matID, T)"]
    A --> D["Temperature gradient<br/>(normal + tangential)"]
    D --> E["Metric correction<br/>∇T = M · ∇T_comp"]
    C --> F["Heat flux<br/>q = κ A (∇T · n̂)"]
    E --> F

    style A fill:#263238,stroke:#90a4ae,color:#eceff1
    style F fill:#1b5e20,stroke:#a5d6a7,color:#fff
```

---

## Data Structures

### Simulation container

The top-level data type is `FUSS_simulation_type`, which holds an
array of grid levels (for multigrid) and an I/O container:

```mermaid
classDiagram
    class FUSS_simulation_type {
        +domain(:) : FUSS_domain_type
        +IOfield(:) : ORION_data
    }
    class FUSS_domain_type {
        +blk(:) : FUSS_block_type
        +nb : integer
        +iter : integer
        +time : real64
        +dtglobal : real64
    }
    class FUSS_block_type {
        +T(:,:,:) : real64
        +TO(:,:,:) : real64
        +R(:,:,:) : real64
        +dtlocal(:,:,:) : real64
        +qvol(:,:,:) : real64
        +matID(:,:,:) : real64
        +dim(3) : integer
    }

    FUSS_simulation_type "1" --> "*" FUSS_domain_type : domain
    FUSS_domain_type "1" --> "*" FUSS_block_type : blk
```

| Array | Shape | Content |
|-------|-------|---------|
| `T` | `(ni, nj, nk)` | Temperature at current time step |
| `TO` | `(ni, nj, nk)` | Temperature at previous time step |
| `R` | `(ni, nj, nk)` | Residual (flux accumulator) |
| `dtlocal` | `(ni, nj, nk)` | Local time step per cell |
| `qvol` | `(ni, nj, nk)` | Volumetric heat source term |
| `matID` | `(ni, nj, nk)` | Material ID per cell |

---

## Build System

### CMake targets

```mermaid
graph LR
    FUSSL["<b>FUSSL</b><br/>(Fortran library)"]
    FUSS_EXE["<b>FUSS</b><br/>(executable)"]
    DocGen["<b>DocGen</b><br/>(doc generator)"]
    ORION_LIB["ORION"]
    FiNeR_LIB["FiNeR"]

    FUSS_EXE --> FUSSL
    DocGen --> FUSSL
    FUSSL --> ORION_LIB
    FUSSL --> FiNeR_LIB

    style FUSSL fill:#37474f,stroke:#cfd8dc,color:#fff
    style FUSS_EXE fill:#1565c0,stroke:#90caf9,color:#fff
    style DocGen fill:#1565c0,stroke:#90caf9,color:#fff
    style ORION_LIB fill:#6a1b9a,stroke:#ce93d8,color:#fff
    style FiNeR_LIB fill:#6a1b9a,stroke:#ce93d8,color:#fff
```

| Target | Type | Description |
|--------|------|-------------|
| `FUSSL` | Static library | Core solver + all physics/numerics |
| `FUSS` | Executable | Standalone solver (`src/app/main.f90`) |
| `DocGen` | Executable | Input-parameter docs generator (`src/app/docgen.f90`) |

### Build workflow

```bash
# Option A: install.sh (recommended for first build)
./install.sh build --compilers=gnu

# Option B: CMake presets (for iterative development)
./install.sh compile          # uses CMakePresets.json
# or equivalently:
cmake --preset default && cmake --build build
```

Key CMake options:

| Option | Default | Effect |
|--------|:-------:|--------|
| `USE_OPENMP` | OFF | Enable OpenMP threading |
| `USE_MPI` | OFF | Enable MPI parallelism |
| `USE_TECIO` | OFF | Enable TecIO output format |

### Dependency paths

External library paths are set via CMake cache variables or the
`install.sh --include-*` flags:

| Variable | Default | Library |
|----------|---------|---------|
| `ORION_PATH` | `lib/ORION/` | ORION I/O |
| `FINER_PATH` | `lib/third_party/FiNeR/` | FiNeR INI parser |

---

## Naming Conventions

| Convention | Example | Meaning |
|------------|---------|---------|
| `FUSS_` prefix | `FUSS_Global_m` | Public Fortran module |
| `obj_` prefix | `obj_sim_param` | Global configuration singleton |
| `Lib_` prefix | `Lib_Solid` | Computational routine library |
| `Mod_` prefix | `Mod_Fluxes` | Module with types + pointers |
| `Wrap_` prefix | `Wrap_Solve` | Driver-level wrapper |
| `_m` suffix | `Config_Types_m` | Fundamental type module |
| `_type` suffix | `FUSS_block_type` | Derived type |
