# Code Structure

This page describes the repository layout, the modules of the FUSS library and the order in which the solver calls them.

---

## Repository layout

```
FUSS/
├── CMakeLists.txt        # top-level build: ORION, FiNeR, src/lib, src/app, test
├── install.sh            # build / compile / update
├── mkdocs.yml            # documentation site
├── cmake/                # compiler flags, OpenMP/MPI detection
├── src/
│   ├── app/
│   │   ├── main.f90      # solver executable (bin/FUSS)
│   │   └── docgen.f90    # writes docs/user/registry.md (bin/DocGen)
│   └── lib/              # FUSS library (FUSSL), see below
├── lib/
│   ├── ORION/            # submodule: Tecplot/VTK I/O, bundles TecIO
│   └── third_party/
│       ├── FiNeR/        # submodule: INI parser
│       └── PENF/ FACE/ FLAP/ BeFoR64/ StringiFor/   # submodules: FiNeR's libraries, added by CMake before FiNeR
├── test/
│   ├── CMakeLists.txt    # CTest definitions
│   ├── run_case.sh       # runs one case for CTest
│   ├── steady-state/     # verification cases
│   ├── transient/        # verification cases
│   └── numerics-features/# chimera, residual convergence, parallel timing
└── docs/                 # MkDocs sources
```

`bin/`, `build/` and `CMakePresets.json` are created by the build and ignored by git.

---

## Library sources

| Directory | Files | Content |
|-----------|-------|---------|
| `base/` | `Parameters_m`, `Global_m`, `Base_Types_m`, `Advanced_Types_m`, `Series_Data_m`, `Finite_m` | Constants, ghost layers `gc` and dimensions `ndir`, block/domain/BC types, time tables of BCs 301/302, test of a value for NaN and infinities (`Not_Finite`) |
| `config/` | `Config_Types_m`, `Registry`, `Backend_INI`, `Read_Ini`, `Register_Sim_Param`, `Register_IO`, `Register_Numerics`, `Assign_Setup` | Settings objects (`obj_*`), `input.ini` key registry and reader, scheme selection |
| `driver/` | `Procedures_m`, `Wrap_Setup`, `Wrap_Solve`, `Wrap_Postprocess`, `Mod_Allocate_Data` | `FUSS_type` with `setup`, `solve`, `postprocess`; allocation |
| `io/` | `Load_Table`, `IO_Solution`, `IO_BC`, `IO_Wall`, `IO_Probes` | Property tables, initial condition and solution files, `bc.txt`, wall output, probes |
| `numerics/fluxes/` | `Mod_Fluxes`, `Lib_Diffusive` | Interior diffusive fluxes |
| `numerics/fluxes/bc/` | `Mod_BC_Fluxes`, `Lib_BC_Fluxes`, `Lib_BC_Fluxes_Connection`, `Lib_BC_Fluxes_Wall` | Boundary fluxes: dispatch, face helpers, block connection and chimera (101, 102), walls (301–305) |
| `numerics/space/` | `Mod_Metrics`, `Lib_Metrics`, `Lib_Stencil_Metrics`, `Lib_Ghost` | Metrics, face normals and areas; stencil face metrics (`diffusive-metric = stencil`); ghost-cell fill for every BC type |
| `numerics/time/` | `Mod_dt` | Local and global time step |
| `numerics/time/explicit/` | `Mod_Explicit`, `Mod_Newstate`, `Lib_Newstate`, `Lib_RK`, `Lib_IRS` | Iteration loop, state update (`cons`/`prim`), Runge–Kutta stages, residual smoothing |
| `numerics/multigrid/` | `Mod_Multigrid`, `Lib_Multigrid`, `Mod_FAS` | Coarse levels, restriction, prolongation; FAS V and W cycles |
| `parallel/` | `Mod_MPI`, `Mod_GhostExchange`, `Mod_Timers` | MPI set-up, block partition, reductions, the stops of the errors (`mpi_stop_all`, `mpi_stop_if_any`, `stop_from_thread`); ghost and chimera exchange, gathers to rank 0; wall-clock timers |
| `physics/` | `Lib_Solid` | Property interpolation `co_Kappa`, `co_Rho`, `co_Cp`, `co_H` and inversion `co_T` |
| `diagnostic/` | `Mod_Diagnostic` | Residual norm, `diagnostic.tec` |

All files are `.f90`. `src/lib/CMakeLists.txt` collects every `.f90` under `src/lib/` with `GLOB_RECURSE`.

---

## Dependencies

```mermaid
graph TD
    FUSS["<b>FUSSL</b><br/>FUSS library"]
    ORION["<b>ORION</b><br/>Tecplot / VTK I/O"]
    FiNeR["<b>FiNeR</b><br/>INI parser"]
    FDeps["PENF · StringiFor · FACE · FLAP · BeFoR64"]
    TecIO["TecIO (optional, USE_TECIO)"]
    OMP["OpenMP (optional, USE_OPENMP)"]
    MPI["MPI (optional, USE_MPI)"]

    FUSS --> ORION
    FUSS --> FiNeR
    FiNeR --> FDeps
    ORION -.-> TecIO
    FUSS -.-> OMP
    FUSS -.-> MPI
```

---

## Module dependencies

The main `use` relations between the modules (every module also uses `base/` and most use `Config_Types_m`; module names carry the prefix `FUSS_`, omitted here):

```mermaid
graph TB
    subgraph driver ["driver/"]
        Procedures["Procedures_m<br/><i>FUSS_type</i>"]
        WSetup["Wrap_Setup"]
        WSolve["Wrap_Solve"]
        WPost["Wrap_Postprocess"]
        Alloc["Mod_Allocate_Data"]
    end

    subgraph config ["config/"]
        ReadIni["Read_Ini"]
        Assign["Assign_Setup"]
    end

    subgraph io ["io/"]
        LoadTable["Load_Table"]
        IOSol["IO_Solution"]
        IOBC["IO_BC"]
        IOWall["IO_Wall"]
        IOProbes["IO_Probes"]
    end

    subgraph numerics ["numerics/"]
        Explicit["Mod_Explicit"]
        DT["Mod_dt"]
        Ghost["Lib_Ghost"]
        Fluxes["Mod_Fluxes"]
        Diff["Lib_Diffusive"]
        BCF["Mod_BC_Fluxes"]
        Wall["Lib_BC_Fluxes_Wall"]
        Conn["Lib_BC_Fluxes_Connection"]
        Newstate["Mod_Newstate / Lib_Newstate"]
        RK["Lib_RK"]
        IRS["Lib_IRS"]
        MG["Mod_Multigrid"]
        Metrics["Mod_Metrics"]
    end

    subgraph other ["physics/ · parallel/ · diagnostic/"]
        Solid["Lib_Solid"]
        GX["Mod_GhostExchange"]
        MPI["Mod_MPI"]
        Diag["Mod_Diagnostic"]
    end

    Procedures --> WSetup & WSolve & WPost
    WSetup --> LoadTable & ReadIni & Assign & IOSol & Alloc & MG & IOBC & Metrics & IOProbes & IOWall & Ghost & GX
    WSolve --> Explicit & MG
    WPost --> IOSol & IOWall & IOProbes & Diag & MG & ReadIni
    Assign --> Newstate
    Explicit --> DT & Ghost & Fluxes & BCF & Newstate & Diag
    Fluxes --> Diff
    Diff --> Solid
    BCF --> Wall & Conn
    Ghost --> Wall & GX
    IOWall --> Wall
    Newstate --> RK & IRS & Solid
    IRS --> Ghost
    GX --> MPI
```

---

## Naming

| Pattern | Use | Example |
|---------|-----|---------|
| `FUSS_` | Prefix of every module name | `FUSS_Mod_Fluxes` in `Mod_Fluxes.f90` |
| `Mod_*` | Mostly domain-level drivers and module state (procedure pointers, schedules) | `Mod_Fluxes`, `Mod_Newstate`, `Mod_GhostExchange` |
| `Lib_*` | Mostly routines on a block, face or cell | `Lib_Diffusive`, `Lib_BC_Fluxes_Wall`, `Lib_Solid` |
| `*_m` | Types and parameters | `Global_m`, `Config_Types_m` |
| `Wrap_*` | The three phases of `FUSS_type` | `Wrap_Setup` |
| `IO_*` | File input/output | `IO_Wall` |
| `obj_*` | Settings objects in `Config_Types_m` | `obj_time_scheme`, `obj_io` |
| `*_type` | Derived types | `FUSS_block_type` |

The module name is `FUSS_` plus the file name, except in `config/`: `Registry.f90` holds `FUSS_Input_Registry`, and `Register_IO.f90`, `Register_Numerics.f90`, `Register_Sim_Param.f90` hold `FUSS_Read_IO`, `FUSS_Read_Numerics`, `FUSS_Read_Sim_Param`. Routine names have no prefix.

---

## Solver pipeline

`src/app/main.f90` creates a `FUSS_type` object and a `FUSS_simulation_type`, calls `setup`, then alternates `solve` (one iteration) and `postprocess` while `obj_sim_param%TODO` ≤ 2, and calls `postprocess` once more at the end. `TODO` is 1 for a plain iteration, 2 when output is due or a multigrid level ends, 3 when the run ends.

```mermaid
sequenceDiagram
    participant Main as main.f90
    participant Setup as Wrap_Setup
    participant Solve as Wrap_Solve
    participant Post as Wrap_Postprocess

    Main->>Setup: FUSS%setup(simulation)
    activate Setup
    Setup->>Setup: Load_Table (phase.txt, properties.dat)
    Setup->>Setup: Read_Inifile (input.ini)
    Setup->>Setup: Assign_Setup (schemes, procedure pointers)
    Setup->>Setup: Read_IC (ic.tec, st.tec or restart file)
    Setup->>Setup: Setup_Data_Structure
    Setup->>Setup: Setup_Multigrid (if levels > 1)
    Setup->>Setup: Setup_BC (bc.txt, and bc2.txt ... on coarse levels)
    Setup->>Setup: partition_blocks, exchange schedules, Fill_matIDg
    Setup->>Setup: Setup_Metrics (every level)
    Setup->>Setup: Setup_Probes, Setup_Output_Solution, Initialize_Wall_File
    Setup->>Setup: Fill_Ghost_Cell, Restriction to coarse levels
    deactivate Setup

    loop while TODO ≤ 2
        Main->>Solve: FUSS%solve(simulation, External_Function)
        activate Solve
        Solve->>Solve: Explicit_Step (current level)
        Solve->>Solve: Prolongation (when a coarse level ends)
        deactivate Solve
        Main->>Post: FUSS%postprocess(simulation)
        Note right of Post: residual history, probes, screen,<br/>runtime input.ini, output files
    end

    Main->>Post: FUSS%postprocess(simulation)
    Note right of Post: final output and timing summary
```

The external subroutine passed to `solve` is called at every stage; `main.f90` passes an empty one, hydra uses it for the coupling.

### Explicit step

```mermaid
flowchart TB
    Start([Explicit_Step]) --> DT["Compute_dt<br/>local Δt, global minimum if time-accurate"]
    DT --> Copy["Copy_State (TO = T)"]

    subgraph RK ["each Runge–Kutta stage"]
        direction TB
        Ghost["Fill_Ghost_Cell<br/>(MPI ghost and chimera exchange)"] --> Diff["Fluxes<br/>interior faces"]
        Diff --> BCF["BC_Fluxes<br/>101/102, 301–305, 103/104"]
        BCF --> Ext["External_Function"]
        Ext --> Update["RK_Newstate<br/>source qvol, IRS, state update"]
    end

    Copy --> RK
    RK --> Res["Compute_Residual (every res-diter)<br/>MPI sum"]
    Res --> Flags["TODO and level switch on rank 0,<br/>broadcast"]
    Flags --> End([return])
```

### Interior face flux

`Fluxes_blk` (`Mod_Fluxes`) loops over the faces of each direction and calls `Diffusive_Flux` (`Lib_Diffusive`):

```mermaid
flowchart LR
    A["Stencil<br/>2 cells across the face<br/>+ 8 tangential neighbours"] --> D["Computational gradient<br/>normal: T₂ − T₁<br/>tangential: central differences"]
    A --> B["Face temperature<br/>½(T₁ + T₂)"]
    B --> C["k = co_Kappa(matID, T)<br/>(material of cell 1)"]
    D --> E["Face metric<br/>inverse of mean M⁻¹"]
    E --> G["∇T = Σ_d ∂T/∂ξ_d ∇ξ_d"]
    C --> F["Flux = k A ∇T·n<br/>R₁ −= Flux, R₂ += Flux"]
    G --> F
```

See [Spatial Discretisation](../theory/numerics.md) for the formulas.

---

## Data structures

```mermaid
classDiagram
    class FUSS_simulation_type {
        +domain(:) FUSS_domain_type
        +IOfield(:) ORION_data
    }
    class FUSS_domain_type {
        +time, dtglobal : real64
        +iter, itermax, mg_level : integer
        +nb, nbound : integer
        +n_bf(:,:) : integer
        +blk(:) FUSS_block_type
        +bc(:) FUSS_bc_type
    }
    class FUSS_block_type {
        +dim : integer
        +node, M, dl, vol, dir, fm, bm
        +T, TO, R, qvol, matID : real64
        +dtlocal : real64
        +RS0, RS1, RS2 : real64
    }
    class FUSS_bc_type {
        +b, i, j, k, f, type : integer
        +bs, is, js, ks, fs, d11..d22 : integer
        +Tg, matIDg, Mg
        +ni, donorID, volume_fraction
        +ext_flux
        +qw, Tw, hconv, qrad, Tref, eps
        +BCtime : time_series_type
    }

    FUSS_simulation_type "1" --> "levels" FUSS_domain_type : domain
    FUSS_domain_type "1" --> "nb" FUSS_block_type : blk
    FUSS_domain_type "1" --> "nbound" FUSS_bc_type : bc
```

`domain(m)` is multigrid level `m` (1 = finest). `bc(:)` has one entry per boundary cell face, in the order of `bc.txt`.

| Array | Bounds | Content |
|-------|--------|---------|
| `T`, `TO` | `1-gc:n+gc` | Temperature, and its value at the start of the iteration |
| `R` | `1-gc:n+gc` | Residual (sum of face fluxes) |
| `qvol`, `matID` | `1-gc:n+gc` | Volumetric source; material index (stored as real) |
| `dtlocal` | `1:n` | Local time step |
| `RS0`, `RS1`, `RS2` | `1-gc:n+gc` | IRS work arrays (allocated with IRS only) |
| `M`, `dl`, `vol` | `1-gc:n+gc` | Metric (rows ∇ξ_d), cell lengths, volume |
| `fm(d)%Mf` | `1:n-1` along `d`, `1:n` in the others | Stencil face metric (3×3) of the interior face between a cell and the next one along `d` (`diffusive-metric = stencil`) |
| `bm(f)` | six block faces | Stencil metrics of the walls (`wM`) and of the block connections (`cM`) that the record owns |
| `node` | `0:n` | Grid nodes |
| `dir(d)%f` | faces normal to `d` | Unit normal `N` and area `A` |

---

## Build system

| Target | Type | Source |
|--------|------|--------|
| `FUSSL` (alias `FUSSL::FUSSL`) | Library (static by default) | `src/lib/**/*.f90` |
| `FUSS` | Executable, `bin/FUSS` | `src/app/main.f90` |
| `DocGen` | Executable, `bin/DocGen` | `src/app/docgen.f90` |

| Variable | Default | Effect |
|----------|---------|--------|
| `USE_OPENMP` | `OFF` | OpenMP |
| `USE_MPI` | `OFF` | MPI; defines the preprocessor symbol `USE_MPI` |
| `USE_TECIO` | `ON` (ORION option) | TecIO; `install.sh` passes `false` unless `--use-tecio` |
| `CMAKE_BUILD_TYPE` | `RELEASE` | `RELEASE`, `DEBUG` or `TESTING` (`cmake/SetFortranFlags.cmake`) |
| `ORION_PATH`, `FINER_PATH`, `FINER_DEPS_PATH` | `lib/ORION/`, `lib/third_party/FiNeR/`, `lib/third_party/` | Dependency source trees |

When FUSS is added to another CMake project (as in hydra) only `FUSSL` is built; the parent defines ORION and FiNeR. The executables and the tests are built only when FUSS is the top-level project. See [Installation](../getting-started/installation.md).
