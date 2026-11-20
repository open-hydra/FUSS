# Installation

FUSS is built with CMake, either through the `install.sh` script or by calling CMake directly. A build produces the static library `libFUSSL.a` and two executables in `bin/`: `FUSS` (the solver) and `DocGen` (writes the [input parameter table](../user/registry.md)).

## Prerequisites

- **CMake** 3.23 or newer.
- **Fortran compiler**: GNU (`gfortran`) or Intel oneAPI (`ifx`).
- **C and C++ compilers**: the project enables C and C++ (ORION and TecIO).
- **git**: FiNeR and MORPH are Git submodules (`git clone --recurse-submodules`, or `git submodule update --init` after pulling a change of their pins); the configuration clones FiNeR's own libraries and checks out a release of ORION (over https, the first time).
- **OpenMP** and **MPI**: optional.

### Dependencies

| Path | Repository | Purpose |
|------|------------|---------|
| `lib/ORION` | [MarcoGrossi92/ORION](https://github.com/MarcoGrossi92/ORION) (a release from 2.1 on, below 3.0, checked out by the configuration; not tracked) | Tecplot/VTK reading and writing; bundles TecIO |
| `lib/FiNeR` | [szaghi/FiNeR](https://github.com/szaghi/FiNeR) (submodule, v2.3.0) | `input.ini` parser |
| `lib/FiNeR/src/third_party/{PENF,FACE,FLAP,BeFoR64,StringiFor}` | [szaghi](https://github.com/szaghi) (cloned by the configuration) | FiNeR's own libraries, at their default branches (`FUSS_FINER_DEPS=FREE`) or at fixed releases (`BACKUP`) |
| `lib/MORPH` | [open-hydra/MORPH](https://github.com/open-hydra/MORPH) (submodule) | Mesh motion; FUSS links its library `MORPHL` |

The configuration checks ORION and FiNeR against the releases FUSS requires (see [Versioning](../development/versioning.md)); `lib/MORPH` is not checked.

!!! warning "ifx 2024.1 needs IPO"
    ifx 2024.1 miscompiles an ordinary Fortran loop pattern (found in MOSE) at `-O2` and `-O3`
    without interprocedural optimization (IPO). No such loop is known in FUSS, but CMake refuses
    the affected ifx 2024.1 builds, as MOSE does:

    - refused: RELEASE with `--no-ipo` (`-DUSE_IPO=OFF`) or with a toolchain that cannot do IPO,
      and TESTING (`-O2`);
    - accepted: RELEASE with IPO (the default) and DEBUG (`-O0`).

    `-DFUSS_ALLOW_IFX_NOIPO=ON` turns the refusal into a warning. Other ifx versions are not
    tested.

## Build with `install.sh`

```bash
git clone https://github.com/open-hydra/FUSS.git
cd FUSS
./install.sh build --compilers=gnu --use-openmp
```

```bash
./install.sh [GLOBAL_OPTIONS] COMMAND [COMMAND_OPTIONS]
```

Global options: `-h`, `--help` (usage) and `-v`, `--verbose` (log the configuration).

### `build`

1. Syncs the submodule URLs and initialises `lib/FiNeR` and `lib/MORPH` (each unless an external clone is given). The configuration clones FiNeR's five libraries.
2. Sets the compilers, deletes `build/`, configures CMake in `build/` with `CMAKE_BUILD_TYPE=RELEASE` and builds.
3. Writes `CMakePresets.json` with the configuration just used.

| Option | Effect |
|--------|--------|
| `--compilers=gnu` | `gfortran`, `gcc`, `g++` |
| `--compilers=intel` | `ifx`, `icx`, `icpx` |
| (none) | CMake chooses the compilers (`FC`, `CC`, `CXX` from the environment if set) |
| `--use-openmp` | OpenMP |
| `--use-mpi` | MPI |
| `--use-tecio` | TecIO, for Tecplot binary (`.szplt`) files |
| `--no-ipo` | Without interprocedural (link-time) optimization in the RELEASE build (`-DUSE_IPO=OFF`) |
| `--include-orion=PATH` | Build a local ORION checkout instead of a release (`-DFETCHCONTENT_SOURCE_DIR_ORION`) |
| `--include-finer=PATH` | Build an external FiNeR clone instead of the submodule `lib/FiNeR` (`-DFETCHCONTENT_SOURCE_DIR_FINER`), checked against `FUSS_FINER_TAG` |
| `--include-morph=PATH` | Use an external MORPH tree instead of the submodule |
| `--update-deps` | Select again the newest release of ORION (`-DSEMVER_UPDATE=ALL`; see below) |

```bash
./install.sh build --compilers=intel --use-openmp --use-mpi --use-tecio
```

### `compile`

Rebuilds without deleting `build/`, with the configuration stored in `CMakePresets.json` (so it needs one earlier `build`):

```bash
./install.sh compile        # cmake --preset default && cmake --build build
```

### Build output

| File | Content |
|------|---------|
| `bin/FUSS` | Solver |
| `bin/DocGen` | Writes `docs/user/registry.md` (run it from the repository root) |
| `build/src/lib/libFUSSL.a` | FUSS library |
| `build/modules/` | Fortran module files |

## Build with CMake

`install.sh build` runs the equivalent of:

```bash
git submodule sync --quiet && git submodule update --init lib/FiNeR lib/MORPH
cmake -B build \
      -DUSE_OPENMP=true -DUSE_MPI=false -DUSE_TECIO=false \
      -DCMAKE_BUILD_TYPE=RELEASE
cmake --build build
```

- **ORION** is required as a 2.x release from 2.1 on (`cmake/SemVerDeps.cmake`). The first configuration checks out the newest one in `lib/ORION`, and `deps.lock` (not tracked) records it, so later configurations need no network.
    - `-DSEMVER_UPDATE=ALL` selects the newest release again (once).
    - `-DFETCHCONTENT_SOURCE_DIR_ORION=PATH` builds a local checkout instead.
    - `-DFETCHCONTENT_FULLY_DISCONNECTED=ON` forbids any clone or fetch.
- **FiNeR** is the submodule `lib/FiNeR` (initialised by CMake when it is not), or the clone of `-DFETCHCONTENT_SOURCE_DIR_FINER=PATH`. The configuration clones its five libraries into `lib/FiNeR/src/third_party/`. A FiNeR off the commit or the release FUSS pins stops the configuration ([Versioning](../development/versioning.md)).
- `MORPH_PATH` defaults to the submodule `lib/MORPH/`; only its library `MORPHL` is built.
- **`USE_IPO`** (default ON) adds interprocedural optimization to the RELEASE flags (`-ipo` with Intel, `-flto` with GNU). The archiver must index the optimised objects: with ifx it is the compiler's `llvm-ar` and `llvm-ranlib` (`-DIPO_AR`, `-DIPO_RANLIB` for another installation layout), with gfortran `gcc-ar` and `gcc-ranlib`. Without them IPO is dropped with a warning. `install.sh` stores `CMAKE_AR`, `CMAKE_RANLIB` and `USE_IPO` in `CMakePresets.json`.
- **`USE_TECIO`** is an ORION option and defaults to ON there: pass `-DUSE_TECIO=false` to build without TecIO.
- **`CMAKE_BUILD_TYPE`** is `RELEASE` (default), `DEBUG` (`-O0`, warnings, traceback, bounds checking) or `TESTING` (`-O2`); see `cmake/SetFortranFlags.cmake`.
- Use the `build/` directory: FiNeR and ORION are built in `_deps/` of the build directory (`-B` elsewhere leaves the source tree clean), and the executables always go to `bin/`.

The configuration stored by `install.sh` can be reused with:

```bash
cmake --preset default
cmake --build build
```

### The `mpi` preset: a second build tree

With gfortran and an Intel MPI in the environment (`I_MPI_ROOT` set, the gfortran modules in
`$I_MPI_ROOT/include/gfortran/11.1.0`), `install.sh build` writes a second preset, `mpi`. It inherits `default`,
configures its own tree, `build-mpi/`, and builds executables of its own, `bin/FUSS-mpi` and `bin/DocGen-mpi`
(cache variable `FUSS_EXE_SUFFIX`), so the serial `bin/FUSS` and the MPI binary coexist:

```bash
./install.sh compile --preset=mpi     # or: cmake --preset mpi && cmake --build build-mpi
```

- The preset bypasses Intel's `mpif90` wrapper, which gives gfortran 12 no module directory, and links Intel MPI's
  Fortran library directly. `CMAKE_BUILD_RPATH` lets `bin/FUSS-mpi` find it without `LD_LIBRARY_PATH`.
- With MPI the Fortran sources are compiled without `-fno-underscoring` (the MPI library exports the underscored
  names) and with `-fallow-argument-mismatch` (gfortran 10 and later reject one MPI routine called with a scalar and
  with an array).
- In this build and in the serial one, `src/lib/numerics/space/Lib_Stencil_Metrics.f90` is compiled with
  `-ffp-contract=off` (gfortran): the two builds contracted its face-metric arithmetic into fused multiply-adds
  differently, and the MPI binary on one rank then differed from the serial binary in the last bits on skewed meshes.
- `test/check_build_config.sh --preset mpi` checks `build-mpi/` against the preset; without `--preset` it checks
  `build/` against `default`.

## Optional components

### TecIO

TecIO writes Tecplot binary files (`.szplt`, with `sol-format = tecplot binary`). Its sources ship with ORION, whose CMake builds it at configure time (the MPI variant when `USE_MPI` is on).

### OpenMP and MPI

OpenMP shares the cell loops of each block among threads; MPI distributes whole blocks over the ranks. See [Parallel Execution](../user/parallel.md).

## Testing the build

```bash
ctest --test-dir build -L fast
```

See [Testing](../development/testing.md) for the labels and the cases.

## Library linking (advanced)

`src/lib/CMakeLists.txt` has no install rules (they are commented out), so `cmake --install` and `find_package(FUSS)` are not available. A program that uses FUSS as a library adds the FUSS source tree to its own CMake project, as hydra does: the parent project defines the ORION and FiNeR targets (FUSS adds them itself only when it is the top-level project) and links the `FUSSL` target.

```cmake
add_subdirectory(lib/ORION)
add_subdirectory(lib/FiNeR)
add_subdirectory(FUSS)
target_link_libraries(myapp PRIVATE FUSSL FiNeR::FiNeR)
```

`src/app/main.f90` is the reference driver: it calls `setup`, `solve` and `postprocess` of `FUSS_type` (module `FUSS_Procedures_m`).

## Next steps

- [Quick Start](quick-start.md): run and check a first case.
