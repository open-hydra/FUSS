# Installation

FUSS is built with CMake, either through the `install.sh` script or by calling CMake directly. A build produces the static library `libFUSSL.a` and two executables in `bin/`: `FUSS` (the solver) and `DocGen` (writes the [input parameter table](../user/registry.md)).

## Prerequisites

- **CMake** 3.23 or newer.
- **Fortran compiler**: GNU (`gfortran`) or Intel oneAPI (`ifx`).
- **C and C++ compilers**: the project enables C and C++ (ORION and TecIO).
- **git**: FiNeR and its libraries are Git submodules (`git clone --recurse-submodules`, or `git submodule update --init` after pulling a change of the pins); the configuration checks out a release of ORION (over https, the first time).
- **OpenMP** and **MPI**: optional.

### Dependencies

| Path | Repository | Purpose |
|------|------------|---------|
| `lib/ORION` | [MarcoGrossi92/ORION](https://github.com/MarcoGrossi92/ORION) (a release from 2.1 on, below 3.0, checked out by the configuration; not tracked) | Tecplot/VTK reading and writing; bundles TecIO |
| `lib/third_party/FiNeR` | [szaghi/FiNeR](https://github.com/szaghi/FiNeR) (submodule) | `input.ini` parser |
| `lib/third_party/{PENF,FACE,FLAP,BeFoR64,StringiFor}` | [szaghi](https://github.com/szaghi) (submodules) | FiNeR's own libraries, added by CMake before FiNeR; pinned at the tags of MOSE (PENF v2.0.1, FACE v1.1.17, FLAP v2.5.3, BeFoR64 v1.2.4, StringiFor v1.8.0) |

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

1. Initialises the submodules `lib/third_party/FiNeR` and FiNeR's five libraries under `lib/third_party/` (unless an external FiNeR path is given).
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
| `--include-finer=PATH` | Use an external FiNeR tree instead of the submodule |
| `--update-deps` | Select again the newest release of ORION (`-DSEMVER_UPDATE=ALL`; see below) |

```bash
./install.sh build --compilers=intel --use-openmp --use-mpi --use-tecio
```

### `compile`

Rebuilds without deleting `build/`, with the configuration stored in `CMakePresets.json` (so it needs one earlier `build`):

```bash
./install.sh compile        # cmake --preset default && cmake --build build
```

### `update`

Runs `git submodule update --init` (with `--remote`: moves each submodule to the latest commit of its remote branch instead of the commit recorded in FUSS), FiNeR's libraries included.

```bash
./install.sh update
./install.sh update --remote
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
./install.sh update          # submodules (FiNeR and its five libraries)
cmake -B build \
      -DFINER_PATH=$PWD/lib/third_party/FiNeR/ \
      -DUSE_OPENMP=true -DUSE_MPI=false -DUSE_TECIO=false \
      -DCMAKE_BUILD_TYPE=RELEASE
cmake --build build
```

- **ORION** is required as a 2.x release from 2.1 on (`cmake/SemVerDeps.cmake`). The first configuration checks out the newest one in `lib/ORION`, and `deps.lock` (not tracked) records it, so later configurations need no network.
    - `-DSEMVER_UPDATE=ALL` selects the newest release again (once).
    - `-DFETCHCONTENT_SOURCE_DIR_ORION=PATH` builds a local checkout instead.
    - `-DFETCHCONTENT_FULLY_DISCONNECTED=ON` forbids any clone or fetch.
- **FiNeR**: `FINER_PATH` defaults to `lib/third_party/FiNeR/`, and `FINER_DEPS_PATH` (default `lib/third_party/`) holds PENF, FACE, FLAP, BeFoR64 and StringiFor. CMake initialises a submodule that is not initialised. FiNeR older than v2.1.0 stops the configuration.
- **`USE_IPO`** (default ON) adds interprocedural optimization to the RELEASE flags (`-ipo` with Intel, `-flto` with GNU). With ifx the archiver is the compiler's `llvm-ar` and `llvm-ranlib` (`-DIPO_AR`, `-DIPO_RANLIB` for another installation layout); without them IPO is dropped with a warning. `install.sh` stores `CMAKE_AR`, `CMAKE_RANLIB` and `USE_IPO` in `CMakePresets.json`.
- **`USE_TECIO`** is an ORION option and defaults to ON there: pass `-DUSE_TECIO=false` to build without TecIO.
- **`CMAKE_BUILD_TYPE`** is `RELEASE` (default), `DEBUG` (`-O0`, warnings, traceback, bounds checking) or `TESTING` (`-O2`); see `cmake/SetFortranFlags.cmake`.
- Use the `build/` directory: FiNeR and its libraries are built in `lib/` of the build directory, ORION in `_deps/` (`-B` elsewhere leaves the source tree clean), and the executables always go to `bin/`.

The configuration stored by `install.sh` can be reused with:

```bash
cmake --preset default
cmake --build build
```

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
add_subdirectory(lib/third_party/FiNeR)
add_subdirectory(FUSS)
target_link_libraries(myapp PRIVATE FUSSL FiNeR::FiNeR)
```

`src/app/main.f90` is the reference driver: it calls `setup`, `solve` and `postprocess` of `FUSS_type` (module `FUSS_Procedures_m`).

## Next steps

- [Quick Start](quick-start.md): run and check a first case.
