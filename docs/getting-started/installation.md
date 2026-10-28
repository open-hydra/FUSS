# Installation

FUSS is built with CMake, either through the `install.sh` script or by calling CMake directly. A build produces the static library `libFUSSL.a` and two executables in `bin/`: `FUSS` (the solver) and `DocGen` (writes the [input parameter table](../user/registry.md)).

## Prerequisites

- **CMake** 3.23 or newer.
- **Fortran compiler**: GNU (`gfortran`) or Intel oneAPI (`ifx`).
- **C and C++ compilers**: the project enables C and C++ (ORION and TecIO).
- **git**: the dependencies are Git submodules (`git clone --recurse-submodules`, or `git submodule update --init` after pulling a change of the pins).
- **OpenMP** and **MPI**: optional.

### Dependencies

| Path | Repository | Purpose |
|------|------------|---------|
| `lib/ORION` | [MarcoGrossi92/ORION](https://github.com/MarcoGrossi92/ORION) (submodule) | Tecplot/VTK reading and writing; bundles TecIO |
| `lib/third_party/FiNeR` | [szaghi/FiNeR](https://github.com/szaghi/FiNeR) (submodule) | `input.ini` parser |
| `lib/third_party/{PENF,FACE,FLAP,BeFoR64,StringiFor}` | [szaghi](https://github.com/szaghi) (submodules) | FiNeR's own libraries, added by CMake before FiNeR; pinned at the tags of MOSE (PENF v2.0.1, FACE v1.1.17, FLAP v2.5.3, BeFoR64 v1.2.4, StringiFor v1.8.0) |

### ifx 2024.1 and IPO

ifx 2024.1 (tested: 2024.1.0, build 20240308) compiles a loop of the ghost metrics of MOSE wrongly at `-O2` and `-O3` without interprocedural optimization: a minimal program in standard Fortran gives wrong values at `-O2` and `-O3`, also at the default level with no `-O` option, and the right ones at `-O0`, `-O1` and with `-ipo`. FUSS built without IPO (the RELEASE build of ifx up to 2026-10-28, whose `-ipo` was dropped silently; see the CHANGELOG) passes every case of the stress battery and of the patch tests, so no miscompiled loop is known in FUSS; CMake refuses that build all the same, as MOSE does, because the pattern is ordinary Fortran. The refused builds are the RELEASE build with `-DUSE_IPO=OFF` (`--no-ipo`) or whose toolchain cannot do IPO, and the TESTING build (`-O2`); the RELEASE build with IPO (the default) and the DEBUG build (`-O0`) are accepted. `-DFUSS_ALLOW_IFX_NOIPO=ON` turns the refusal into a warning. Other versions of ifx have not been tested.

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

1. Initialises the submodules `lib/ORION`, `lib/third_party/FiNeR` and FiNeR's five libraries under `lib/third_party/` (ORION and FiNeR unless an external path is given).
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
| `--include-orion=PATH` | Use an external ORION tree instead of the submodule |
| `--include-finer=PATH` | Use an external FiNeR tree instead of the submodule |

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
./install.sh update          # submodules (ORION, FiNeR and its five libraries)
cmake -B build \
      -DORION_PATH=$PWD/lib/ORION/ \
      -DFINER_PATH=$PWD/lib/third_party/FiNeR/ \
      -DUSE_OPENMP=true -DUSE_MPI=false -DUSE_TECIO=false \
      -DCMAKE_BUILD_TYPE=RELEASE
cmake --build build
```

- `ORION_PATH` and `FINER_PATH` default to `lib/ORION/` and `lib/third_party/FiNeR/`; `FINER_DEPS_PATH` (default `lib/third_party/`) is where CMake finds PENF, FACE, FLAP, BeFoR64 and StringiFor, which it adds before FiNeR. A submodule that is not initialised is initialised by CMake (`git submodule update --init`); FiNeR older than v2.1.0 stops the configuration.
- `USE_IPO` (default ON) adds the interprocedural optimization (`-ipo` with Intel; `-flto` with GNU) to the RELEASE flags (results of a RELEASE build with ifx change at round-off against builds before 2026-10-28, see the CHANGELOG); with ifx the archiver is the `llvm-ar` and `llvm-ranlib` of the compiler (`-DIPO_AR`, `-DIPO_RANLIB` give them for another layout of the installation), and IPO is dropped with a warning when they are not found. `install.sh` writes `CMAKE_AR`, `CMAKE_RANLIB` and `USE_IPO` into `CMakePresets.json`.
- `USE_TECIO` is an ORION option and defaults to ON there: pass `-DUSE_TECIO=false` to build without TecIO.
- `CMAKE_BUILD_TYPE` is `RELEASE` (default), `DEBUG` (`-O0`, warnings, traceback, bounds checking) or `TESTING` (`-O2`); see `cmake/SetFortranFlags.cmake`.
- Use the `build/` directory: ORION, FiNeR and its libraries are built in `lib/` of the build directory (`-B` elsewhere leaves the source tree clean), and the executables always go to `bin/`.

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
