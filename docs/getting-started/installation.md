# Installation

This document describes how to obtain and build **FUSS**. The instructions describe the current `install.sh` script, CMake configuration and Git submodule layout.

!!! note
    FUSS has a dual nature: it is both a library and an executable. The installation process described here produces both the static library `libFUSSL.a` and the main executable `bin/FUSS`. If you are only interested in using FUSS as a library, you can link against `libFUSSL.a` without caring about the executable.

## Prerequisites

Before attempting to build FUSS make sure your system provides the following external tools and compilers:

- **CMake** – 3.23 or newer.
- **Fortran compiler** – either the GNU toolchain (`gfortran`) or Intel/oneAPI (`ifort`/`ifx`) are supported.
- **C / C++ compiler** – required by the ORION I/O library and the optional TecIO component.
- **OpenMP** – needed for optional shared-memory parallelisation.
- **MPI** – needed for optional distributed-memory parallelisation.

### Git submodules

FUSS depends on three repositories that are included as Git submodules.

| Path | Repository URL | Purpose |
|------|----------------|---------|
| `lib/ORION` | `https://github.com/MarcoGrossi92/ORION.git` | I/O routines (TecIO, VTK) |
| `lib/third_party/FiNeR` | `https://github.com/szaghi/FiNeR.git` | INI file parser |
| `lib/MORPH` | `https://github.com/open-hydra/MORPH.git` | Mesh motion (MORPH; FUSS links its library `MORPHL`) |

## Build methods

First clone the repository with submodules:

```bash
git clone https://github.com/open-hydra/FUSS.git
cd FUSS
# initialise submodules
git submodule update --init --recursive
```

To fully install FUSS, you may either use the bundled install script or invoke CMake manually. The script is convenient and is the preferred route for most users.

### Build with `install.sh` (recommended)

The script exposes three commands: `build`, `compile`, and `update`. It also maintains a `CMakePresets.json` file that records the configuration used for the most recent `build` invocation.

```bash
./install.sh [GLOBAL_OPTIONS] COMMAND [COMMAND_OPTIONS]
```

**Global options**

* `-v`, `--verbose` – enable verbose logging.

**`build` command**

Performs a clean configure+build cycle. Example usage:

```bash
# minimal GNU build with OpenMP enabled
./install.sh build --compilers=gnu --use-openmp

# full configuration with Intel compilers and all optional features
./install.sh build --compilers=intel \
                   --use-openmp --use-mpi \
                   --use-tecio
```

Options accepted by `build`:

* `--compilers=<gnu|intel>` – select the compiler family (default: `gnu`).
* `--use-openmp` – enable OpenMP parallelisation.
* `--use-mpi` – enable MPI parallelisation.
* `--use-tecio` – enable TecIO support (requires a C++ compiler).
* `--include-orion=PATH` – use an external ORION tree instead of the submodule.
* `--include-finer=PATH` – same for FiNeR.
* `--include-morph=PATH` – same for MORPH.

The script sets up environment variables for the chosen compilers and then invokes CMake. After a successful build, a `CMakePresets.json` file is written in the source root so that subsequent compilations can reuse the configuration.

**`compile` command**

Re-runs CMake using the previously generated preset and rebuilds the project without clearing the build directory. This is useful during development when only the source has changed. Example usage:

```bash
./install.sh compile
```

**`update` command**

Synchronises the Git submodules. By default it checks out the commit recorded in `.gitmodules`; passing `--remote` will fetch the latest commit from each remote branch.

```bash
./install.sh update            # sync to recorded commit
./install.sh update --remote   # update to newest remote commit
```

### Build with CMake

If you prefer fine-grained control, perform the configuration yourself. This is essentially what `install.sh` does under the hood.

```bash
mkdir build && cd build
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_Fortran_COMPILER=gfortran \  # or ifx
    -DUSE_OPENMP=ON \                    # or OFF
    -DUSE_MPI=OFF \                      # optional
    -DUSE_TECIO=ON \                     # optional (needs C++ compiler)
    -DORION_PATH=/path/to/ORION \        # optional
    -DFINER_PATH=/path/to/FiNeR \        # optional
    -DMORPH_PATH=/path/to/MORPH          # optional (default: lib/MORPH)
cmake --build . --parallel
```

The resulting artifacts are placed in `build/` by default. The static library is `lib/libFUSSL.a` and the main executable is `bin/FUSS` (inside the build directory unless you set `CMAKE_INSTALL_PREFIX`).

## CMake presets

The file `CMakePresets.json` produced by the install script records all of the cache variables that were used during configuration. You can build the project later simply with

```bash
cmake --preset default
cmake --build build
```

or using the `compile` command of the install script as described above.

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
- `test/check_build_config.sh --preset mpi` checks `build-mpi/` against the preset; without `--preset` it checks
  `build/` against `default`.

## Optional components

### TecIO

FUSS can be built with support for TecIO, a library for writing Tecplot binary files. This is an optional feature shipped with ORION. If enabled, FUSS will be able to write solution output in Tecplot binary format (`.szplt`), which is useful for post-processing and visualisation of large datasets where ASCII Tecplot files would be impractical. Enabling TecIO requires a working C++ compiler.

### MPI

When MPI is enabled, FUSS distributes whole blocks across ranks and exchanges ghost-cell layers across block boundaries at every time step. MPI is most useful for multi-block cases where the block count is greater than or equal to the rank count.

## Library linking (advanced)

To use FUSS from an external Fortran program you can compile as follows:

```bash
gfortran -I/path/to/FUSS/include \
         -L/path/to/FUSS/lib \
         -lFUSSL \
         your_program.f90 -o your_program
```

or, from a CMake project:

```cmake
find_package(FUSS REQUIRED)
add_executable(myapp main.f90)
target_link_libraries(myapp FUSS::FUSSL)
```

Installation prefix and other details may be customised via standard CMake variables such as `CMAKE_INSTALL_PREFIX` and by running `cmake --install`.

## Next steps

* **[Quick Start Tutorial](quick-start.md)** – build and verify the installation.

---