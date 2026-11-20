#!/bin/bash

set -e  # Exit on any command failure
set -u  # Treat unset variables as an error

PROGRAM=$(basename "$0")
readonly DIR=$(pwd)
VERBOSE=false
BUILD_DIR="$DIR/build"
project=FUSS

function usage() {
    cat <<EOF

Install script for $project

Usage:
  $PROGRAM [GLOBAL_OPTIONS] COMMAND [COMMAND_OPTIONS]

Global Options:
  -h       , --help         Show this help message and exit
  -v       , --verbose      Enable verbose output

Commands:
  build                     Perform a full build
    --include-orion=<path>  Build a local ORION checkout instead of a release
    --include-finer=<path>  Build an external FiNeR clone instead of the submodule lib/FiNeR
    --include-morph=<path>  Set external MORPH path
    --update-deps           Select again the newest release of ORION (deps.lock)
    --compilers=<name>      Set compilers (intel, gnu)
    --use-openmp            Use OpenMP
    --use-mpi               Use MPI
    --use-tecio             Use TecIO
    --no-ipo                Disable interprocedural (link-time) optimization

  compile                   Compile the program using the CMakePresets file

EOF
    exit 1
}

log() {
    if [ "$VERBOSE" = true ]; then
        # Bold and dim gray (ANSI escape: bold + color 90)
        echo -e "\033[1;90m$1\033[0m"
    fi
}

error() {
    # Bold red + [ERROR] tag, output to stderr
    echo -e "\033[1;31m[ERROR] $1\033[0m" >&2
}

task() {
    # Bold yellow + ==> tag, output to stdout
    echo -e "\033[1;38;5;186m==> $1\033[0m"
}


# Create default CMakePresets.json if it doesn't exist
function write_presets() {
  FC=$(grep -E '^CMAKE_Fortran_COMPILER:(FILEPATH|STRING)=' "$BUILD_DIR/CMakeCache.txt" | cut -d= -f2-)
  CC=$(grep -E '^CMAKE_C_COMPILER:(FILEPATH|STRING)=' "$BUILD_DIR/CMakeCache.txt" | cut -d= -f2-)
  CXX=$(grep -E '^CMAKE_CXX_COMPILER:(FILEPATH|STRING)=' "$BUILD_DIR/CMakeCache.txt" | cut -d= -f2-)
  AR=$(grep '^CMAKE_AR:FILEPATH=' "$BUILD_DIR/CMakeCache.txt" | cut -d= -f2-)
  RANLIB=$(grep '^CMAKE_RANLIB:FILEPATH=' "$BUILD_DIR/CMakeCache.txt" | cut -d= -f2-)
  # the local ORION and FiNeR checkouts, when given
  local sources=""
  if [[ -n $ORION_PATH ]]; then
    sources+="        \"FETCHCONTENT_SOURCE_DIR_ORION\": \"${ORION_PATH}\","$'\n'
  fi
  if [[ -n $FINER_PATH ]]; then
    sources+="        \"FETCHCONTENT_SOURCE_DIR_FINER\": \"${FINER_PATH}\","$'\n'
  fi

  # The mpi preset (plan 11 section 2.1): a second tree, build-mpi/, with its own executable
  # bin/FUSS-mpi. gfortran + Intel MPI with the compiler wrapper bypassed: Intel's mpif90
  # maps gfortran 12 to no module directory and `use mpi` reads the ifx mpi.mod; the gfortran
  # 11.1 modules load under gfortran 12 (measured 2026-10-06 and 2026-10-15 on sprop2; FUSS
  # needs `use mpi` only: include/ itself holds the ifx modules, never on the path). Written
  # only with gfortran and an Intel MPI ($I_MPI_ROOT) that ships them.
  MPI_PRESET=""
  if [[ -n "${I_MPI_ROOT:-}" && -d "$I_MPI_ROOT/include/gfortran/11.1.0" && "$FC" == *gfortran* ]]; then
    MPI_PRESET=$(cat <<EOM
,
    {
      "name": "mpi",
      "description": "MPI preset: Intel MPI, wrapper bypassed (build-mpi/, bin/FUSS-mpi)",
      "inherits": "default",
      "binaryDir": "\${sourceDir}/build-mpi",
      "cacheVariables": {
        "USE_MPI": "true",
        "FUSS_EXE_SUFFIX": "-mpi",
        "MPI_SKIP_COMPILER_WRAPPER": "TRUE",
        "MPI_Fortran_ADDITIONAL_INCLUDE_DIRS": "$I_MPI_ROOT/include/gfortran/11.1.0",
        "MPI_Fortran_LIB_NAMES": "mpifort;mpi",
        "MPI_mpifort_LIBRARY": "$I_MPI_ROOT/lib/libmpifort.so",
        "MPI_mpi_LIBRARY": "$I_MPI_ROOT/lib/release/libmpi.so",
        "MPIEXEC_EXECUTABLE": "$I_MPI_ROOT/bin/mpiexec",
        "CMAKE_BUILD_RPATH": "$I_MPI_ROOT/lib/release;$I_MPI_ROOT/lib"
      }
    }
EOM
)
  else
    log "No mpi preset: it needs gfortran and an Intel MPI (I_MPI_ROOT) with include/gfortran/11.1.0"
  fi

  cat <<EOF > CMakePresets.json
{
  "version": 3,
  "cmakeMinimumRequired": {
    "major": 3,
    "minor": 23
  },
  "configurePresets": [
    {
      "name": "default",
      "description": "Default preset",
      "binaryDir": "\${sourceDir}/build",
      "cacheVariables": {
${sources}        "MORPH_PATH": "${MORPH_PATH}",
        "CMAKE_BUILD_TYPE": "${BUILD_TYPE}",
        "CMAKE_Fortran_COMPILER": "${FC}",
        "CMAKE_C_COMPILER": "${CC}",
        "CMAKE_CXX_COMPILER": "${CXX}",
        "CMAKE_AR": "${AR}",
        "CMAKE_RANLIB": "${RANLIB}",
        "USE_TECIO": "${USE_TECIO}",
        "USE_OPENMP": "${USE_OPENMP}",
        "USE_MPI": "${USE_MPI}",
        "USE_IPO": "${USE_IPO}"
      }
    }${MPI_PRESET}
  ]
}
EOF
}


# Default global values
COMMAND=""
COMPILERS=""
# ORION: the release of deps.lock, or the newest one that FUSS accepts (see CMakeLists.txt);
# FiNeR: the submodule lib/FiNeR, or an external clone
ORION_PATH=""
UPDATE_DEPS="false"
FINER_PATH=""
MORPH_PATH=$(pwd)'/lib/MORPH/'
USE_OPENMP="false"
USE_MPI="false"
USE_TECIO="false"
# IPO is on by default; CMake self-disables it if the toolchain cannot do it.
USE_IPO="true"
REMOTE="false"
BUILD_TYPE="RELEASE"

# Define allowed options for each command using regular arrays
CMD=("build" "compile")
CMD_OPTIONS_build=("--compilers --include-orion --include-finer --include-morph --update-deps --use-openmp --use-mpi --use-tecio --no-ipo")
CMD_OPTIONS_compile=("--preset=NAME")

# Parse global options
while getopts "hv-:" opt; do
    case "$opt" in
        -)
            case "$OPTARG" in
                verbose) VERBOSE=true ;;
                help) usage ;;
                *) error "Unknown global option '--$OPTARG'"; usage ;;
            esac
            ;;
        h) usage ;;
        v) VERBOSE=true ;;
        ?) error "Unknown global option '-$OPTARG'"; usage ;;
    esac
done
shift $((OPTIND -1))

# Ensure a command was provided
if [[ $# -eq 0 ]]; then
    error "No command provided!"
    usage
fi

COMMAND="$1"
# Check if the command is valid
if [[ ! " ${CMD[@]} " =~ " ${COMMAND} " ]]; then
    error "Unknown command '$COMMAND'"
    usage
fi
shift

# Parse command-specific options
while [[ $# -gt 0 ]]; do
    case "$1" in
        --include-orion=*)
            [[ "$COMMAND" == "build" ]] || { error " --include-orion is only valid for 'build' command"; exit 1; }
            # absolute: CMake takes the folder from the build directory
            ORION_PATH=$(cd "${1#*=}" && pwd) || { error "No folder ${1#*=}"; exit 1; }
            ;;
        --include-finer=*)
            [[ "$COMMAND" == "build" ]] || { error " --include-finer is only valid for 'build' command"; exit 1; }
            # absolute: CMake takes the folder from the build directory
            FINER_PATH=$(cd "${1#*=}" && pwd) || { error "No folder ${1#*=}"; exit 1; }
            ;;
        --include-morph=*)
            [[ "$COMMAND" == "build" ]] || { error " --include-morph is only valid for 'build' command"; exit 1; }
            MORPH_PATH="${1#*=}"
            ;;
        --update-deps)
            [[ "$COMMAND" == "build" ]] || { error " --update-deps is only valid for 'build' command"; exit 1; }
            UPDATE_DEPS="true"
            ;;
        --compilers=*)
            [[ "$COMMAND" == "build" ]] || { error " --compilers is only valid for 'build' command"; exit 1; }
            if [[ ! "$1" =~ ^--compilers=(intel|gnu)$ ]]; then
                error "Invalid value for --compilers. Valid values are 'intel' or 'gnu'."
                exit 1
            fi
            COMPILERS="${1#*=}"
            ;;
        --use-openmp)
            [[ "$COMMAND" == "build" ]] || { error " --use-openmp is only valid for 'build' command"; exit 1; }
            USE_OPENMP="true"
            ;;
        --use-mpi)
            [[ "$COMMAND" == "build" ]] || { error " --use-mpi is only valid for 'build' command"; exit 1; }
            USE_MPI="true"
            ;;
        --use-tecio)
            [[ "$COMMAND" == "build" ]] || { error " --use-tecio is only valid for 'build' command"; exit 1; }
            USE_TECIO="true"
            ;;
        --no-ipo)
            [[ "$COMMAND" == "build" ]] || { error " --no-ipo is only valid for 'build' command"; exit 1; }
            USE_IPO="false"
            ;;
        --preset=*)
            [[ "$COMMAND" == "compile" ]] || { error " --preset is only valid for 'compile' command"; exit 1; }
            PRESET="${1#*=}"
            ;;
        *)
            eval "opts=(\"\${CMD_OPTIONS_${COMMAND}[@]}\")"
            error "Unknown option '$1' for command '$COMMAND'. Valid options: ${opts[@]}"
            usage
            exit 1
            ;;
    esac
    shift
done


# Execute the selected command
case "$COMMAND" in
    build)
        task "Building $project"

        task "Cloning submodules"
        # sync first: a clone made before FiNeR moved (lib/third_party/FiNeR to lib/FiNeR) lacks its new entry
        git submodule sync --quiet
        [[ -z $FINER_PATH ]] && git submodule update --init lib/FiNeR
        [[ $MORPH_PATH == $(pwd)'/lib/MORPH/' ]] && git submodule update --init lib/MORPH
        # ORION (a release, in lib/) and FiNeR's own dependencies are cloned by CMake at configure time.

        task "Configuring and building $project"
        if [[ $COMPILERS == "intel" ]]; then
          export FC="ifx"
          export CC="icx"
          export CXX="icpx"
        elif [[ $COMPILERS == "gnu" ]]; then
          export FC="gfortran"
          export CC="gcc"
          export CXX="g++"
        fi
        log "Build dir: $BUILD_DIR"
        log "Build type: $BUILD_TYPE"
        log "ORION: ${ORION_PATH:-release}"
        log "Update ORION: $UPDATE_DEPS"
        log "FiNeR: ${FINER_PATH:-lib/FiNeR}"
        log "MORPH path: $MORPH_PATH"
        log "Use OpenMP: $USE_OPENMP"
        log "Use MPI: $USE_MPI"
        log "Use TecIO: $USE_TECIO"
        log "Use IPO: $USE_IPO"
        if [[ -z "${FC+x}" || -z "${CXX+x}" || -z "${CC+x}" ]]; then
          log "Compilers not set. CMake will decide."
        else
          log "Compilers: FC=$FC, CXX=$CXX, CC=$CC"
        fi
        CMAKE_ARGS=(-DMORPH_PATH=$MORPH_PATH -DUSE_TECIO=$USE_TECIO -DUSE_OPENMP=$USE_OPENMP -DUSE_MPI=$USE_MPI -DUSE_IPO=$USE_IPO -DCMAKE_BUILD_TYPE=$BUILD_TYPE)
        if [[ -n $ORION_PATH ]]; then
          CMAKE_ARGS+=(-DFETCHCONTENT_SOURCE_DIR_ORION="$ORION_PATH")
        fi
        if [[ -n $FINER_PATH ]]; then
          CMAKE_ARGS+=(-DFETCHCONTENT_SOURCE_DIR_FINER="$FINER_PATH")
        fi
        if [[ $UPDATE_DEPS == "true" ]]; then
          CMAKE_ARGS+=(-DSEMVER_UPDATE=ALL)
        fi
        rm -rf $BUILD_DIR
        cmake -B $BUILD_DIR "${CMAKE_ARGS[@]}" || exit 1
        cmake --build $BUILD_DIR || exit 1
        log "[OK] Compilation successful"

        task "Write CMakePresets.json"
        write_presets
        log "[OK] CMakePresets.json created"
        ;;
    compile)
        task "Compiling $project using CMakePresets (preset ${PRESET:-default})"
        # each preset its own tree: default -> build/, mpi -> build-mpi/
        PDIR=$BUILD_DIR
        [[ "${PRESET:-default}" == "mpi" ]] && PDIR="$DIR/build-mpi"
        cmake --preset ${PRESET:-default} || exit 1
        cmake --build $PDIR || exit 1
        log "[OK] Compilation successful"
        ;;
    *)
        error "Unknown command '$COMMAND'"
        usage
        ;;
esac
