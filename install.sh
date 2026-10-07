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
    --include-orion=<path>  Set external ORION path
    --include-finer=<path>  Set external FiNeR path
    --compilers=<name>      Set compilers (intel, gnu)
    --use-openmp            Use OpenMP
    --use-mpi               Use MPI
    --use-tecio             Use TecIO

  compile                   Compile the program using the CMakePresets file

  update                    Download git submodules
    --remote                Use the latest remote commit

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
  FC=$(grep '^CMAKE_Fortran_COMPILER:FILEPATH=' "$BUILD_DIR/CMakeCache.txt" | cut -d= -f2-)
  CC=$(grep '^CMAKE_C_COMPILER:FILEPATH=' "$BUILD_DIR/CMakeCache.txt" | cut -d= -f2-)
  CXX=$(grep '^CMAKE_CXX_COMPILER:FILEPATH=' "$BUILD_DIR/CMakeCache.txt" | cut -d= -f2-)

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
        "ORION_PATH": "${ORION_PATH}",
        "FINER_PATH": "${FINER_PATH}",
        "CMAKE_BUILD_TYPE": "${BUILD_TYPE}",
        "CMAKE_Fortran_COMPILER": "${FC}",
        "CMAKE_C_COMPILER": "${CC}",
        "CMAKE_CXX_COMPILER": "${CXX}",
        "USE_TECIO": "${USE_TECIO}",
        "USE_OPENMP": "${USE_OPENMP}",
        "USE_MPI": "${USE_MPI}"
      }
    }
  ]
}
EOF
}


# FiNeR (v2.2.0+) no longer tracks its five dependencies as submodules, yet its CMakeLists.txt still
# add_subdirectory's them from src/third_party/. Put each one there at a pinned commit (the same five
# commits hydra, IGLOO, ICE and ATLAS pin: hydra/scripts/utils/finer_deps.sh). FiNeR's own
# src/third_party/.gitignore ignores the clones, so the submodule stays clean.
FINER_DEPS=(
    "PENF       34e10af852f81822bac0df7c5f73936e0655f54f"
    "StringiFor 77b75b1d7f7012d2984135b231008f751ab32dda"
    "FACE       95226cd6135b9daeb8ff2c748d2065deee230358"
    "FLAP       7f5ec8e1f26be5e2a263fdd2118501b340c0d461"
    "BeFoR64    40d19c5f6127088126606b90448e99e02d8b641c"
)
function fill_finer_deps() {
  local finer="${1%/}" entry dep sha dir url
  [[ -f "$finer/CMakeLists.txt" ]] || { error "FiNeR not found at $finer"; return 1; }
  for entry in "${FINER_DEPS[@]}"; do
    read -r dep sha <<< "$entry"
    dir="$finer/src/third_party/$dep"; url="https://github.com/szaghi/$dep"
    if [[ ! -e "$dir" ]] || [[ -d "$dir" && -z "$(ls -A "$dir")" ]]; then
      git clone -q "$url" "$dir" || return 1
    elif ! git -C "$dir" rev-parse --git-dir > /dev/null 2>&1; then
      error "finer_deps: $dir exists but is not a git checkout; remove it and rerun"; return 1
    fi
    if [[ "$(git -C "$dir" rev-parse HEAD 2>/dev/null)" != "$sha" ]]; then
      git -C "$dir" cat-file -e "$sha^{commit}" 2>/dev/null || git -C "$dir" fetch -q "$url" "$sha" || return 1
      git -C "$dir" checkout -q "$sha" || return 1
    fi
  done
  log "[OK] FiNeR dependencies at their pinned commits"
}

# Default global values
COMMAND=""
COMPILERS=""
ORION_PATH=$(pwd)'/lib/ORION/'
FINER_PATH=$(pwd)'/lib/third_party/FiNeR/'
USE_OPENMP="false"
USE_MPI="false"
USE_TECIO="false"
REMOTE="false"
BUILD_TYPE="RELEASE"

# Define allowed options for each command using regular arrays
CMD=("build" "compile" "update")
CMD_OPTIONS_build=("--compilers --include-orion --include-finer --use-openmp --use-mpi --use-tecio")
CMD_OPTIONS_update=("--remote")

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
            ORION_PATH="${1#*=}"
            ;;
        --include-finer=*)
            [[ "$COMMAND" == "build" ]] || { error " --include-finer is only valid for 'build' command"; exit 1; }
            FINER_PATH="${1#*=}"
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
        --remote)
            [[ "$COMMAND" == "update" ]] || { error " --remote is only valid for 'update' command"; exit 1; }
            REMOTE="true"
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
        [[ $ORION_PATH == $(pwd)'/lib/ORION/' ]] && git submodule update --init lib/ORION
        [[ $FINER_PATH == $(pwd)'/lib/third_party/FiNeR/' ]] && git submodule update --init --recursive lib/third_party/FiNeR
        fill_finer_deps "$FINER_PATH" || exit 1

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
        log "ORION path: $ORION_PATH"
        log "FINER path: $FINER_PATH"
        log "Use OpenMP: $USE_OPENMP"
        log "Use MPI: $USE_MPI"
        log "Use TecIO: $USE_TECIO"
        if [[ -z "${FC+x}" || -z "${CXX+x}" || -z "${CC+x}" ]]; then
          log "Compilers not set. CMake will decide."
        else
          log "Compilers: FC=$FC, CXX=$CXX, CC=$CC"
        fi
        rm -rf $BUILD_DIR
        cmake -B $BUILD_DIR -DORION_PATH=$ORION_PATH -DFINER_PATH=$FINER_PATH -DUSE_TECIO=$USE_TECIO -DUSE_OPENMP=$USE_OPENMP -DUSE_MPI=$USE_MPI -DCMAKE_BUILD_TYPE=$BUILD_TYPE || exit 1
        cmake --build $BUILD_DIR || exit 1
        log "[OK] Compilation successful"

        task "Write CMakePresets.json"
        write_presets
        log "[OK] CMakePresets.json created"
        ;;
    compile)
        task "Compiling $project using CMakePresets"
        cmake --preset default || exit 1
        cmake --build $BUILD_DIR || exit 1
        log "[OK] Compilation successful"
        ;;
    update)
        task "Updating git submodules"
        if [[ "$REMOTE" == "true" ]]; then
          log "Updating submodules to latest remote commit"
          log "NOTE: FiNeR's five dependencies stay at the commits pinned in this script (FINER_DEPS); re-pin them when FiNeR moves"
          git submodule update --init --remote
        else
          log "Updating submodules to current commit"
          git submodule update --init
        fi
        fill_finer_deps "$FINER_PATH" || exit 1
        log "[OK] Submodules updated"
        ;;
    *)
        error "Unknown command '$COMMAND'"
        usage
        ;;
esac
