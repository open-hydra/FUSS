# CMake project file for 'FUSS'

cmake_minimum_required(VERSION 3.23)

# Name the code
set(solver FUSS)

project (${solver} LANGUAGES Fortran CXX C)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)

  set(CMAKE_VERBOSE_MAKEFILE FALSE)
  # seach path for cmake modules
  set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

  # Set some options the user may choose
  # Uncomment the below if you want the user to choose a parallelization library
  OPTION(USE_MPI "Use the MPI library for parallelization" OFF)
  OPTION(USE_OPENMP "Use OpenMP for parallelization" OFF)
  # Interprocedural (link-time) optimization. On by default because it is worth ~10% -- see
  # docs/development/performance_notes.md -- but switchable off for toolchains where it misbehaves. It also self-disables
  # if unsupported.
  OPTION(USE_IPO "Use interprocedural optimization in RELEASE builds" ON)
  # Lets the check on ifx 2024.1 below pass without IPO (to test the compiler, see docs/getting-started/installation.md).
  OPTION(FUSS_ALLOW_IFX_NOIPO "Do not refuse ifx 2024.1 at -O2/-O3 without IPO" OFF)

  # ------------------------------------------------------------------------
  # LTO-aware archiver (GNU)
  #
  # SetFortranFlags.cmake adds "-flto" to the RELEASE flags for GNU. The object
  # files then carry LTO bytecode instead of ordinary symbols, and plain
  # binutils ar/ranlib cannot read them ("plugin needed to handle lto object").
  # The result is a libFUSSL.a with an EMPTY symbol index and ~200 undefined
  # __fuss_*_MOD_* references when linking bin/FUSS -- after the build has
  # already reported 97% progress. gcc-ar/gcc-ranlib/gcc-nm are thin wrappers
  # that load the LTO plugin and index the archive correctly.
  #
  # Two reasons this lives HERE and not in SetFortranFlags.cmake:
  #   1. That file returns early when the flags are already cached, so on an
  #      incremental re-configure (install.sh compile / cmake --preset) the
  #      block would be skipped and CMAKE_AR would silently revert.
  #   2. project() loads CMakeFortranCompiler.cmake, which does a plain
  #      "set(CMAKE_AR /usr/bin/ar)". A normal variable shadows the cache entry
  #      of the same name, so BOTH must be set -- cache for persistence, normal
  #      to win in this scope and to propagate into add_subdirectory() for the
  #      ORION / FiNeR sub-projects.
  # gcc-ar is a superset of ar, so this is applied for GNU regardless of build
  # type rather than being made conditional on the flag string.
  # ------------------------------------------------------------------------
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_Fortran_COMPILER_AR AND CMAKE_Fortran_COMPILER_RANLIB)
      # version-matched wrappers CMake already detected (e.g. gcc-ar-12)
      set(LTO_AR     "${CMAKE_Fortran_COMPILER_AR}")
      set(LTO_RANLIB "${CMAKE_Fortran_COMPILER_RANLIB}")
    else()
      find_program(LTO_AR     NAMES gcc-ar     DOC "LTO-aware ar")
      find_program(LTO_RANLIB NAMES gcc-ranlib DOC "LTO-aware ranlib")
    endif()
    find_program(LTO_NM NAMES gcc-nm DOC "LTO-aware nm")

    if(LTO_AR AND LTO_RANLIB)
      set(CMAKE_AR     "${LTO_AR}"     CACHE FILEPATH "LTO-aware archiver" FORCE)
      set(CMAKE_RANLIB "${LTO_RANLIB}" CACHE FILEPATH "LTO-aware ranlib"   FORCE)
      set(CMAKE_AR     "${LTO_AR}")
      set(CMAKE_RANLIB "${LTO_RANLIB}")
      if(LTO_NM)
        set(CMAKE_NM "${LTO_NM}" CACHE FILEPATH "LTO-aware nm" FORCE)
        set(CMAKE_NM "${LTO_NM}")
      endif()
      message(STATUS "GNU toolchain: archiving with ${LTO_AR} / ${LTO_RANLIB}")
    else()
      message(WARNING
        "No LTO-aware ar/ranlib found (gcc-ar / gcc-ranlib). If the RELEASE "
        "build enables -flto, the static library will get an empty symbol "
        "index and linking bin/FUSS will fail with undefined references.")
    endif()
  endif()

  # This INCLUDE statement executes code that sets the compile flags for DEBUG,
  # RELEASE, and TESTING. You should  review this file and make sure the flags 
  # are to your liking.
  INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake) 

  # ifx 2024.1 compiles a loop wrongly at -O2 and -O3 without interprocedural optimization: in the ghost metrics of MOSE (the
  # sister code), a minimal program in standard Fortran (an automatic array of a derived type, bounds from a module variable,
  # filled element by element and passed to a procedure) gives wrong values at -O2 and -O3, also at the default level with no
  # -O option, and the right ones at -O0, -O1 and with -ipo. FUSS built so (the RELEASE build of ifx before the
  # port of this check: its -ipo was dropped by the probe) passes the stress battery, so no miscompiled loop is known in FUSS; the
  # build is refused all the same, as in MOSE, because the pattern is ordinary Fortran. The level is the last -O option of the
  # flags of the build type (-O2, the default of ifx, when there is none); IPO is -ipo in those flags (the RELEASE build adds
  # it with USE_IPO, the default, and drops it when the toolchain cannot do it) or CMAKE_INTERPROCEDURAL_OPTIMIZATION. This
  # only reads the flags: it changes none. -DFUSS_ALLOW_IFX_NOIPO=ON turns the refusal into a warning.
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM" AND CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 2024.1
     AND CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 2024.2)
    string(TOUPPER "${CMAKE_BUILD_TYPE}" _ifx_bt)
    separate_arguments(_ifx_flags UNIX_COMMAND "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${_ifx_bt}}")
    set(_ifx_level "no -O option (-O2)")
    set(_ifx_high TRUE)
    set(_ifx_ipo FALSE)
    if(CMAKE_INTERPROCEDURAL_OPTIMIZATION OR CMAKE_INTERPROCEDURAL_OPTIMIZATION_${_ifx_bt})
      set(_ifx_ipo TRUE)
    endif()
    foreach(_ifx_f IN LISTS _ifx_flags)
      if(_ifx_f MATCHES "^[-/]O[01]$")
        set(_ifx_level "${_ifx_f}")
        set(_ifx_high FALSE)
      elseif(_ifx_f MATCHES "^[-/]O([23s]|fast)?$")
        set(_ifx_level "${_ifx_f}")
        set(_ifx_high TRUE)
      elseif(_ifx_f MATCHES "^[-/]Q?ipo$")
        set(_ifx_ipo TRUE)
      endif()
    endforeach()
    if(_ifx_high AND NOT _ifx_ipo)
      set(_ifx_msg "ifx ${CMAKE_Fortran_COMPILER_VERSION} with ${_ifx_level} and without IPO compiles a loop of the ghost "
                   "metrics of MOSE wrongly (an error of the compiler: a minimal program in standard Fortran gives wrong "
                   "values at -O2 and -O3 and the right ones at -O0, -O1 and with -ipo); no such loop is known in FUSS, "
                   "but the build is refused as in MOSE. Build with IPO (RELEASE with -DUSE_IPO=ON, the default), with -O0 "
                   "or -O1 (DEBUG), or with gfortran (-DFUSS_ALLOW_IFX_NOIPO=ON accepts this build). See "
                   "docs/getting-started/installation.md, \"ifx 2024.1 and IPO\".")
      if(FUSS_ALLOW_IFX_NOIPO)
        message(WARNING ${_ifx_msg})
      else()
        message(FATAL_ERROR ${_ifx_msg})
      endif()
    endif()
  endif()
  # Locate and set parallelization libraries. There are some CMake peculiarities
  # taken care of here, such as the fact that the FindOpenMP routine doesn't know
  # about Fortran.
  INCLUDE(${CMAKE_MODULE_PATH}/SetParallelizationLibrary.cmake)

  # Set directories
  include(GNUInstallDirs)
  set (CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")

  # -------- Path fallbacks (standalone builds) --------
  if (NOT DEFINED FINER_PATH)
    set(FINER_PATH "${CMAKE_CURRENT_SOURCE_DIR}/lib/third_party/FiNeR/")
  endif()
  if (NOT DEFINED FINER_DEPS_PATH)
    set(FINER_DEPS_PATH "${CMAKE_CURRENT_SOURCE_DIR}/lib/third_party/")
  endif()

  # -------- Helper: initialize a git submodule from the correct git root --------
  macro(ensure_submodule_initialized SUBMOD_PATH)
    if (NOT EXISTS "${SUBMOD_PATH}/CMakeLists.txt")
      message(STATUS "Submodule at ${SUBMOD_PATH} not initialized. Attempting git submodule update --init...")
      get_filename_component(_sm_parent "${SUBMOD_PATH}" DIRECTORY)
      execute_process(
        COMMAND git rev-parse --show-toplevel
        WORKING_DIRECTORY "${_sm_parent}"
        OUTPUT_VARIABLE _sm_git_root
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE _sm_git_result)
      if(_sm_git_result EQUAL 0)
        file(RELATIVE_PATH _sm_rel "${_sm_git_root}" "${SUBMOD_PATH}")
        execute_process(
          COMMAND git submodule update --init "${_sm_rel}"
          WORKING_DIRECTORY "${_sm_git_root}"
          RESULT_VARIABLE _sm_init_result)
        if(NOT _sm_init_result EQUAL 0)
          message(FATAL_ERROR "Failed to initialize submodule: ${SUBMOD_PATH}")
        endif()
      else()
        message(FATAL_ERROR "Could not find git root for submodule: ${SUBMOD_PATH}")
      endif()
    endif()
  endmacro()

  # ORION: the newest release that satisfies the requirement below, checked out in lib/ORION (cmake/SemVerDeps.cmake). The
  # release stays that of the lock file deps.lock (not tracked), -DSEMVER_UPDATE=ALL selects the newest one again,
  # -DFETCHCONTENT_SOURCE_DIR_ORION=<path> builds a local checkout instead. A project that adds FUSS (hydra) chooses the
  # release of its build tree
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/SemVerDeps.cmake)
  semver_dependency(ORION GIT_REPOSITORY https://github.com/MarcoGrossi92/ORION.git VERSION 2.1 SOURCE_DIR lib/ORION)
  FetchContent_MakeAvailable(ORION)

  # Libraries of FiNeR: FiNeR 2.0.7 and later do not carry them as submodules (its CMake build expects them in its
  # src/third_party/, fetched with FoBiS) and add each one only when the target is not defined yet: FUSS provides
  # them, in the order of their dependencies (FLAP needs PENF and FACE, StringiFor needs PENF and BeFoR64)
  foreach (finer_dep PENF FACE FLAP BeFoR64 StringiFor)
    if (NOT TARGET ${finer_dep})
      ensure_submodule_initialized("${FINER_DEPS_PATH}/${finer_dep}")
      add_subdirectory("${FINER_DEPS_PATH}/${finer_dep}" "${CMAKE_BINARY_DIR}/lib/third_party/${finer_dep}")
    endif()
  endforeach()

  # FiNeR
  if (NOT TARGET FiNeR)
    ensure_submodule_initialized("${FINER_PATH}")
    # FUSS needs FiNeR v2.1.0 or later: before it, a comment line that holds "=" or a bracket was taken for an option or
    # a section header, and the reading of its section of input.ini ended there. A clone that pulls without
    # "git submodule update" keeps the FiNeR it had checked out: ensure_submodule_initialized
    # leaves it as it is, and it would build in silence. FiNeR names its version in a VERSION file ("vX.Y.Z") since
    # v2.0.7; a FiNeR without it is older
    set(finer_version "")
    if (EXISTS "${FINER_PATH}/VERSION")
      file(STRINGS "${FINER_PATH}/VERSION" finer_version LIMIT_COUNT 1)
      string(REGEX REPLACE "^[ \t]*v?([0-9][0-9.]*).*$" "\\1" finer_version "${finer_version}")
    endif()
    if (NOT "${finer_version}" VERSION_GREATER_EQUAL "2.1.0")
      if ("${finer_version}" STREQUAL "")
        set(finer_found "a FiNeR without a VERSION file, older than v2.0.7")
      else()
        set(finer_found "FiNeR v${finer_version}")
      endif()
      message(FATAL_ERROR "FiNeR at ${FINER_PATH} is ${finer_found}: FUSS needs FiNeR v2.1.0 or later (it pins v2.3.0), "
                          "as before v2.1.0 a comment line of input.ini that holds \"=\" ends the reading of its section. "
                          "Run: git submodule update --init")
    endif()
    set(finerD "${FINER_PATH}" "${CMAKE_BINARY_DIR}/lib/third_party/FiNeR")
    add_subdirectory(${finerD})
  endif()

endif()

# MORPH, the mesh-motion component (plan 09): the git submodule lib/MORPH, or the tree
# MORPH_PATH names (install.sh build --include-morph=). Only its library MORPHL is
# configured, with this project's flags. A parent that already provides MORPHL keeps its
# own; one that does not (hydra embeds FUSS) gets FUSS's submodule.
if (NOT TARGET MORPHL)
  set ( MORPH_PATH "${CMAKE_CURRENT_SOURCE_DIR}/lib/MORPH" CACHE PATH "MORPH source tree (default: the submodule lib/MORPH)")
  if (NOT EXISTS "${MORPH_PATH}/CMakeLists.txt")
    message(FATAL_ERROR "MORPH not found at ${MORPH_PATH}: run ./install.sh update (git submodule update --init lib/MORPH), or pass -DMORPH_PATH=<MORPH tree>")
  endif()
  set ( morphD "${MORPH_PATH}" "${CMAKE_CURRENT_BINARY_DIR}/lib/MORPH")
  add_subdirectory(${morphD})
endif()

# Build library
add_subdirectory(src/lib)

# Build application and tests (only if this file has not been called by another CMakeLists)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  add_subdirectory(src/app)
  enable_testing()
  add_subdirectory(test)
endif()
