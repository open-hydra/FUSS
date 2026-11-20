# CMake project file for 'FUSS'

cmake_minimum_required(VERSION 3.23)

# Name the code
set(solver FUSS)

project (${solver} LANGUAGES Fortran CXX C)

# -------- Libraries and the releases FUSS is tested with (docs/development/versioning.md) --------
# - ORION. Required by a range of releases with semver_dependency of cmake/SemVerDeps.cmake. Standalone, the release is
# checked out in lib/ORION, out of the build folder, and kept in deps.lock (both git ignored); -DSEMVER_UPDATE=ALL
# selects the newest one again.
set(FUSS_ORION_VERSION 2.1)
set(FUSS_ORION_URL https://github.com/MarcoGrossi92/ORION.git)
# - FiNeR. Directly pinned at one release, the submodule lib/FiNeR (bump the gitlink and the tag together).
set(FUSS_FINER_TAG v2.3.0)
# FiNeR's own dependencies (cmake/FiNeRDependencies.cmake). FREE fetches them as FiNeR's author does
# (`fobis fetch`: the default branch of each); BACKUP checks out the releases below, which FUSS builds and passes its
# tests with. When an upstream change breaks FUSS, set BACKUP here and commit it, or configure one build with
# -DFUSS_FINER_DEPS=BACKUP (kept in its cache until -UFUSS_FINER_DEPS).
if (NOT DEFINED FUSS_FINER_DEPS)
  set(FUSS_FINER_DEPS "FREE")
endif()
set(FUSS_FINER_BACKUP_PINS
  "PENF=v2.0.1"
  "StringiFor=v1.8.0"
  "FACE=v1.1.17"
  "FLAP=v2.5.3"
  "BeFoR64=v1.2.4")

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

  # Libraries building
  # An external clone given by a relative path is taken from the source folder
  foreach(_name ORION FINER)
    if (FETCHCONTENT_SOURCE_DIR_${_name})
      get_filename_component(FETCHCONTENT_SOURCE_DIR_${_name} "${FETCHCONTENT_SOURCE_DIR_${_name}}" ABSOLUTE
                             BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
    endif()
  endforeach()
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/SemVerDeps.cmake)
  semver_dependency(ORION GIT_REPOSITORY ${FUSS_ORION_URL} VERSION ${FUSS_ORION_VERSION} SOURCE_DIR lib/ORION)
  # ORION and FiNeR with its dependencies, the checks of both, and FetchContent_MakeAvailable
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/AddLibraries.cmake)

else()

  # Part of a larger project (hydra), which gives FUSS its libraries. When that project selected ORION with SemVerDeps,
  # its release is checked against FUSS's requirement; nothing is added here.
  get_property(_fuss_selected GLOBAL PROPERTY SEMVER_ORION_SOURCE SET)
  if (_fuss_selected)
    include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/SemVerDeps.cmake)
    semver_dependency(ORION GIT_REPOSITORY ${FUSS_ORION_URL} VERSION ${FUSS_ORION_VERSION})
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
