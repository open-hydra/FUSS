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
  # Locate and set parallelization libraries. There are some CMake peculiarities
  # taken care of here, such as the fact that the FindOpenMP routine doesn't know
  # about Fortran.
  INCLUDE(${CMAKE_MODULE_PATH}/SetParallelizationLibrary.cmake)

  # Set directories
  include(GNUInstallDirs)
  set (CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")

  # ORION
  if (NOT TARGET ORION)
    set ( orionD "${ORION_PATH}" "${CMAKE_CURRENT_SOURCE_DIR}/build/lib/ORION/")
    add_subdirectory(${orionD})
  endif()

  # FiNeR
  if (NOT TARGET FiNeR)
    set ( finerD "${FINER_PATH}" "${CMAKE_CURRENT_SOURCE_DIR}/build/lib/third_party/FiNeR/")
    add_subdirectory(${finerD})
  endif()

endif()

# Build library
add_subdirectory(src/lib)

# Build application (only if this file has not been called by another CMakeLists)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  add_subdirectory(src/app)
endif()
