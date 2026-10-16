# Turns on either OpenMP and/or MPI

IF (USE_OPENMP)
    # Find OpenMP
    FIND_PACKAGE (OpenMP REQUIRED)
    message(STATUS "OpenMP enabled")
ENDIF()

IF (USE_MPI)
    # Find MPI
    IF (NOT MPI_Fortran_FOUND)
        # Fortran only: FUSS links MPI::MPI_Fortran, and the mpi preset bypasses the
        # compiler wrappers for Fortran alone (C and C++ would want their own variables)
        FIND_PACKAGE (MPI REQUIRED COMPONENTS Fortran)
        message(STATUS "MPI enabled, using compilers: ${MPI_Fortran_COMPILER}")
    ENDIF (NOT MPI_Fortran_FOUND)
    # The MPI modules give the choice-buffer routines (MPI_REDUCE, MPI_BCAST, ...) no explicit
    # interface, so gfortran >= 10 rejects one routine called with a scalar and with an array
    # in the same file. Intel's own mpif90 adds this flag for gfortran; the macro adds it only
    # where the compiler knows it (plan 11, 2026-10-16)
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                     Fortran "-fallow-argument-mismatch")
ENDIF()

IF (NOT USE_OPENMP AND NOT USE_MPI)
    # Turn off both OpenMP and MPI
    SET (OMP_NUM_PROCS 0 CACHE
         STRING "Number of processors OpenMP may use" FORCE)
    UNSET (OpenMP_Fortran_FLAGS CACHE)
    UNSET (GOMP_Fortran_LINK_FLAGS CACHE)
    UNSET (MPI_FOUND CACHE)
    UNSET (MPI_COMPILER CACHE)
    UNSET (MPI_LIBRARY CACHE)
ENDIF ()
