function(add_fortran_executable target source)
  
  add_executable(${target} ${source})
  # Link with the Fortran driver: TecIO brings C++ into the link, and a C++ driver
  # (icpx) leaves out the Fortran runtime object that provides main().
  set_target_properties(${target} PROPERTIES LINKER_LANGUAGE Fortran)
  target_link_libraries(${target} PRIVATE FUSSL FiNeR::FiNeR)

  if (USE_TECIO)
    target_link_libraries(${target} PRIVATE -pthread)
  endif()

  target_include_directories(${target} INTERFACE
      $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)

  if (USE_OPENMP)
    target_compile_options(${target} PRIVATE ${OpenMP_Fortran_FLAGS})
    target_link_options(${target} PRIVATE ${OpenMP_Fortran_FLAGS})
  endif()

  if (USE_MPI)
    target_compile_definitions(${target} PRIVATE USE_MPI)
    target_link_libraries(${target} PRIVATE MPI::MPI_Fortran)
  endif()

endfunction()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# Every build tree writes into bin/: the suffix keeps two trees from overwriting one
# executable (the mpi preset: bin/FUSS-mpi; the default tree keeps bin/FUSS).
set(FUSS_EXE_SUFFIX "" CACHE STRING "Suffix of the executables built by this tree (mpi preset: -mpi)")

add_fortran_executable(${solver} "main.f90")
add_fortran_executable(DocGen "docgen.f90")
set_target_properties(${solver} PROPERTIES OUTPUT_NAME ${solver}${FUSS_EXE_SUFFIX})
set_target_properties(DocGen PROPERTIES OUTPUT_NAME DocGen${FUSS_EXE_SUFFIX})
