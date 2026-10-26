file(GLOB_RECURSE FUSS_SRC
    CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/*.f90
)

add_library(FUSSL ${FUSS_SRC})
add_library(FUSSL::FUSSL ALIAS FUSSL)
# MORPHL: the mesh-motion component (lib/MORPH), called through numerics/space/Adapter_Morph.f90 only
target_link_libraries(FUSSL PRIVATE FiNeR::FiNeR ORION MORPHL)

target_include_directories(FUSSL INTERFACE $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)

IF(USE_OPENMP)
    SET_TARGET_PROPERTIES(FUSSL PROPERTIES
                          COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
                          LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
ENDIF(USE_OPENMP)

IF(USE_MPI)
    target_compile_definitions(FUSSL PUBLIC USE_MPI)
    target_link_libraries(FUSSL PUBLIC MPI::MPI_Fortran)
ENDIF(USE_MPI)

# The TecIO switch is ORION's (USE_TECIO). Without TecIO ORION stops with exit
# status 0 at the first binary Tecplot read or write, so FUSS refuses the input.
IF(DEFINED USE_TECIO AND NOT USE_TECIO)
    target_compile_definitions(FUSSL PUBLIC FUSS_NO_TECIO)
ENDIF()


# # Install library
# install(TARGETS FUSSL
#     EXPORT FUSSLTargets
#     ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
#     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
#     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
# )

# # Install Fortran modules
# install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
#         DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# # Install export set
# install(EXPORT FUSSLTargets
#     FILE FUSSLTargets.cmake
#     NAMESPACE FUSS::
#     DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/FUSS
# )

# write_basic_package_version_file(
#     "${CMAKE_CURRENT_BINARY_DIR}/FUSSConfigVersion.cmake"
#     VERSION 1.0.0
#     COMPATIBILITY AnyNewerVersion
# )

# configure_package_config_file(
#     "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FUSSConfig.cmake.in"
#     "${CMAKE_CURRENT_BINARY_DIR}/FUSSConfig.cmake"
#     INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/FUSS
# )

# install(FILES
#     "${CMAKE_CURRENT_BINARY_DIR}/FUSSConfig.cmake"
#     "${CMAKE_CURRENT_BINARY_DIR}/FUSSConfigVersion.cmake"
#     DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/FUSS
# )