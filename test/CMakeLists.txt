# CTest for the FUSS test cases: each runs in a copy under the build tree (run_case.sh).
# Labels:
#   fast        under a minute each, run + verify.py      (ctest -L fast)
#   validation  every case with a verify.py
#   smoke       cases without a verifier, 300 iterations: the run must end cleanly
#   needs-mpi   2 ranks must give the same output files as 1 process

set(_RUN "${CMAKE_CURRENT_SOURCE_DIR}/run_case.sh")
set(_WORK "${CMAKE_CURRENT_BINARY_DIR}/runs")

function(fuss_add_test NAME DIR TIMEOUT LABELS)
  add_test(NAME ${NAME}
           COMMAND bash ${_RUN} ${CMAKE_CURRENT_SOURCE_DIR}/${DIR} ${_WORK}/${NAME} $<TARGET_FILE:${solver}> ${ARGN})
  set_tests_properties(${NAME} PROPERTIES LABELS "${LABELS}" TIMEOUT ${TIMEOUT})
endfunction()

# Verification cases (analytic or reference solutions)
fuss_add_test(QvolPlate          steady-state/qvol_Plate          300 "fast;validation" --verify)
fuss_add_test(TemperaturePlate   steady-state/temperature_Plate   300 "fast;validation" --verify)
fuss_add_test(RadiativeRod       steady-state/radiative_Rod       300 "fast;validation" --verify)
fuss_add_test(ConvectiveQradRod  steady-state/convective_qrad_Rod 300 "fast;validation" --verify)
fuss_add_test(SkewedPlate        steady-state/skewed_Plate        300 "fast;validation" --verify)
fuss_add_test(MultimatPlate      transient/multimat_Plate         300 "fast;validation" --verify)
fuss_add_test(VariableCpPlate    transient/variable_cp_Plate      300 "fast;validation" --verify)
fuss_add_test(OscillatingRod     transient/oscillating_Rod        300 "fast;validation" --verify)
fuss_add_test(ConvectivePlate    steady-state/convective_Plate   1200 "validation"      --verify)

# Smoke runs of the cases without a verifier
foreach(_c nominal irs mg1 mg2 mg3 mg_irs fas_v fas_w fas_w_irs)
  fuss_add_test(ResidualConvergence_${_c} numerics-features/residual_convergence/${_c} 300 "smoke" --iters 300)
endforeach()
fuss_add_test(Chimera          numerics-features/chimera                300 "smoke" --iters 300)
fuss_add_test(ChimeraMG        numerics-features/chimera_mg             300 "smoke" --iters 300)
fuss_add_test(AnisotropicPlate numerics-features/anisotropic_Plate      300 "smoke" --iters 300)
fuss_add_test(TimeConvergence  numerics-features/time_convergence/simdir 300 "smoke" --iters 100)

# The FAS keys must be inert outside a cycle: cycle = none with every FAS key set
# must reproduce the plain multigrid run byte for byte
fuss_add_test(MultigridNoneInert numerics-features/residual_convergence/mg3 300 "smoke" --iters 300 --none-inert)

# dt-face-k: the multimaterial liner runs on the face-conductivity time step (with the cell k it
# dies within a few iterations), and on one material at constant k the key must not change a byte
# (walls 301-305 on the skewed plate, chimera faces on three multigrid levels)
fuss_add_test(MultimatLiner        numerics-features/multimat_liner 300 "smoke" --iters 300)
fuss_add_test(DtFaceKInertWalls    steady-state/skewed_Plate        300 "smoke" --dt-face-k-inert)
fuss_add_test(DtFaceKInertChimera  numerics-features/chimera_mg     300 "smoke" --iters 300 --dt-face-k-inert)

# MPI: 1 process vs 2 ranks
if (USE_MPI)
  fuss_add_test(ChimeraMPI   numerics-features/chimera                      300 "needs-mpi" --iters 300 --mpi 2)
  fuss_add_test(ChimeraMGMPI numerics-features/chimera_mg                   300 "needs-mpi" --iters 300 --mpi 2)
  fuss_add_test(MultigridMPI numerics-features/residual_convergence/mg_irs  300 "needs-mpi" --iters 300 --mpi 2)
  fuss_add_test(FASMPI       numerics-features/residual_convergence/fas_w_irs 300 "needs-mpi" --iters 300 --mpi 2)
  fuss_add_test(MultimatLinerMPI numerics-features/multimat_liner          300 "needs-mpi" --iters 300 --mpi 2)
endif()
