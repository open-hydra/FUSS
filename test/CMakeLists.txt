# CTest for the FUSS test cases: each runs in a copy under the build tree (run_case.sh).
# Labels:
#   fast        under a minute each, run + verify.py      (ctest -L fast)
#   validation  every case with a verify.py
#   smoke       cases without a verifier, 300 iterations: the run must end cleanly
#   needs-mpi   2 ranks must give the same output files as 1 process; ranks x threads layouts the same bits

set(_RUN "${CMAKE_CURRENT_SOURCE_DIR}/run_case.sh")
set(_WORK "${CMAKE_CURRENT_BINARY_DIR}/runs")

function(fuss_add_test NAME DIR TIMEOUT LABELS)
  add_test(NAME ${NAME}
           COMMAND bash ${_RUN} ${CMAKE_CURRENT_SOURCE_DIR}/${DIR} ${_WORK}/${NAME} $<TARGET_FILE:${solver}> ${ARGN})
  set_tests_properties(${NAME} PROPERTIES LABELS "${LABELS}" TIMEOUT ${TIMEOUT})
endfunction()

# Scripts of input-io/ (input checks, stops, restart, decomposition map): fuss_add_script(NAME SCRIPT TIMEOUT LABELS [args]);
# each runs `bash SCRIPT <FUSS binary> <work folder> [args]` on copies of the test cases and prints one CASE line per check
function(fuss_add_script NAME SCRIPT TIMEOUT LABELS)
  add_test(NAME ${NAME}
           COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/input-io/${SCRIPT} $<TARGET_FILE:${solver}> ${_WORK}/${NAME} ${ARGN})
  set_tests_properties(${NAME} PROPERTIES LABELS "${LABELS}" TIMEOUT ${TIMEOUT})
endfunction()

# Verification cases (analytic or reference solutions)
fuss_add_test(QvolPlate          steady-state/qvol_Plate          300 "fast;validation" --verify)
fuss_add_test(TemperaturePlate   steady-state/temperature_Plate   300 "fast;validation" --verify)
fuss_add_test(RadiativeRod       steady-state/radiative_Rod       300 "fast;validation" --verify)
fuss_add_test(ConvectiveQradRod  steady-state/convective_qrad_Rod 300 "fast;validation" --verify)
fuss_add_test(SkewedPlate        steady-state/skewed_Plate        300 "fast;validation" --verify)
fuss_add_test(MultilayerSlab     steady-state/multilayer_Slab     300 "fast;validation" --verify)
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

# Input checks and stops: the values of input.ini are checked and every error ends the run with exit code 1; input.ini read
# again during the run (ini-error); the errors of bc.txt, time tables, probes, the initial condition and the registry
fuss_add_script(InputChecks  input-checks.sh  300 "fast")
fuss_add_script(RuntimeInput runtime-input.sh 300 "fast")
fuss_add_script(ExitCodes    exit-codes.sh    300 "fast")

# decomposition-map: the rank of every block from a file of INPUT/ (one rank here; MPIDecompositionMap with MPI)
fuss_add_script(DecompositionMap decomposition-map.sh 300 "fast")

# Restart: the run split in two writes the files of the run in one go, byte for byte, in every format (the counters stored
# with the solution); the .szplt case needs TecIO
fuss_add_script(Restart       restart.sh 600 "fast" core)
if (USE_TECIO)
  fuss_add_script(RestartBinary restart.sh 300 "fast" szplt)
endif()

# Material jumps inside a block: the 4-layer liner at vnn 1.0 (series face conductivities and the
# face-conductivity time step)
fuss_add_test(MultimatLiner    numerics-features/multimat_liner     300 "smoke" --iters 300)

# Chimera edges: a linear field on a sheared chimera pair (two chimera faces meeting at an edge in each block), walls
# consistent with it, is preserved to round-off (the ghost row beyond an edge is extrapolated, not copied)
add_test(NAME ChimeraEdge
         COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/mpi/bit-identity/chimera_edge.sh $<TARGET_FILE:${solver}> ${_WORK}/ChimeraEdge)
set_tests_properties(ChimeraEdge PROPERTIES LABELS "fast;validation" TIMEOUT 300)

# MPI: 1 process vs 2 ranks
if (USE_MPI)
  fuss_add_test(ChimeraMPI   numerics-features/chimera                      300 "needs-mpi" --iters 300 --mpi 2)
  fuss_add_test(ChimeraMGMPI numerics-features/chimera_mg                   300 "needs-mpi" --iters 300 --mpi 2)
  fuss_add_test(MultigridMPI numerics-features/residual_convergence/mg_irs  300 "needs-mpi" --iters 300 --mpi 2)
  fuss_add_test(FASMPI       numerics-features/residual_convergence/fas_w_irs 300 "needs-mpi" --iters 300 --mpi 2)
  fuss_add_test(MultimatLinerMPI numerics-features/multimat_liner          300 "needs-mpi" --iters 300 --mpi 2)

  # The same output bits on 1x1, 1x4, 2x1, 4x1 and 2x2 (ranks x threads) on a box of 4 blocks and a chimera pair
  add_test(NAME MPIBitIdentity
           COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/mpi/bit-identity/run.sh $<TARGET_FILE:${solver}> ${_WORK}/MPIBitIdentity)
  set_tests_properties(MPIBitIdentity PROPERTIES LABELS "needs-mpi" TIMEOUT 900)
  fuss_add_script(MPIRuntimeInput mpi-runtime-input.sh 300 "needs-mpi")
  fuss_add_script(MPIExitCodes    mpi-exit-codes.sh    300 "needs-mpi")
  fuss_add_script(MPIRestart      mpi-restart.sh       300 "needs-mpi")
  fuss_add_script(MPIDecompositionMap mpi-decomposition-map.sh 300 "needs-mpi")
endif()
