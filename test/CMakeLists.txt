# CTest for the FUSS test cases: each runs in a copy under the build tree (run_case.sh).
#
# Tiers are expressed through LABELS, as in MOSE (test/TESTING-PLAN.md of MOSE):
#   fast          under a minute each (the tier in about a minute at -j4): the selector of the pre-push check, ctest -L fast
#   validation    comparison with an analytic or reference solution (also the cases of the fast tier that have one)
#   verification  full forms of verification cases (families of grids, long runs): ConvergenceOrder, MultigridSolution; run on
#                 request, ctest -L verification
#   smoke         cases without a verifier, 300 iterations: the run must end cleanly (FUSS only; MOSE has none)
#   needs-mpi     registered only with USE_MPI=ON: 2 ranks, or ranks x threads layouts, must give the same output files
# Second-level tags, as MOSE's: physics (a model or a boundary condition against its solution), numerics (a property of the
# discretisation), implementation (inputs, stops, restart, parallel runs), unit (a program that calls a routine).
#
# The PR tier runs everything except the MPI tier and the full forms of the verification cases:
#   ctest -LE "needs-mpi|verification"
# and, in a build with USE_MPI=ON, the MPI tier:  ctest -L needs-mpi
# The restart test of the Tecplot binary format (.szplt) is registered only with USE_TECIO=ON (ORION's option).

set(_RUN "${CMAKE_CURRENT_SOURCE_DIR}/run_case.sh")
set(_WORK "${CMAKE_CURRENT_BINARY_DIR}/runs")

# fuss_register(NAME TIMEOUT LABELS command...): any test, with its labels and timeout
function(fuss_register NAME TIMEOUT LABELS)
  add_test(NAME ${NAME} COMMAND ${ARGN})
  set_tests_properties(${NAME} PROPERTIES LABELS "${LABELS}" TIMEOUT ${TIMEOUT})
endfunction()

# fuss_add_test(NAME DIR TIMEOUT LABELS [run_case.sh options]): a case of DIR (relative to test/) run by run_case.sh
function(fuss_add_test NAME DIR TIMEOUT LABELS)
  fuss_register(${NAME} ${TIMEOUT} "${LABELS}"
                bash ${_RUN} ${CMAKE_CURRENT_SOURCE_DIR}/${DIR} ${_WORK}/${NAME} $<TARGET_FILE:${solver}> ${ARGN})
endfunction()

# fuss_add_script(NAME SCRIPT TIMEOUT LABELS [args]): `bash SCRIPT <FUSS binary> <work folder> [args]`, SCRIPT relative to test/.
# The scripts of input-io/ and mpi/ run the solver on copies of the test cases and print one CASE line per check
function(fuss_add_script NAME SCRIPT TIMEOUT LABELS)
  fuss_register(${NAME} ${TIMEOUT} "${LABELS}"
                bash ${CMAKE_CURRENT_SOURCE_DIR}/${SCRIPT} $<TARGET_FILE:${solver}> ${_WORK}/${NAME} ${ARGN})
endfunction()

# fuss_add_folder_script(NAME DIR SCRIPT TIMEOUT LABELS): `bash DIR/SCRIPT <DIR> <work folder> <FUSS binary>` (the folder first,
# DIR relative to test/): the scripts of the patch test and of the connection split
function(fuss_add_folder_script NAME DIR SCRIPT TIMEOUT LABELS)
  fuss_register(${NAME} ${TIMEOUT} "${LABELS}"
                bash ${CMAKE_CURRENT_SOURCE_DIR}/${DIR}/${SCRIPT}
                ${CMAKE_CURRENT_SOURCE_DIR}/${DIR} ${_WORK}/${NAME} $<TARGET_FILE:${solver}>)
endfunction()

# fuss_add_mpi_script(NAME SCRIPT TIMEOUT LABELS [args]): as fuss_add_script for a script that launches up to 4 ranks
# itself with `mpirun` (as MOSE's mose_add_mpi_script): CTest keeps 4 processors for it under -j
function(fuss_add_mpi_script NAME SCRIPT TIMEOUT LABELS)
  fuss_add_script(${NAME} ${SCRIPT} ${TIMEOUT} "${LABELS}" ${ARGN})
  set_tests_properties(${NAME} PROPERTIES PROCESSORS 4)
endfunction()

# ===========================================================================
#  Verification cases against analytic or reference solutions (verify.py)
# ===========================================================================
fuss_add_test(QvolPlate          steady-state/qvol_Plate          300 "fast;validation;physics" --verify)
fuss_add_test(TemperaturePlate   steady-state/temperature_Plate   300 "fast;validation;physics" --verify)
fuss_add_test(RadiativeRod       steady-state/radiative_Rod       300 "fast;validation;physics" --verify)
fuss_add_test(ConvectiveQradRod  steady-state/convective_qrad_Rod 300 "fast;validation;physics" --verify)
fuss_add_test(SkewedPlate        steady-state/skewed_Plate        300 "fast;validation;numerics" --verify)
fuss_add_test(MultilayerSlab     steady-state/multilayer_Slab     300 "fast;validation;physics" --verify)
fuss_add_test(MultimatPlate      transient/multimat_Plate         300 "fast;validation;physics" --verify)
fuss_add_test(VariableCpPlate    transient/variable_cp_Plate      300 "fast;validation;physics" --verify)
fuss_add_test(OscillatingRod     transient/oscillating_Rod        300 "fast;validation;physics" --verify)
# about two minutes at one thread: outside the fast tier
fuss_add_test(ConvectivePlate    steady-state/convective_Plate   1200 "validation;physics"      --verify)

# 3-D cells skewed in all three pairs of directions: the time step takes the factor (2 + m + 1/m)/4 of the spectral radius
# of the diffusion operator (verify.py: the step of every cell against a replica, in both passes)
fuss_add_test(Skewed3D           numerics-features/skewed_3D      300 "fast;validation;numerics" --verify)

# Chimera edges: a linear field on a sheared chimera pair (two chimera faces meeting at an edge in each block), walls
# consistent with it, is preserved to round-off (the ghost row beyond an edge is extrapolated, not copied)
fuss_register(ChimeraEdge 300 "fast;validation;numerics"
              bash ${CMAKE_CURRENT_SOURCE_DIR}/mpi/bit-identity/chimera_edge.sh $<TARGET_FILE:${solver}> ${_WORK}/ChimeraEdge)

# Linear-field patch test of the diffusive fluxes (stencil metric exact on non-smooth meshes, 22 families, the metric-mean
# part, and the conservation of the connection fluxes), and the same mesh as one block and as two blocks joined by a
# connection (one flux per face)
fuss_add_folder_script(PatchTest       numerics-features/patch_test       run_patch.sh 300 "fast;validation;numerics")
fuss_add_folder_script(ConnectionSplit numerics-features/connection_split run_split.sh 300 "fast;numerics")

# Stability limit of the time step (vnn 0.5 for euler, 0.628 for RK3, with the skew factor on rhombohedral cells): a noisy
# constant field on Cartesian and skewed 16^3 boxes stays bounded for 500 iterations just under the limits and stops just over
# them (test/numerics-features/vnn_stability; its decks come from the generator of test/verification/convergence-order)
fuss_add_script(VnnStability numerics-features/vnn_stability/run.sh 300 "fast;numerics")

# ===========================================================================
#  Verification tier: families of grids, run on request (ctest -L verification)
# ===========================================================================
# Observed order of the discretisation on meshes that are not smooth (pois2 and harm2 on rand, jump, shear at N = 16, 32, 64:
# p(L2) >= 1.8, p(Linf) >= 1.6 on the finest pair), and the FAS W and W + IRS solution against the single-grid one
# (<= 1e-6 K) (test/verification/convergence-order/run.sh)
fuss_add_script(ConvergenceOrder  verification/convergence-order/run.sh 1800 "verification;numerics" order)
fuss_add_script(MultigridSolution verification/convergence-order/run.sh  900 "verification;numerics" fas)
set_tests_properties(ConvergenceOrder PROPERTIES PROCESSORS 6)   # its six cases run in parallel

# ===========================================================================
#  Set-up checks, inputs, stops, restart (fast, implementation)
# ===========================================================================
# Checks of the mesh at set-up and the type of the mesh (test/mesh-checks/run.sh): a block one cell thick in k joined to a
# 3-D block, a wedge next to a plane block, a cell folded at a corner (2-D and 3-D), the first ghost layer behind a chimera
# face turned over, and the same mesh with a wall there
fuss_register(MetricChecks 600 "fast;validation;implementation"
              bash ${CMAKE_CURRENT_SOURCE_DIR}/mesh-checks/run.sh $<TARGET_FILE:${solver}> ${_WORK}/MetricChecks)

# A block one cell across stops the set-up (connection to a block of one cell across its face; one cell in a direction whose
# lines are not orthogonal to the others), the exact thin blocks run (test/mesh-checks/thin-blocks.sh)
fuss_add_script(ThinBlocks mesh-checks/thin-blocks.sh 300 "fast;implementation")

# Input checks and stops: the values of input.ini are checked and every error ends the run with exit code 1; input.ini read
# again during the run (ini-error); the errors of bc.txt, time tables, probes, the initial condition and the registry
fuss_add_script(InputChecks  input-io/input-checks.sh  300 "fast;implementation")
fuss_add_script(RuntimeInput input-io/runtime-input.sh 300 "fast;implementation")
fuss_add_script(ExitCodes    input-io/exit-codes.sh    300 "fast;implementation")

# decomposition-map: the rank of every block from a file of INPUT/ (one rank here; MPIDecompositionMap with MPI)
fuss_add_script(DecompositionMap input-io/decomposition-map.sh 300 "fast;implementation")

# Restart: the run split in two writes the files of the run in one go, byte for byte, in every format (the counters stored
# with the solution); the .szplt case needs TecIO
fuss_add_script(Restart input-io/restart.sh 600 "fast;implementation" core)
if (USE_TECIO)
  fuss_add_script(RestartBinary input-io/restart.sh 300 "fast;implementation" szplt)
endif()

# ===========================================================================
#  Unit programs
# ===========================================================================
# Cell volume (CellVolume): a unit check of the volume of Compute_Metric_Tensor, the integral of det J of the trilinear cell,
# against its exact value on cells with twisted faces, a non-convex cell, a collapsed edge, 2-D, wedge and left-handed cells,
# and on a twisted mesh against the volume enclosed by its boundary faces; a cell whose volume has not the sign of det A at
# its centre raises the flag of degenerate or inverted cells
add_fortran_executable(cell_volume_check unit/cell_volume_check.f90)
add_test(NAME CellVolume COMMAND cell_volume_check)
set_tests_properties(CellVolume PROPERTIES LABELS "fast;numerics;unit" TIMEOUT 60)

# Order of rounding (VolatileOrder): a difference stored in a VOLATILE variable and read back is rounded as written, also
# with -ffast-math, which lets the compiler evaluate (q - c1) - (c2 - c1) as q - c2; Compute_Metric_Tensor takes its node
# differences so. The values come from the command line (no folding).
add_fortran_executable(volatile_order_check unit/volatile_order_check.f90)
add_test(NAME VolatileOrder COMMAND volatile_order_check 1 18014398509481984 3)
set_tests_properties(VolatileOrder PROPERTIES LABELS "fast;numerics;unit" TIMEOUT 60)

# ===========================================================================
#  Smoke runs of the cases without a verifier
# ===========================================================================
foreach(_c nominal irs mg1 mg2 mg3 mg_irs fas_v fas_w fas_w_irs)
  fuss_add_test(ResidualConvergence_${_c} numerics-features/residual_convergence/${_c} 300 "smoke" --iters 300)
endforeach()
fuss_add_test(Chimera          numerics-features/chimera                300 "smoke" --iters 300)
fuss_add_test(ChimeraMG        numerics-features/chimera_mg             300 "smoke" --iters 300)
fuss_add_test(AnisotropicPlate numerics-features/anisotropic_Plate      300 "smoke" --iters 300)
fuss_add_test(TimeConvergence  numerics-features/time_convergence/simdir 300 "smoke" --iters 100)
# Material jumps inside a block: the 4-layer liner at vnn 0.5 (series face conductivities and the face-conductivity time step)
fuss_add_test(MultimatLiner    numerics-features/multimat_liner         300 "smoke" --iters 300)

# The FAS keys must be inert outside a cycle: cycle = none with every FAS key set
# must reproduce the plain multigrid run byte for byte
fuss_add_test(MultigridNoneInert numerics-features/residual_convergence/mg3 300 "smoke" --iters 300 --none-inert)

# ===========================================================================
#  MPI tier (needs-mpi): registered only in builds with USE_MPI=ON
#  1 process against 2 ranks, and ranks x threads layouts; the scripts launch the ranks themselves (mpirun)
# ===========================================================================
if (USE_MPI)
  fuss_add_test(ChimeraMPI   numerics-features/chimera                      300 "needs-mpi;implementation" --iters 300 --mpi 2)
  fuss_add_test(ChimeraMGMPI numerics-features/chimera_mg                   300 "needs-mpi;implementation" --iters 300 --mpi 2)
  fuss_add_test(MultigridMPI numerics-features/residual_convergence/mg_irs  300 "needs-mpi;implementation" --iters 300 --mpi 2)
  fuss_add_test(FASMPI       numerics-features/residual_convergence/fas_w_irs 300 "needs-mpi;implementation" --iters 300 --mpi 2)
  fuss_add_test(MultimatLinerMPI numerics-features/multimat_liner          300 "needs-mpi;implementation" --iters 300 --mpi 2)

  # The same output bits on 1x1, 1x4, 2x1, 4x1 and 2x2 (ranks x threads) on a box of 4 blocks and a chimera pair
  fuss_add_mpi_script(MPIBitIdentity       mpi/bit-identity/run.sh        900 "needs-mpi;implementation")
  # MPI calls only from the main thread (MPI_THREAD_FUNNELED): a probe library preloaded into the solver counts the calls made
  # by the other threads (2x2, 2x4 and 2x2 with IRS, ranks x threads); Linux only (LD_PRELOAD)
  if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
    find_package(Threads REQUIRED)
    add_library(mpi_thread_probe MODULE mpi/thread-funnel/mpi_thread_probe.c)
    target_link_libraries(mpi_thread_probe PRIVATE ${CMAKE_DL_LIBS} Threads::Threads)
    set_target_properties(mpi_thread_probe PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
    fuss_add_mpi_script(MPIThreadFunnel mpi/thread-funnel/run.sh 600 "needs-mpi;implementation")
    set_property(TEST MPIThreadFunnel APPEND PROPERTY ENVIRONMENT "FUSS_MPI_THREAD_PROBE=$<TARGET_FILE:mpi_thread_probe>")
  endif()
  fuss_add_mpi_script(MPIRuntimeInput      input-io/mpi-runtime-input.sh  300 "needs-mpi;implementation")
  fuss_add_mpi_script(MPIExitCodes         input-io/mpi-exit-codes.sh     300 "needs-mpi;implementation")
  fuss_add_mpi_script(MPIRestart           input-io/mpi-restart.sh        300 "needs-mpi;implementation")
  fuss_add_mpi_script(MPIDecompositionMap  input-io/mpi-decomposition-map.sh 300 "needs-mpi;implementation")
endif()
