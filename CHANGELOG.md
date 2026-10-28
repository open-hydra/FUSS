# Changelog

All notable changes to FUSS are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

This entry covers the sync of October 2026. The 78 commits between `main` (830f228) and the first commit of the sync
(FAS cycles and their keys `cycle`, `pre-iter`, `post-iter`, `coarse-iter`, `prolongated-diter`; the series face
conductivity and the face-conductivity time step; wall boundary conditions 301-305; probes; the MPI exchange) are not
listed here: the pages of `docs/` describe them.

### Added

- The key `diffusive-metric` of `[FUSS-Numerics]` (read at set-up only; an edit during the run is ignored) chooses the face
  metric of the conduction fluxes: `stencil`, the default, takes it from the same differences of the cell centres that the
  face gradient takes of the temperature, so that a linear field has the exact flux on any mesh; `mean` takes the arithmetic
  mean of the two cell metrics at the interior faces, the metric of FUSS before this release, which overstates the flux by
  (1+r)²/(4r) − 1 at a cell-size ratio r (4.2 % at r = 1.5, 8.9 % at 1.8, at any resolution). Under `mean` the block
  connections keep the two-sided balance with the cell metrics and the walls the cell metric. A deck that names
  `inverse-mean` stops at set-up (see Removed).
- With MPI, rank 0 reads `input.ini`, at set-up and at every re-read, and passes its text to every rank. A file that cannot be
  opened, is empty or cannot be read stops every rank with exit code 1 (every rank read its own copy, and a missing one ran on
  defaults until the first required key).
- The key `ini-error` of `[FUSS-IO]` (`warn`, the default, or `stop`) decides what an error of `input.ini` read again during
  the run does: with `warn` the reading is not taken and the run goes on with every value it had, on every rank; with `stop`
  every rank stops with the message and exit code 1. At set-up every error stops the run.
- The key `decomposition-map` of `[FUSS-IO]` names a file of `INPUT/` that gives the MPI rank of every block; rank 0 reads it
  and every rank checks it; ranks without blocks idle with a warning. `none`, the default, assigns the blocks as before. The map
  changes the rank of each block, not the solution.
- At set-up FUSS reports, with the line number, every `key = value` of `input.ini` that the INI reader does not give in its
  section (a line after a malformed header, for instance; such a key is not read and its parameter keeps its default). A
  misspelled name in a well-formed section is not reported. A value that does not read as the type of its key, a number that is not finite, a
  list of the wrong length or a logical other than true/false is an error of the reading, and the rules of the registry
  (bounds, intervals, lists, logicals) are checked in every form. Inputs whose keys are all read do not change.
- Every field solution stores the counters of the run (the global iteration count `CYCLE`; `FUSS_ITERATION`,
  `FUSS_OUTPUT_NUMBER`, `FUSS_GRID_LEVEL`, `FUSS_LAST_OUTPUT_TIME`) as vtk field data or Tecplot data set auxiliary data, and a
  restart reads them back: a run split in two writes the same files as the run in one go, byte for byte, in every format
  (`.szplt` included), within the limits the guide lists (a FAS run cut inside a cycle is not exactly restartable). A Tecplot
  ASCII field ends with five `DATASETAUXDATA` records; a script that takes every number after the last zone must stop at them.
  Results of a run do not change; the field files gain those records.
- The solution due at the last iteration (by `sol-diter` or `sol-dtime`) is written with its number, wall and diagnostic files
  before the final `field.*`. Every row of `residual-history.dat` and of the probe files is flushed when written, and the files
  are closed at the end.
- The set-up summary prints the type of the mesh (1-D, 2-D plane, 2-D axisymmetric wedge of an angle, 3-D) and the largest
  skew factor of the explicit step (`skew_factor_max`; exactly 1.000000 on every mesh of the test suite).
- CTest: the tiers are MOSE's labels `fast`, `validation`, `verification` (full forms of verification cases; none registered
  yet) and `needs-mpi`, plus `smoke` (`ctest -L <tier>`), with the tags `physics`, `numerics`, `implementation` and `unit`; the
  registrations go through `fuss_register`, `fuss_add_test`, `fuss_add_script`, `fuss_add_mpi_script` and
  `fuss_add_folder_script` of `test/CMakeLists.txt`. New tests: `InputChecks`, `RuntimeInput`, `ExitCodes`, `DecompositionMap`, `Restart`, `RestartBinary` (needs TecIO) and, with
  `USE_MPI=ON`, `MPIRuntimeInput`, `MPIExitCodes`, `MPIRestart`, `MPIDecompositionMap` (input, stops, restart, decomposition
  map); `MPIBitIdentity` (1x4, 2x1, 4x1 and 2x2 ranks x threads write the same files as 1x1, on a sheared box of four blocks
  plus a chimera pair, steady with IRS and time-accurate) and `ChimeraEdge` (parallel); `Skewed3D`, `MetricChecks`,
  `CellVolume` and `VolatileOrder` (metrics); `PatchTest` (linear-field patch test of the conduction fluxes, 8 mesh families
  and the `metric-mean` part) and `ConnectionSplit` (the same mesh as one block and as two blocks joined by a connection)
  (face metric). 18 tests; the suite had 30 (25 without MPI).
- A GitHub workflow `.github/workflows/ctest.yml`, modelled on MOSE's: the pull-request tier under gfortran (`ctest -LE
  "needs-mpi|verification"`, a build without MPI and TecIO) and the MPI tier under gfortran + Open MPI (`ctest -L needs-mpi`);
  the labels carry MOSE's second-level tags `physics`, `numerics`, `implementation`, `unit`; the MPI script tests reserve 4
  processors under `ctest -j`. This file.
- Set-up stop for a block one cell across, on every grid level, before the ghost metrics are built (`Check_Thin_Blocks`,
  exit code 1 on every rank): a block connection (101) whose block across has fewer than `gc` = 2 cells across its face (as MOSE
  2e31d8a; a build with bounds checking stopped in `BC_Connect_Metrics` with an index out of bounds, a release build read past
  the node array and ran on with a wrong metric, or stopped on a cell with zero corner nodes), and a block with one cell in
  an index direction (not k of a 1-D or 2-D mesh) whose lines are not orthogonal to the other directions (|cos| > 1e-3: the
  face gradient takes no difference along a direction with one cell, so a sheared block one cell across gave a wrong wall heat
  flux, 0.16 to 0.30 of k|g|, with no message). Orthogonal thin blocks run as before; no deck of the test suite changes.

### Changed

- **The diffusion number `vnn` follows the convention of MOSE**: the local step of a cell is
  Δt = vnn ρc_p / Σ_d k_d |∇ξ_d|² (cell pass: k_d = k; face pass: k_d the larger of the two face conductivities of direction
  d), without the factor 2 that FUSS carried. The stability limits are vnn ≤ 0.5 for forward Euler and RK2 and 0.628 for SSP
  RK3 (they were 1 and about 1.25: the same steps). **A deck written before this change runs at twice its step unless its
  `vnn` is halved** (vnn_new = vnn_old / 2; halving is exact in binary, so a halved deck reproduces its old run bit for bit):
  0.5 → 0.25, 0.8 → 0.4, 1.0 → 0.5, 1.5 → 0.75, 2.0 → 1.0. All 24 decks of the test suite and the test scripts were halved
  with this change (the stability checks of the stress battery are 0.6 and 0.65, were 1.2 and 1.3). The default of `vnn` is
  0.15 (it was 0.3: a deck without the key keeps its step). A runtime re-read of `vnn` is in the new convention. `vnn` appears
  nowhere else in the code (`vnn-rise-threshold` is an iteration count, `irs-beta` is unrelated). **Decks outside this
  repository** follow the same rule: halve every `vnn` before running them.
- Interior faces and block connections (101) take the stencil metric J⁻¹, J = [D_a x, D_b x, D_c x] of the cell centres (see
  Added): exact for a linear field on any mesh, skewed, stretched, with jumps of size, folded at a connection. The linear-field
  patch test fell from 1.2e-2 … 1.2e-1 of k|g|A to ≤ 7e-14 on every family. The walls 301-305 take the metric of the
  differences their gradient takes: wall flux and `qw` from up to 1.2e-1 to ≤ 1e-13 of k|g|. Chimera faces (102) are unchanged.
  Results change at round-off on Cartesian meshes: measured on the integrated tree against the 17-digit reference of the
  input/io theme, max |ΔT| ≤ 1.2e-11 K on the 15 short decks of the test suite (one unit of the last printed digit at 2000 K),
  3.0e-11 … 4.3e-11 K on the eight runs of 27 000-295 000 iterations (`nominal` is capped at 100 000 by the gate: `nominal`, `irs`, `mg1`, `mg2`, `mg_irs`, `fas_v`,
  `fas_w`, `fas_w_irs`: round-off accumulated over the iterations) and 3.4e-10 K on `mg3`, all within the declared 1e-7 K; and
  wherever cells are skewed, stretched or non-smooth.
- A block connection carries one flux, computed by the owner record (the cell on the low-index side of the face) and given with
  opposite signs to both cells: the two shares cancel exactly, also between MPI ranks (one double per face per stage,
  `TAG_FC`), and a mesh split into blocks gives the fluxes of the unsplit mesh when the index directions agree and the face
  nodes have the same bits (`ConnectionSplit`: 0.0 relative, it was 4.6e-5). The conductivity at a connection follows the
  interior rule (k at the mean temperature, series k between two materials) instead of the two-sided balance of the temperature
  at the face: the same flux for constant k, O(dT² k″) different with k(T) and one material.
- The volume of a cell is the integral of det J of its trilinear map, V = 1/6 Σ_d e_d · (S_d⁻ + S_d⁺), instead of five
  tetrahedra in absolute value: exact on every cell (twisted faces included), and the volumes of a block sum to the volume
  enclosed by its boundary faces. The node differences go through a VOLATILE array, so that the optimiser
  (`ifx -fp-model fast`, `gfortran -ffast-math`) cannot rewrite a difference of differences at the size of the coordinates. On
  cells with plane faces (all the cells of the suite) the two volumes are equal in exact arithmetic and differ by at most 4e-16
  relative in rounding; the volume enters every update, so **every case changes in its last bits** (measured with the ghost
  change below on the integrated tree against the 17-digit reference: the temperature field of 12 of the 24 decks moves by at
  most 2.5e-12 K and `chimera` by 1.1e-11 K, 9 decks keep their field and move in `diagnostic.tec` and the residual history
  only, 2 transient decks are byte-identical). On cells with twisted faces the exact volume differs from the tetrahedra by up to 10.8 % per cell (3-D patch mesh
  `box3d`).
- The ghost nodes behind the faces that extrapolate (everything but block connections and symmetry faces; chimera included) are
  extrapolated linearly along the index lines, N(−g) = 2N(−g+1) − N(−g+2), on every face. The cubic rule used until now read
  out of bounds on faces 5-6 of a 3-D mesh and with one cell in i. The only reader of an extrapolated ghost metric is the flux
  of chimera faces (102): where the spacing across the face is uniform, as at both chimera faces of the test cases, the two
  rules agree in exact arithmetic and the first ghost metric moves at round-off (|ΔM|/|M| ≤ 2.1e-13 on `chimera`, 2.2e-14 on
  `chimera_mg`: those two cases change in the last bits of the temperature, 1.1e-11 K and 1.8e-12 K on the integrated tree;
  every other case is unchanged by this change). Ghost metrics are computed on the rank that owns the block only.
- The size of every cell, interior or ghost, is checked at set-up (|det A| > 1e-12 |A1||A2||A3|), and the first ghost layer
  behind a chimera face must keep the sign and the size of the first cell. A failing cell writes block, cell, face, the ratio
  dl₂/dl₁ and the corner nodes (once per rank), takes the identity as metric, and stops every rank with exit code 1. A chimera
  face whose edges along it grow 3 times or more across the first cell now stops at set-up and says so; a wall there runs. A
  cell folded at a corner or turned against its block, or whose volume has not the sign of det A at its centre, stops the run
  at set-up on every grid level, naming block, cell and corner node (the former `stop` of a volume ≤ 0 had exit code 0). The 24
  cases of the suite have none; non-convex cells of a generated mesh (a node moved twice by a generator) do, and did not stop
  before. A degenerate stencil (the three differences of a face do not span the space) stops at set-up naming block, cell, face
  and the differences. A block one cell wide across a block connection (101) now stops at set-up (its second ghost layer read
  nodes out of bounds, silently before).
- The residual written to `residual-history.dat` (and tested against `res-threshold`) is summed in a fixed order, the same for
  any number of OpenMP threads, MPI ranks and any split of the blocks over the ranks: each thread sums whole i lines of a block,
  the lines are added in turn, and the root adds the blocks in order. The last printed digit (E20.10) of a few rows of a run may
  differ from the previous version, also on one thread and one rank (2 rows of 100001 in `residual_convergence/nominal`, 1 of
  79818 in `mg_irs`; relative difference below 1e-9). A run with OpenMP reduction-order effects (`nominal`, `mg_irs` on 4
  threads) now writes the same `residual-history.dat` as on one thread. Fields are unchanged.
- The MPI exchange tags are the constants `TAG_T`, `TAG_TG`, `TAG_CHIMERA`, `TAG_R`, `TAG_FC`, `TAG_IO`, `TAG_MATID` (0-6), all
  below 32767, the least `MPI_TAG_UB` the standard guarantees; the sender rank, rank + 10000…40000 and the block number are no
  longer used as tags. No result changes. The boundary fluxes of a rank run over the records of the blocks it owns only (lists
  by block and face built at set-up); they scanned the whole file on every rank and stage. Bit-identical for a `bc.txt` sorted
  by block and face.
- Every stop on an error ends with exit code 1 on every rank: the checks of rank 0 broadcast their verdict, decisions taken on
  shared data stop together after a barrier and `MPI_FINALIZE`, errors found by one rank are reduced, and the output units are
  flushed before `MPI_ABORT` (two stops ended with exit code 0, `error stop 'text'` with 128 under ifx, and the other ranks went
  on into the set-up after a root-only check). bc.txt errors stop every rank with the file, the line and the record (a bad line
  was printed and ignored, the property lines were never checked); time tables with a non-increasing time or a value that is
  not finite stop together. `Admissible_T` tests the bits of T (`Not_Finite`): the same answer on every finite T.
- A steady solution is stamped with minus the global iteration count in every format (Tecplot `SOLUTIONTIME`, vtk `TIME`, which
  was −1), and a restart takes a solution as steady when its time is negative or its Tecplot zones say `STRANDID = 0` (files
  written before). A time-accurate restart from a steady solution starts at t = 0 with new counters (it started at t = the
  iteration count from Tecplot, at t = −1 from vtk); a steady restart from a time-accurate solution counts its iterations from
  the integer part of the time, as MOSE. A restart from a solution without a time or whose time or counters are not finite
  stops at set-up with the file, the value and the remedy. A multigrid restart from `field.*` goes on at level 1 (it started
  again from the coarsest level); the numbering of `field<n>.*` is restored from the solution (the scan of `OUTPUT/` stays as
  the legacy path for a solution without counters). Time-based probe rows fall on the exact multiples of `dtime` (an integer
  counter). The rows of `residual-history.dat` and of the probe files written after the solution a restart reads are dropped,
  a row cut in the middle too, and the files end with a newline (rows were doubled after a restart from an earlier `field<n>`).
  On a restart `INPUT/st.*` is read with `ic-format` (it was read with `sol-format`: a restart in `vtk raw` of a case with
  `st.tec` ran with qvol = 0 in silence). The `.vtm` of a field lists its block files relative to its own folder
  (`vtk/fieldBlock<n>.vts`; it listed `../vtk//fieldBlock<n>.vts`); a restart reads them next to the `.vtm`, or under
  `OUTPUT/vtk/` for files written before.
- ORION main a978930 (PR #7 merged; was 5a1234c), the pin of MOSE: the multi-block VTK reader and writer handle field data
  (`TIME`, `CYCLE`, named scalars), any number of blocks and surface blocks; the Tecplot reader and writer handle the solution
  time with its sign, strands and data set auxiliary data, and the `.szplt` reader returns the solution time. The text formats
  (`tecplot ascii`, `vtk ascii`) write every 64-bit value with 17 significant digits, the fewest that give it back bit for bit:
  every text output changes, each real two characters longer, and a restart from `tecplot ascii` reproduces the field exactly
  (it reproduced it to 1e-11 K). The doubles written are the same as before: the old text is their 15-digit rounding (proved on
  the 662 `.tec` files of the gate runs: 11 585 956 values, 0 mismatches; `vtk raw` and `vtk binary` are unchanged). The Tecplot
  zones of a steady solution carry minus the global iteration count as `SOLUTIONTIME` and are static (no `STRANDID`), those of
  a time-accurate solution the strand of their block (`STRANDID = b`).
- PENF, FACE, FLAP, BeFoR64 and StringiFor, the libraries of FiNeR v2.3.0, are submodules of FUSS under `lib/third_party/` at
  the tags of MOSE (PENF 4cb9d32, FACE da7b740, FLAP a3741e6, BeFoR64 c90ce95, StringiFor b12f62b), added by CMake before FiNeR
  from `FINER_DEPS_PATH` (`lib/third_party/` by default). After pulling this change into an existing clone run
  `git submodule update --init`. CMake builds ORION, FiNeR and the libraries in `${CMAKE_BINARY_DIR}/lib/` (they were built in
  `build/lib/` of the source tree whatever `-B` said) and stops with a remedy when FiNeR is older than v2.1.0. No value of any
  output changes by this (the Fortran sources of the two sets are identical).
- Integration: `origin/main` 830f228 (PR #2, Giacomi: bc 104) is merged. Its one hunk, `Mod_Allocate_Data.f90` (the remote-donor
  mask no longer keeps the donors of 104, which belong to the other solver), is already in this line in another form
  (`case (102)` with `is_local_block`), so the merge keeps this line's text and changes no code (tree identical to the line
  before the merge).
- Integration: after the stops were moved to the routines of every rank, `Check_Mesh_Type` calls `mpi_stop_all`,
  `Stop_On_Bad_Metric` `mpi_stop_if_any`, and `Check_Folded_Cells` returns a flag consumed by `mpi_stop_if_any` in
  `Setup_Metrics`; the `vnn` of the test scripts of the parallel and input/io themes was halved (bit-identity and chimera-edge
  1.0 → 0.5, divergence probe 1.4 → 0.7, transient-from-steady restart 0.5 → 0.25, runtime bad-value 0.4 → 0.2). `TAG_FC` is
  declared once, in the tag table of `Mod_GhostExchange`.

- An update of a cell temperature that is not finite or leaves the properties table, also with IRS, stops the run with exit
  code 1, as an explicit update does and as MOSE (decision D2, 354bf28). The rank that finds the cell writes
  `[RANK r] STOP:` and ends with exit code 1 (the check runs inside an OpenMP loop, where no MPI call is allowed) and the
  launcher ends the other ranks. Until now an IRS-smoothed update that failed was halved up to 10 times and else dropped, with an
  `IRS:` report every `shell-diter` iterations; no case of the test suite or of the FAS acceptance set used it. The two guards of
  the FAS cycles stay (a prolongated correction that leaves the table, a coarse smoothing update under a forcing: halved or
  dropped, `FAS:` reports); the counter of the halvings is `FAS_HALVINGS` (it was `IRS_HALVINGS`).

### Removed

- The inverse-mean face metric (the `inverse-mean` value of `diffusive-metric`, a24b17c, and the metric that replaced the key,
  88fb0e1; neither released) and the stored inverse cell metric `Minv` (9 doubles per ghost-inclusive cell; the ghost layers
  were never touched, so resident memory does not fall: +116 B/cell in 2-D, +114 in 3-D for the boundary metrics per block
  face). The widths of the series conductivity and of the face time step are computed where they are needed (same bits).
- `install.sh` no longer clones PENF, FACE, FLAP, BeFoR64 and StringiFor (they are submodules, see Changed).

### Fixed

- A 3-D mesh whose block 1 is one cell thick in k ran as a 2-D mesh (`ndir` 2): the k-normal gradient of every connection on a
  k face was zeroed, so block 1 was decoupled from its neighbours across it, and the time step ignored the k direction in every
  block. The type of the mesh comes from the nodes of every block: `ndir` is 3 when any block has more than one cell in k, 2
  when any has more than one cell in j, otherwise 1. The 2-D plane or wedge decision looks at every node of every block (the
  k = 1 node plane is the k = 0 plane translated, or turned about the x axis by the same angle in every block, to 1 % of the
  thickness of the block): blocks nearer a turn than a translation (or the converse), wedges of different angles, and wedge and
  plane blocks together stop the run at set-up and say so. No case of the suite changes (all 2-D plane).
- The explicit step of 3-D cells skewed in all three pairs of directions (c12² + c13² + c23² + 2 c12 c13 c23 > 1, e.g. a
  rhombohedron with edges 115° apart) takes the spectral radius of the operator, which exceeds 4α Σ|∇ξ_d|² there by up to 1/8:
  the sum of the step is multiplied by (2 + m + 1/m)/4 (exact on rhombohedra), in the cell pass and in the face-conductivity
  pass. At the stability limit such a mesh ran at an effective vnn 3.7 % above the limit and diverged. The step of every other
  cell, and of every 2-D mesh, is unchanged to the last bit (the factor is exactly 1).
- The metrics of the ghost cells of the symmetry faces (300) are filled one layer at a time: on a direction with one cell (the k
  faces of a 2-D mesh) the second layer of a face copied the first-layer ghost cell of the opposite face, which another thread
  wrote in the same parallel loop (a data race at set-up). Nothing reads those layers: no result changes.
- MPI is called from the primary thread only, as `MPI_THREAD_FUNNELED` requires: the exchange was inside `!$omp single`, which
  any thread of the team may run. `Fill_Ghost_Cell` and the residual extrapolation of IRS use `!$omp master` plus a barrier, and
  on a single rank skip the exchange and its barriers. No result changes.
- Chimera faces that meet at an edge of a block no longer write the same corner ghost cell from different records (a data
  race), and the ghost row beyond the end of a chimera face is extrapolated along the row (2a − b) instead of copied, which
  halved the tangential difference there: a linear field on a sheared chimera pair was disturbed by 0.39 K at the cells at the
  ends of the chimera faces after 5 iterations, it is now preserved (`ChimeraEdge`). Cartesian chimera meshes (the suite's) are
  unchanged or change in the last bits at those cells.
- gfortran ≥ 10 builds of the MPI code: the reductions of `Mod_MPI` take a copy as send buffer (no `MPI_IN_PLACE`, which
  gfortran does not recognise with `-fno-underscoring`), and `Mod_MPI.f90` is compiled with `-fallow-argument-mismatch`.
- A restart from a `tecplot binary` (`.szplt`) solution continues from its iteration count or time; it counted from 0 (ORION's
  binary reader dropped the solution time); now exact (see Added).
- Two stops ended with exit code 0 and, under OpenMP threads, `error stop` inside a parallel region gave exit codes 0 and 134
  and lost messages: all go through the routines of every rank (see Changed).
- `test/transient/variable_cp_Plate/verify.py` ignores `DATASETAUXDATA` records after the last zone.
- The vtk writer of the wall files read `Ni`, `Nj`, `Nk` of the wall blocks, which
  `Initialize_Wall_File` never set: with ORION a978930 a vtk wall output of a multigrid run (and, nondeterministically, of a
  single-grid run: `wallB2F2.vts` differed between two identical runs) crashed with SIGSEGV in `reshape` under gfortran-12, or
  wrote an empty file: under ifx the unset sizes read as zero, so every wall `.vts` of a `vtk` run had the extent 0 0 0 0 0 0,
  one point and no data (838 bytes), while the Tecplot wall zone of the same face has its full extent (nominal: 1 × 11 × 2
  nodes). The three sizes are now set from the mesh array, as MOSE does: the wall `.vts` files carry the face (1 × 11 × 2,
  1490 bytes on that case); field files, `.vtm`, Tecplot and every other output are byte-identical to before (proof: 3 decks,
  `nominal`, `chimera_mg`, `anisotropic_Plate`, in `vtk raw` with both binaries, 77 files, only the 34 wall `.vts` differ).
