# Changelog

All notable changes to FUSS are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Breaking

- **`vnn` follows MOSE's convention**, without the factor 2 FUSS carried: the stability limits are 0.5 (euler, RK2)
  and 0.628 (RK3). **Halve every `vnn` of an existing deck** to keep its time step (exact: the halved deck reproduces
  its old run bit for bit). The default is 0.15.
- **`res-threshold` is relative** to the residual of the first iteration (of each grid level), as in MOSE; the default
  1e-10 is a fall of ten orders. To keep the stop of an old deck, divide its `res-threshold` by its first residual. A
  time-accurate run no longer ends on `res-threshold`.
- An update that leaves the property table or is not finite stops the run (exit code 1), also with IRS: smoothed
  updates are no longer halved.

### Added

- `diffusive-metric` (`stencil`, default, or `mean`): the stencil face metric makes a linear field exact on any mesh and
  the scheme second order on non-smooth meshes.
- `ini-error` (`warn` or `stop`): what an error in a runtime reading of `input.ini` does.
- `decomposition-map`: the MPI rank of every block from a file of `INPUT/`.
- Exact restart: every solution stores the run counters and the first residual, so a run split in two writes the files
  of the run in one go, in every format.
- Set-up checks with exit code 1 on every rank: unreadable or out-of-rule `input.ini` values, unread keys (warning),
  `bc.txt` indices out of range, invalid property tables, degenerate, folded or turned cells, blocks one cell across.
- The log prints the mesh type and the largest skew factor of the time step.
- CTest tiers `fast`, `validation`, `verification`, `smoke`, `needs-mpi`, 20 new tests (patch test, convergence
  order, multigrid solution, stability limits, restart, exit codes, MPI bit identity), and a GitHub CI workflow.
- `USE_IPO` (default ON, `install.sh --no-ipo`): interprocedural optimization of the RELEASE build, with the archiver
  that can index the objects (`llvm-ar` for ifx, `gcc-ar` for gfortran). CMake refuses ifx 2024.1 at `-O2`/`-O3`
  without IPO (`-DFUSS_ALLOW_IFX_NOIPO=ON` overrides).
- Versioning control: ORION is required as 2.1 or a later 2.x, FiNeR is pinned at v2.3.0; a library off its
  requirement stops the configure (`-DFUSS_PIN_CHECK=WARNING` overrides). See `docs/development/versioning.md`.

### Changed

- Block connections carry one flux per face, given with opposite signs to both cells: a mesh split into blocks gives
  the fluxes of the unsplit mesh.
- A face between two materials inside a block takes the series conductivity of its two half cells, and the time step
  takes the face conductivities.
- The cell volume is the exact volume of the trilinear cell; ghost nodes are extrapolated linearly. Results change in
  the last bits.
- The residual is summed in a fixed order: results do not depend on the number of threads or ranks.
- MPI runs in `FUNNELED` mode with fixed tags; every error stops every rank with exit code 1.
- Text outputs write 17 significant digits; steady solutions carry minus the iteration count as their time.
- ORION is a release checked out by the configuration (in `lib/ORION`, kept in `deps.lock`), not a submodule.
- FiNeR is the submodule `lib/FiNeR`; the configuration clones its five libraries.
- The documentation is built with Zensical; V&V setup and run commands are in the README of each test folder.

### Removed

- The submodules `lib/ORION` and `lib/third_party/*`, the CMake variables `ORION_PATH`, `FINER_PATH`,
  `FINER_DEPS_PATH`, and `install.sh update`.

### Fixed

- A 3-D mesh whose block 1 is one cell thick ran as 2-D.
- The time step of 3-D cells skewed in all three directions exceeded the stability limit by up to 1/8.
- Wall boundary conditions on skewed meshes (tangential term, ghost cells, signs of BC 301 and 302); BC 303 dropped
  `qrad`; the `prim` update used the wrong heat capacity with temperature-dependent `c_p`.
- The wall output now matches the applied flux, so a converged run closes its energy balance.
- Chimera corner ghosts (data race) and the ghosts at the ends of chimera faces.
- Restart from `.szplt` solutions, multigrid restarts, probe rows after a restart, `st.*` read with `ic-format`.
- VTK wall output of multigrid runs (crash or empty file).
- Build of the MPI code with gfortran ≥ 10.
- The gfortran RELEASE build with IPO did not link (plain `ar` cannot index `-flto` objects).
