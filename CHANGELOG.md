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
- Periodic faces, BC `201` (the record of a `101`, as ATLAS/BCB's `type = periodic` and AFUS write it): a block
  connection whose cells of the other face are moved by the rigid motion of the pair of faces, fitted at set-up to
  the records (no axis or angle in input), for both `diffusive-metric` values. Faces that are not images by a rigid
  motion, or periodic one way only, stop the run. Tests `PeriodicSector` and `PeriodicSectorMPI`.

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

### Pyrolysis, ablating walls and mesh motion (branch `worktree-porous-pyrolysis`)

In this branch only; upstream `main` up to `4f5673c` is merged in (2026-11-17 to 2026-11-21).

- Pyrolysing porous materials, `[FUSS-Pyrolysis]`: `gas-closure = type1` (quasi-steady one-dimensional outflow of the
  pyrolysis gas) or `type2` (Darcy pressure, `[FUSS-Darcy]`).
- Ablating wall, BC `306`: surface energy balance with B′ tables, `[FUSS-Ablation]`; the wall recedes with
  `[FUSS-MeshMotion] law = nodeshift`. Coupled ablating wall, BC `307`: the same balance, with the exchange values given
  by hydra.
- Mesh motion through MORPH (the submodule `lib/MORPH`): `[FUSS-MeshMotion] law` = `static`, `prescribed`,
  `translation` or `nodeshift`; the solution is remapped onto the moved mesh (conservative, first order).
- MPI for all of the above: pyrolysis and Darcy, the ablating walls, the receding mesh.
- The combinations refused at set-up, each with a message (pyrolysis with IRS or with more than one grid level,
  type-1 pyrolysis with periodic faces, a receding wall with multigrid, among others), are listed in
  `docs/user/input.md`.
- Periodic faces (`201`) with these features: a moving mesh whose motion takes a node of a periodic face and its image
  to images (translation along the axis of the pair, prescribed motions of that kind, the receding wall of
  `law = nodeshift`, whose nodes on a periodic seam are joined with their images), checked at set-up and every step;
  type-2 pyrolysis across periodic faces when the permeability tensors are invariant under the pair's rotation.
- The hook that a coupling driver (hydra) passes to the solve has an explicit interface: one argument,
  `type(FUSS_domain_type), intent(inout) :: domain`.
- A coupled face (`103`) is accepted under pyrolysis: type 1 may end its outflow lines there, type 2 opens it by default
  (`open-bc-types` default `303 305 306 103`, pressure `p-wall`); the gas leaves it at the cell temperature. The gas is
  not passed to the gas solver: its mass and enthalpy leave FUSS's energy budget only, and a notice at set-up says so.
  Chimera faces (`102`, `104`) stay refused.

Where this branch behaves differently from upstream:

- An unknown key in a FUSS section of `input.ini` stops the run:
  `[ERROR] unknown option in input.ini: [<section>] <key>`.
- A block whose cells are all left-handed stops the run, with a message that names the block (upstream accepts it).
- Under pyrolysis, every face inside a block takes the mean of the two cells' conductivities, also between two
  materials; a block connection takes the mean for one material and the two half cells in series between two
  materials (no test case runs this last branch yet).
- A steady run with a 301/302 time table warns at set-up that the table is held at its value at time −1. A 304/305
  wall whose temperature does not converge in 50 Newton steps warns (the first time on each rank; the count is printed
  at the end of the run); as upstream, the last iterate is used and the run goes on.
- Build: `-fno-underscoring` is set only without MPI (the MPI Fortran library exports the underscored names);
  `install.sh` takes `--include-morph=<path>` with `build` and `--preset=NAME` with `compile`.
