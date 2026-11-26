# Testing

The cases in `test/` run FUSS on problems with analytical or benchmark solutions, and check the
numerical features, the inputs and the parallel runs. CTest runs most of them; each case can also be
run by hand. The pyrolysis, Darcy, ablating-wall and literature-benchmark cases under `test/pyrolysis/`
are graded by the [gate scripts](#gate-scripts), which `test/CMakeLists.txt` registers as the tests
of label `gate`, outside the pull-request check. The setup of each case is in the `README.md` of its
folder.

---

## Test cases

| Case | Checks | Verifier |
|---|---|---|
| `steady-state/temperature_Plate` | Prescribed temperatures (NAFEMS T9 (i)) | `verify.py` |
| `steady-state/qvol_Plate` | Volumetric heat source (NAFEMS T9 (ii)) | `verify.py` |
| `steady-state/convective_Plate` | Convection (BC 303, NAFEMS T4) | `verify.py` |
| `steady-state/convective_qrad_Rod` | Convection with an imposed flux (BC 303), exact solution | `verify.py` |
| `steady-state/skewed_Plate` | Linear patch test on a 45° sheared mesh, walls 301, 302, 303, 305 | `verify.py` |
| `steady-state/radiative_Rod` | Radiation (BC 304, NAFEMS T9 (iii)) | `verify.py` |
| `steady-state/multilayer_Slab` | Three materials in series inside one block, exact solution | `verify.py` |
| `transient/oscillating_Rod` | Oscillating wall temperature (NAFEMS T5), probes | `verify.py` |
| `transient/multimat_Plate` | Three material layers, Fourier-series solution (LANL) | `verify.py` |
| `transient/variable_cp_Plate` | Temperature-dependent heat capacity, exact solution | `verify.py` |
| `numerics-features/skewed_3D` | Time step of 3-D rhombohedral cells, two materials | `verify.py` |
| `numerics-features/chimera` | Chimera (overset) interface | smoke run |
| `numerics-features/chimera_mg` | Chimera donors on 3 multigrid levels, FAS cycle | smoke run |
| `numerics-features/anisotropic_Plate` | Cells of aspect ratio 100, 3 levels, W + IRS | smoke run |
| `numerics-features/multimat_liner` | Material jumps inside a block (four-layer liner, 2 blocks) | smoke run |
| `numerics-features/residual_convergence/*` | Six convergence strategies and three FAS cycles | smoke runs; `plot.py` |
| `numerics-features/time_convergence/simdir` | Parallel timing study | smoke run; `plot.py` |
| `numerics-features/vnn_stability` | Stability limit of the time step | `run.sh` |
| `numerics-features/patch_test`, `connection_split` | Linear-field exactness; a mesh split into blocks | `run_*.sh` |
| `numerics-features/periodic_sector` | Periodic faces (BC 201): one sector of a helical annular wall, as two blocks and as one block periodic with itself, against the full ring of `101` connections, both face metrics, a wrong pairing stopped at set-up | `run_sector.sh` |
| `verification/convergence-order` | Observed order of accuracy on non-smooth meshes; FAS against single grid | `run.sh` |
| `mesh-checks` | Set-up checks: thin blocks, wedges, folded cells, turned ghost layers | `run.sh`, `thin-blocks.sh` |
| `input-io` | Input errors, runtime input, restart, decomposition maps | scripts |
| `unit/*` | Cell volume, rounding order, table reader | unit programs |
| `pyrolysis/kinetics_isothermal`, `adiabatic0d`, `heated0d`, `conduction` | Pyrolysis of a charring solid on a static mesh: the closed-form kinetics propagator, the 0-D decomposition ODE (adiabatic and volumetrically heated), pure conduction on the graded mesh with the kinetics stripped | `verify.py` |
| `pyrolysis/energy_audit` | Discrete energy and solid-mass balances of the pyrolysis update, closed to round-off; the independent state sum and a vnn refinement | `verify.py` |
| `pyrolysis/ATC1`, `ATC2-chemistryOff` | Ablation Test Cases 1 and 2 (chemistry off) on the TACOT tables, against local PATO runs | `verify.py`, `ref/` |
| `pyrolysis/darcy` | The type-2 (Darcy) gas closure: sealed and open slabs, the analytic steady pressure profile (with and without Klinkenberg), an orthotropic permeability tensor, and the 0-D / ATC2 cases rerun under type 2 | `verify.py` per case |
| `pyrolysis/exact` | Exact-solution benchmarks of the type-2 closure: transient pressure diffusion, steady transpiration, a sharp charring front | `verify_exact.py` |
| `pyrolysis/lit` | Literature benchmarks: Baum & Atreya (2007) charring front with Darcy flow, Zhang (2010) three-code decomposing slab | `verify.py` |
| `pyrolysis/ablation` | Ablating wall (BC `306`/`307`): the surface energy balance, B′ tables, the `307`/`306` identity | see `ablation/README.md` |
| `mpi/bit-identity`, `mpi/thread-funnel` | One MPI rank against several ranks/threads: byte-identical output files; MPI calls only from the main thread | `run.sh` |

A case directory holds `input.ini`, `MESH/mesh.tec`, `INPUT/` and `FUSS.sh`. A `verify.py` reads the
run's `OUTPUT/`, prints `Result: PASS` or `FAIL`, and exits with 1 on failure. The studies keep their
data in `SOLUTION/`, from which `plot.py` draws the figures of the [V&V pages](../vv/index.md).
`numerics-features/convergence_figure/` is not a case: it runs and plots the comparison of
[Multigrid Cycles (FAS)](../vv/fas-cycles.md).

---

## CTest

`test/CMakeLists.txt` registers every case. `test/run_case.sh` runs one in a copy under
`build/test/runs/`, with one thread, so the source tree is never written.

```bash
ctest --test-dir build -L fast                        # 25 tests (26 with TecIO), about a minute
ctest --test-dir build -LE "needs-mpi|verification"   # what the pull-request check runs: 41 tests
ctest --test-dir build -L verification                # ConvergenceOrder and MultigridSolution, and the gates
ctest --test-dir build -LE gate -j 4                  # everything but the gates, a few minutes
ctest --test-dir build -R Qvol -V                     # one test, with its output
```

A build without MPI and TecIO registers the 41 tests of the pull-request check (25 `fast`,
15 `validation`, 15 `smoke`; a test can carry several labels), the 2 full verification cases and
the [gate scripts](#gate-scripts) (`ctest --test-dir build -N -L gate` lists them).

### Labels

The labels are those of MOSE, plus `smoke`:

| Label | Tests |
|---|---|
| `fast` | the 9 verification cases under a minute; `ChimeraEdge`, `Skewed3D`, `MetricChecks`, `ThinBlocks`, `PatchTest`, `ConnectionSplit`, `PeriodicSector`, `VnnStability`; the unit programs; the input tests without MPI |
| `validation` | the 10 verification cases (with `convective_Plate`, about two minutes), `ChimeraEdge`, `Skewed3D`, `MetricChecks`, `PatchTest`, `PeriodicSector` |
| `verification` | `ConvergenceOrder` (6 cores) and `MultigridSolution`, and every `gate` test; on request only |
| `gate` | the [gate scripts](#gate-scripts): one test per physics gate, one per other script; all also labelled `verification`, which keeps them out of the pull-request check |
| `mpi-gate` | with `USE_MPI`: the gate scripts on N ranks through Intel MPI's launcher (`$I_MPI_ROOT/bin/mpirun`); *Not Run* without `I_MPI_ROOT` |
| `smoke` | the cases without a verifier: 300 iterations (100 for the timing case) must end with status 0 and a residual history |
| `needs-mpi` | MPI builds: 2 ranks must write the same files as 1 process, bit for bit; plus the MPI tests below |

Each test also carries a tag: `physics` (a model against its solution), `numerics` (a property of
the discretisation), `implementation` (inputs, stops, restart, parallel runs) or `unit`.

### Pass criteria

- **Verification cases**: the run ends with status 0 and `verify.py` passes.
- **`ChimeraEdge`**: a linear field on a sheared chimera pair is preserved to 10⁻⁹ K, including the
  cells at the ends of the chimera faces.
- **`Skewed3D`**: the time step of every cell matches a replica of the formula and its skew factor.
- **`MetricChecks`, `ThinBlocks`**: the set-up checks stop where they must, with exit code 1 and the
  block named; the valid thin blocks run.
- **`MPIBitIdentity`**: 1×4, 2×1, 4×1 and 2×2 (ranks × threads) write the files of 1×1, on a box of
  four blocks and a chimera pair, steady with IRS and time-accurate.
- **`MPIThreadFunnel`**: a probe library preloaded into the solver counts the MPI calls made by any
  thread other than the primary one. It must find zero on every rank at 2×2, 2×4 and 2×2 with IRS.

### Input, restart and stop tests

The tests of `test/input-io/` run scripts on copies of the cases. Each prints one
`CASE <name> PASS|FAIL` line per check; the test fails if one check does.

| Test | Label | What it checks |
|---|---|---|
| `InputChecks` | `fast` | values out of their rule or unreadable, a probe out of range, a gap in the probe numbers, a missing or empty `input.ini`: exit code 1 with the key and value |
| `RuntimeInput` | `fast` | `input.ini` changed during the run: `ini-error = warn` gives one warning and an unchanged history; `stop` gives exit code 1 |
| `ExitCodes` | `fast` | errors in `bc.txt`, time tables, probes, FAS settings, the initial condition, `properties.dat`, a NaN wall temperature, a restart without a solution: exit code 1, the message, no solution |
| `Restart` | `fast` | a run split in two writes the files of the run in one go (Tecplot ASCII and `vtk raw`, multigrid, FAS, probes, numbered solutions); restarts without counters; a transient run from a steady solution |
| `RestartBinary` | `fast`, needs TecIO | the same split with `.szplt` |
| `DecompositionMap` | `fast` | a map with every block on rank 0 equals the run without a map; every error of the map stops the run |
| `MPIRuntimeInput` | `needs-mpi` | a missing or different `input.ini` in the folder of rank 1 changes nothing |
| `MPIExitCodes` | `needs-mpi` | both ranks print the STOP line and exit with 1, on input, `bc.txt`, initial-condition and runtime errors |
| `MPIRestart` | `needs-mpi` | the split restart on 2 ranks, also with two grid levels |
| `MPIDecompositionMap` | `needs-mpi` | reordered blocks on 2 and 4 ranks write the files of the serial run; a map for 3 ranks stops a 2-rank run |

---

## Gate scripts

`test/CMakeLists.txt` registers every gate script below as CTest tests of label `gate`, also
labelled `verification`, so the pull-request check does not run them. Each test runs in its own
folder `build/test/runs/gates/<test>` (`FUSS_GATE_WORK` names another root), on the tree's
executable (`-b $<TARGET_FILE:FUSS>`). Exit code 77 means that the script skipped a check: a
`verify.py` that printed `SKIP` (a reference file missing), or a requirement
`test/ctest_requires.sh` found missing (Intel MPI, the build's own compiler at its version,
`CMakePresets.json`). CTest reports such a test *Not Run*, never *Passed*.

A build tree on a network file system needs `-DFUSS_GATE_WORK=<local dir>`. The `TIMEOUT` of each
gate is at least 5 times its wall with the work folders on a local disk (at least 300 s), and
several gates write at every iteration: on NFS they run far slower and pass their `TIMEOUT`.
Measured with the same executable: `heated0d_type2` (in `Gate_Darcy0d`) took 35.4 s
of solver time with its work folder on NFS and 0.18 s on a local disk, and
`Gate_DarcyAtc2ChemistryOff` stopped at its 300 s `TIMEOUT` on NFS against 31 s on a local disk.

```bash
ctest --test-dir build -N -L gate                     # the gate tests, by name
ctest --test-dir build -L gate -j 4                   # all of them (the longest about 11 minutes at one thread)
ctest --test-dir build -R '^Gate_Pyrolysis$' -V       # one gate, with its output
cmake -DFUSS_GATE_THREADS=4 build                     # thread count of the two cross-thread legs (default 8, their PROCESSORS)
cmake -DFUSS_GATE_WORK=/var/tmp/$USER-gates build     # work folders on a local disk (default build/test/runs/gates)
ctest --test-dir build-mpi -L mpi-gate                # the MPI forms (mpi preset, I_MPI_ROOT set)
```

By hand, each script runs from any folder, on `bin/FUSS` of the tree unless `-b` (or
`FUSS_GATE_BIN`) names another executable, and several call each other.

| Script | What it proves | Usage by hand | CTest tests |
|---|---|---|---|
| `run_physics_gates.sh` | Runs a fresh copy of each physics-claim case and grades it with that case's own `verify.py`: the pyrolysis, Darcy, ablating-wall and benchmark cases of `test/pyrolysis/`, plus the earlier non-pyrolysis cases. A gate whose `verify.py` skipped a check is printed `SKIP`, and the script then exits 77 | `./test/run_physics_gates.sh [-p NTHREADS] [-w WORKDIR] [-c GATE] [-a] [-b BINARY]` (`-h` lists the gates) | `Gate_<Gate>`, one per gate (`-c <gate>`, the name in CamelCase), and `Gate_RecedingAll` (`-a`) |
| `run_restart_gate.sh` | A run split into two halves with a restart in between must reproduce a continuous run; legs include `pyrolysis/energy_audit` (the reaction state `chi_*` round-trips), `pyrolysis/darcy/mass_balance_1d`, `pyrolysis/darcy/heated0d_type2` (also in `vtk raw`), `pyrolysis/lit/baum_atreya2007/qs_N100`, `pyrolysis/ablation/seb_closure` and `pyrolysis/ablation/coupled307/run307` (the restart state `bc306-state.dat`), among other legs | `./test/run_restart_gate.sh [-w WORKDIR] [-m NRANKS] [-c CASE[:HALF[:TOL[:FORMAT]]]]... [-b BINARY]` | `Gate_Restart`; `Gate_RestartMPI` (`-m 4`, one four-block case) |
| `run_mpi_gate.sh` | A case run on N MPI ranks must write the same output files as the same case on one rank, byte for byte (`wall*.tec` compared with −0.0 equal to +0.0) | `./test/run_mpi_gate.sh [-w WORKDIR] [-r] [-c CASE:NRANKS[:verify]]... [-b BINARY]`; runs `check_build_config.sh --preset mpi` and `check_mpi_tags.py` itself before its cases | `Gate_MPI` |
| `run_determinism_gate.sh` | The same build, run twice with identical settings, must write identical manifests (`run_regression.sh` under the hood); `-x` additionally checks that one case gives the same `OUTPUT/field.tec` at two thread counts | `./test/run_determinism_gate.sh [-p NTHREADS] [-m NRANKS] [-a] [-x] [-X] [-c CASE]... [-w WORKDIR] [-b BINARY]` (`-b` with `-X` only) | `Gate_Determinism` (`-X`: the cross-thread runs alone) |
| `run_guard_demos.sh` | Every setup refusal is run against the configuration it is meant to reject, and both the error text and the exit code (1) are checked; covers the pyrolysis section (IRS, multigrid, `prim`, steady, connection and outflow errors, missing or malformed tables), `[FUSS-Darcy]` (the type-1/type-2 combination, bad keys, no open face, chimera, a bad `p_w`), the ablating wall `306`/`307` (a bad `bc.txt` line, the env table, the B′ table, the restart state `bc306-state.dat`, the `307` exchange contract) and periodic faces (`201`, the type-1 refusal), among other sections; `-m N` repeats every demo on N MPI ranks | `./test/run_guard_demos.sh [-w WORKDIR] [-m NRANKS] [-b BINARY]` | `Gate_GuardDemos`; `Gate_GuardDemosMPI` (`-m 2`) |
| `run_runtime_ini_gate.sh` | `input.ini` re-read mid-run applies only the runtime-mutable keys and ignores (and logs) every other change, including one that would reach a refusal checked only at set-up. Every run is bounded (300 s, 900 s on N > 1 ranks): a run the edits fail to stop is a failure | `./test/run_runtime_ini_gate.sh [-w WORKDIR] [-m NRANKS] [-b BINARY]` | `Gate_RuntimeIni`; `Gate_RuntimeIniMPI` (`-m 4`) |
| `run_init_sweep.sh` | Builds the solver twice with opposite garbage-fill policies for uninitialised locals and runs the suite under each (`run_regression.sh` under the hood); any artefact that differs read an uninitialised variable | `./test/run_init_sweep.sh [-p NTHREADS]` | none (below) |
| `run_hook_interface_check.sh` | The per-stage user hook (`FUSS_Hook_Interface_m`, `external_function_if`) is checked at compile time: a procedure with the wrong signature fails to compile, the right one compiles | `./test/run_hook_interface_check.sh` (`FUSS_MODULE_DIR`: another tree's `modules/`; `FC`: the compiler) | `Gate_HookInterface` |
| `run_morph_guards.sh` | The MORPH component FUSS builds (`MORPH_PATH`) uses no FUSS module, stops nowhere and does no I/O, and exactly one FUSS file uses MORPH's modules; MORPH's library compiles `-Wall` clean and its unit programs give the same output at 1 and 8 threads (MORPH's own `test/run_morph_guards.sh`) | `./test/run_morph_guards.sh` | `Gate_MorphGuards` |
| `check_build_config.sh` | `build/` was configured the way `CMakePresets.json` says (a serial or TecIO-linked build would otherwise run silently) | `./test/check_build_config.sh [--quiet] [--preset NAME]` | `Gate_BuildConfig` (`--preset mpi` in a tree with `USE_MPI`, else `default`); also run by `run_regression.sh`, `run_physics_gates.sh`, `run_mpi_gate.sh --preset mpi` |
| `check_mpi_tags.py` | Every point-to-point MPI call in `src/` carries a named, distinct tag constant in `0..32767` (and below `7001`, which hydra-MF's own tags use) | `python3 -I check_mpi_tags.py --dir <src directory>` or `check_mpi_tags.py <repo> <tree-ish>` | `Gate_MPITags` (`--dir src`); also run by `run_mpi_gate.sh`, before its cases |

Two runs are not registered: `run_determinism_gate.sh` without `-X` (its two full passes) and
`run_init_sweep.sh`. Both run every case in place through `run_regression.sh`, which deletes and
writes again the `OUTPUT/` and `bin/` folders of the cases in the source tree (git ignores them, but
the tree is not left as it was), so they cannot run beside other tests. The sweep also builds FUSS
twice with the toolchain its script names (`/usr/bin/gfortran`, `gcc-ar-12`, `gcc-ranlib-12`); each
build writes `bin/FUSS` of the source tree, where `src/app/CMakeLists.txt` puts the executable, and
its exit trap deletes `build/` and builds it again from the preset. Run both by hand, with nothing
else using the tree.

---

## Tests of the discretisation

### `PatchTest`

`test/numerics-features/patch_test/run_patch.sh`. One forward-Euler step from an exact linear field
$T = A + \mathbf{g}\cdot\mathbf{x}$, whose heat flux is uniform: the flux sum of every cell and the
wall flux error must vanish to round-off.

- **Part 1**, `diffusive-metric = stencil`, 22 mesh families: sheared, random, alternating sizes,
  folded connections, symmetry planes, every wall type in 2-D and on skewed 3-D cells, 3 and 5 blocks
  around a point, connections with a 1 : 4 jump and reversed or exchanged index axes. Every family
  must stay below 10⁻¹¹ of $k|\mathbf{g}|A_\text{cell}$; measured 10⁻¹⁷ to 1.4 × 10⁻¹³.
- **Part 2**, `diffusive-metric = mean`: the sheared mesh stays exact, and one jump of the cell size
  (ratio 1.5) shows the expected bias $(1+r)^2/(4r) - 1 = 6.4201 \times 10^{-3}$.
- **Part 3**, conservation: on a random field with adiabatic walls, the residuals of all cells sum to
  zero within 10⁻¹¹ of their absolute sum, so both sides of a connection see the same flux.

`patch_case.py` generates the meshes (`patch_case.py case <family> <dir> [--metric mean]`).

### `ConnectionSplit`

`run_split.sh`. The random mesh as one block and as two blocks joined by a connection, 50 RK3
iterations: the fields agree to 10⁻¹² relative, and both orders of the records give the same bits.

### `VnnStability`

`vnn_stability/run.sh`. The stability limit of the time step (see
[VNN](../theory/time-integration.md#diffusive-stability-condition-vnn)), on 16³ cells with walls at
1000 K and a random ±50 K start, which holds the highest modes. Serial, no IRS, 500 iterations, on
Cartesian and on rhombohedral cells (edges 115° apart).

| Scheme | Must stay bounded | Must stop |
|---|---|---|
| euler | `vnn = 0.50` | `vnn = 0.51` |
| RK3 | `vnn = 0.60` | `vnn = 0.64` |

The stable runs stay within 0.43 K of 1000 K; the unstable ones stop at iterations 80 to 252. A run
that converges before the budget also fails: just under the limit the highest modes do not decay.

### `ConvergenceOrder`

`verification/convergence-order/run.sh order`, label `verification`. The observed order on
non-smooth meshes, against exact solutions:

- Poisson, $T = 1000 + 100\sin\pi x\sin\pi y$ (source in `INPUT/st.tec`);
- harmonic, $T = 1000 + 100\sin\pi x\,\sinh\pi y/\sinh\pi$;

each on three meshes: `rand` (nodes moved by up to a quarter cell), `jump` (cells 4× wider beyond
x = 0.5, sheared) and `shear` (x += 0.5 y), at N = 16, 32, 64.

Each case passes when the finest pair shows $p(L_2) \ge 1.8$ and $p(L_\infty) \ge 1.6$.

| Solution | rand | jump | shear |
|---|---|---|---|
| Poisson, $p(L_2)/p(L_\infty)$ | 2.00 / 1.93 | 2.00 / 1.98 | 1.99 / 1.96 |
| harmonic, $p(L_2)/p(L_\infty)$ | 2.03 / 1.75 | 1.99 / 1.86 | 1.99 / 1.88 |

The order of the wall flux is reported, not asserted: 0.8 to 1.6. With `FUSS_ORDER_METRIC=mean` the
same run uses the `mean` metric and must fail: the control that the test can fail.

### `MultigridSolution`

`run.sh fas`, label `verification`. Two non-smooth plates at N = 32 (2 × 2 blocks, 3 levels) under
FAS W and W + IRS. They converge in 58 to 169 iterations (single grid: about 5 000) with no halved
update, and differ from the single-grid solution by at most 10⁻⁶ K (measured 1.6 to
2.2 × 10⁻⁷ K).

---

The `needs-mpi` tests are registered only in a tree built with `USE_MPI`: with the `mpi` preset, `ctest --test-dir build-mpi -L needs-mpi`.

## Running a case by hand

```bash
cd test/steady-state/temperature_Plate
./FUSS.sh solve               # serial; -p <n> threads, -m <n> MPI ranks, before 'solve'
python3 verify.py

# MPI (2 ranks): the mpi preset's executable through Intel MPI's own launcher
FUSS_EXE=FUSS-mpi FUSS_MPIRUN=$I_MPI_ROOT/bin/mpirun ./FUSS.sh -m 2 solve
```

`FUSS.sh solve` copies the repository's `bin/FUSS` into the case when it is newer.
`FUSS_EXE` selects the executable a case copies from `bin/` (default `FUSS`; `FUSS-mpi` for the `mpi` preset, which
`./FUSS.sh compile` then builds with `install.sh compile --preset=mpi`), `FUSS_MPIRUN` the launcher (default the first
`mpirun` on `PATH`). The regression runner does both itself: `test/run_regression.sh -m N` runs every case with
`bin/FUSS-mpi` on N ranks (N = 1: one process, no launcher; N > 1 under `timeout 900`) and writes its manifest to
`test/.regression/omp<threads>-mpi<N>/`.
`run_case.sh <case> <work dir> <binary> [--iters N] [--verify] [--mpi N]` caps a run at N iterations
by rewriting `iter-threshold` and every `level<m>-iter` of the copy.

---

## Adding a case

1. Create the directory with `input.ini` and the mesh. The `[GRIB-*]`, `[GPB-*]`, `[ICB-*]` and
   `[BCB-*]` sections describe how ATLAS generates `MESH/` and `INPUT/`; the first line of
   `phase.txt` must be `solid-bulk`. A case whose data ATLAS cannot express carries a script that
   writes them (`skewed_Plate/make_inputs.py`).
2. Write `verify.py`: read `OUTPUT/`, compare with the reference, print `Result: PASS` or `FAIL`,
   exit 1 on failure.
3. Write `README.md`: setup, check and run command.
4. Register the case in `test/CMakeLists.txt` with
   `fuss_add_test(<Name> <dir> <timeout> "<labels>" --verify)`; a pyrolysis, Darcy or ablating-wall
   case instead gets a `gate_<name>` function and a case name in `run_physics_gates.sh`, and the
   entry `<Name>:<name>:<timeout>` in the list `_PHYSICS_GATES` of `test/CMakeLists.txt`
   (the test `Gate_<Name>`).
5. Add a page to `docs/vv/` (problem, result, assessment) and to the `nav` of `mkdocs.yml`.

### The physics-gate convention (`test/pyrolysis/`)

A case graded by `run_physics_gates.sh` (a `Gate_` test under CTest) follows the same `verify.py` contract (read
`OUTPUT/`, compare with the reference, print `Result: PASS` or `FAIL`, exit 1 on failure), plus:

- no tolerance is a FUSS result of an earlier run plus a margin. A result checked against an exact
  solution passes by the convergence rule (`test/pyrolysis/tools/conv_rule.py`): on a ladder of runs,
  coarsest first, the error falls on every rung and its observed order on the finest pair lies in
  $[0.8p, 1.2p]$ for the scheme's expected order $p$. A result checked against other codes passes within
  their own disagreement, computed by `verify.py` from the committed reference files each time it runs:
  the PATO-Amaryllis spread per quantity (ATC 2.1, 2.2, 2.3; borrowed from ATC 2.1 where only PATO
  solved the variant), the PATO-FIAT spread (ATC1), the three codes' spread (Zhang 2010). A
  reference-derived tolerance is proven red through FUSS's output, never through the reference: the
  tolerance moves with the reference files it is computed from;
- an `expected.json` beside a verify holds only what the rule needs from outside the run (the expected
  orders, a round-off floor), never a measured FUSS value; there is no `--record`;
- an identity between two FUSS runs (restart, threads, blocks, a mirror) takes a bound from an argument
  where one exists (the run's own round-off, $N_{\rm steps}\times N_{\rm cells}\times\varepsilon$ for a
  budget sum; the Newton stop tolerance per step for the type-2 restart; the surface-balance tolerance
  for the mirrored 307 faces). Where none exists the bound is a floor set from measured draws, and the
  code says so: the type-2 0-D non-uniformity ($10^{-5}$ K), the restart $B'_{pg}$ under type 2
  ($10^{-9}$), the two-block smoothing seam ($10^{-10}$ m, $10^{-7}$ of $T$) and the periodic Darcy
  sector ($10^{-10}$);
- the gate's own `FUSS_GATE_SELFTEST=1` is a proven-red hook, not a correctness check: it runs the
  case with its tolerance deliberately broken (a sedded copy of `verify.py`, or an input key edited
  in a work copy) and requires the run to FAIL, so that a gate nobody has ever seen fail is not
  mistaken for one that works. `verify.py` scripts that take the variable honour
  `FUSS_GATE_SELFTEST=1` by zeroing their own tolerances (and emptying any order window) rather than
  reading it from `run_physics_gates.sh`; `energy_audit`, `ATC1` and `ATC2-chemistryOff` have such a
  hook too.
