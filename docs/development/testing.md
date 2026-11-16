# Testing

The cases in `test/` run FUSS on problems with analytical or benchmark solutions, and check the
numerical features, the inputs and the parallel runs. CTest runs them all; each case can also be run
by hand. The setup of each verification case is in the `README.md` of its folder.

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
| `verification/convergence-order` | Observed order of accuracy on non-smooth meshes; FAS against single grid | `run.sh` |
| `mesh-checks` | Set-up checks: thin blocks, wedges, folded cells, turned ghost layers | `run.sh`, `thin-blocks.sh` |
| `input-io` | Input errors, runtime input, restart, decomposition maps | scripts |
| `unit/*` | Cell volume, rounding order, table reader | unit programs |

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
ctest --test-dir build -L fast                        # 24 tests (25 with TecIO), about a minute
ctest --test-dir build -LE "needs-mpi|verification"   # what the pull-request check runs
ctest --test-dir build -L verification                # ConvergenceOrder and MultigridSolution, about a minute on 6 cores
ctest --test-dir build -j 4                           # everything, a few minutes
ctest --test-dir build -R Qvol -V                     # one test, with its output
```

### Labels

The labels are those of MOSE, plus `smoke`:

| Label | Tests |
|---|---|
| `fast` | the 9 verification cases under a minute; `ChimeraEdge`, `Skewed3D`, `MetricChecks`, `ThinBlocks`, `PatchTest`, `ConnectionSplit`, `VnnStability`; the unit programs; the input tests without MPI |
| `validation` | the 10 verification cases (with `convective_Plate`, about two minutes), `ChimeraEdge`, `Skewed3D`, `MetricChecks`, `PatchTest` |
| `verification` | `ConvergenceOrder` (6 cores) and `MultigridSolution`; on request only |
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

## Running a case by hand

```bash
cd test/steady-state/temperature_Plate
./FUSS.sh solve               # serial; -p <n> threads, -m <n> MPI ranks, before 'solve'
python3 verify.py
```

`FUSS.sh solve` copies the repository's `bin/FUSS` into the case when it is newer.
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
   `fuss_add_test(<Name> <dir> <timeout> "<labels>" --verify)`.
5. Add a page to `docs/vv/` (problem, result, assessment) and to the `nav` of `mkdocs.yml`.
