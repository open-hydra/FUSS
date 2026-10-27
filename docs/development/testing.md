# Testing

The test cases in `test/` run FUSS on problems with analytical or benchmark solutions and check numerical features. CTest runs them all; each case can also be run by hand.

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
| `steady-state/multilayer_Slab` | Three materials in series inside one block, unequal cells, exact solution | `verify.py` |
| `transient/oscillating_Rod` | Oscillating wall temperature (NAFEMS T5), probes | `verify.py` |
| `transient/multimat_Plate` | Three material layers (two materials), Fourier-series solution (LANL) | `verify.py` |
| `transient/variable_cp_Plate` | Temperature-dependent heat capacity, exact solution | `verify.py` |
| `numerics-features/chimera` | Chimera (overset) interface | smoke run |
| `numerics-features/chimera_mg` | Chimera donors on 3 multigrid levels, FAS cycle | smoke run |
| `numerics-features/anisotropic_Plate` | Cells of aspect ratio 100, 3 levels, W + IRS | smoke run |
| `numerics-features/multimat_liner` | Material jumps inside a block (four-layer liner in 2 blocks) at vnn 1.0 | smoke run |
| `numerics-features/residual_convergence/*` | Six convergence strategies (IRS, multigrid) | smoke runs; `plot.py` |
| `numerics-features/time_convergence/simdir` | Parallel timing study | smoke run; `plot.py` |

Each case directory holds `input.ini`, `MESH/mesh.tec`, `INPUT/` (`ic.tec`, `bc.txt`, `phase.txt`, `properties.dat`, optionally `st.tec`) and `FUSS.sh`. A `verify.py` reads the run's `OUTPUT/`, prints `Result: PASS` or `FAIL` and exits with 1 on failure. The two studies keep their data in `SOLUTION/`: `plot.py` draws the figures of the [V&V pages](../vv/index.md) from it. `numerics-features/convergence_figure/` is not a case: its `run_configs.sh` and `plot_convergence.py` run six configurations (none, IRS, V, W, V + IRS, W + IRS) on four of the cases and draw the figures of [Multigrid Cycles (FAS)](../vv/fas-cycles.md).

## CTest

`test/CMakeLists.txt` registers every case; `test/run_case.sh` runs one in a copy under the build tree (`build/test/runs/`), with one thread, so the source tree is never written.

```bash
ctest --test-dir build -L fast        # 9 verification cases, about a minute
ctest --test-dir build -j 4           # everything, a few minutes
ctest --test-dir build -R Qvol -V     # one test, with its output
```

| Label | Tests | What passes |
|---|---|---|
| `fast` | the 9 verification cases and `ChimeraEdge`, under a minute | run ends with status 0, `verify.py` passes; `ChimeraEdge`: a linear field on a sheared chimera pair is preserved to 1e-9 K, also at the cells at the ends of the chimera faces |
| `validation` | the 10 verification cases (also `convective_Plate`) and `ChimeraEdge` | same |
| `smoke` | the cases without a verifier: chimera, chimera_mg, anisotropic_Plate, multimat_liner, residual convergence, time convergence, and the FAS-keys-inert check | 300 iterations (100 for the timing case) end with status 0 and a residual history |
| `needs-mpi` | chimera, chimera_mg, multigrid with IRS, FAS W + IRS, multimat_liner (MPI builds only), `MPIBitIdentity` and the MPI tests below | 2 ranks write the same files as 1 process, bit for bit; `MPIBitIdentity`: 1×4, 2×1, 4×1 and 2×2 (ranks × threads) write the same files as 1×1, on a box of 4 blocks and a chimera pair of overlapping blocks, steady with IRS and time-accurate (`test/mpi/bit-identity/`) |

The tests of `test/input-io/` run scripts on copies of the cases; each prints one `CASE <name> PASS|FAIL` line per check, and the test fails if one case does:

| Test | Label | What it checks |
|---|---|---|
| `InputChecks` | `fast` | values out of their rule or that do not read (`vnn = -1`, `NaN`, `sol-diter = 2.5`, `sol-overwrite = yes`, `levels = 2.0`), a probe index out of range, a gap in the probe numbers, an `input.ini` that is missing, empty or only comments: exit code 1 with the key and the value; a comment with `=` and an unread key |
| `RuntimeInput` | `fast` | `input.ini` changed during the run: with `ini-error = warn` one warning and a residual history byte for byte that of the run without the change, with `stop` exit code 1 at the reading; a reading that changes nothing leaves every file |
| `ExitCodes` | `fast` | `bc.txt` (record line, short or null property line, record count, missing file), time tables (decreasing, missing, empty, bad row, NaN), probes (block, index, variable), FAS with one level, the initial condition (missing file, NaN cell), a wall temperature that is NaN, a restart without a solution: exit code 1, the message, no solution |
| `Restart` | `fast` | the run split in two writes the files of the run in one go (`vtk raw` and Tecplot ASCII; multigrid; FAS; probes and `sol-dtime`; numbered solutions; splits that end on an output); a solution without counters; the rows after the solution; a transient run from a steady solution; `INPUT/st.*` read with `ic-format`; a vtk solution without `TIME` or with a NaN counter |
| `RestartBinary` | `fast`, needs TecIO | the same split with `tecplot binary` (`.szplt`) |
| `DecompositionMap` | `fast` | a map with every block on rank 0 is the run without a map; every error of the map stops the run (exit code 1); a map is read only when the key names it |
| `MPIRuntimeInput` | `needs-mpi` | two ranks in two folders: a missing, different or rewritten `input.ini` in the folder of rank 1 changes nothing |
| `MPIExitCodes` | `needs-mpi` | the STOP line of both ranks, exit code 1, for an input error, a `bc.txt` error, an initial condition that rank 1 alone cannot read or holds a NaN, `ini-error = stop` at a reading during the run |
| `MPIRestart` | `needs-mpi` | the split restart on 2 ranks (`nominal`, and `mg_irs` with two grid levels), 0 differing files |
| `MPIDecompositionMap` | `needs-mpi` | the two blocks reversed on 2 ranks, and on 4 ranks (two idle), write the files of the serial run; a map for 3 ranks run on 2 stops every rank |

`run_case.sh <case> <work dir> <binary> [--iters N] [--verify] [--mpi N]` caps the run at N iterations by rewriting `iter-threshold` and every `level<m>-iter` of the copy.

## Running a case by hand

```bash
cd test/steady-state/temperature_Plate
./FUSS.sh solve               # serial; -p <n> threads, -m <n> MPI ranks, before 'solve'
python3 verify.py
```

`FUSS.sh solve` copies `bin/FUSS` of the repository into the case when it is newer.

## Adding a case

1. Create the directory with `input.ini` and the mesh. The `[GRIB-*]`, `[GPB-*]`, `[ICB-*]` and `[BCB-*]` sections of `input.ini` describe how ATLAS generates `MESH/` and `INPUT/`; the first line of `phase.txt` must be `solid-bulk` for ATLAS. A case whose data ATLAS cannot express carries a script that writes them (`skewed_Plate/make_inputs.py`).
2. Write `verify.py`: read `OUTPUT/`, compare with the reference, print `Result: PASS` or `FAIL`, exit 1 on failure.
3. Register it in `test/CMakeLists.txt` with `fuss_add_test(<Name> <dir> <timeout> "<labels>" --verify)`.
4. Add a page to `docs/vv/` and to the `nav` of `mkdocs.yml`.
