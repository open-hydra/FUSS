# Verification & Validation

The test cases compare FUSS with analytical solutions and benchmark values, and measure how the
numerical features affect convergence and parallel speed. Every number on these pages is a
measurement of FUSS. Each page gives the problem, the result and its assessment; the setup and the
run commands of each case are in the `README.md` of its test folder.

## Test suite

| Test | Dim | Regime | Physics | Verification | Ref |
|---|---|---|---|---|---|
| [Radiative Rod](1D-steady-state-radiative-rod.md) | 1D | Steady | Conduction + radiation (BC 304) | Free-end temperature, 927 K | NAFEMS T9 (iii) |
| [Convective Rod with Imposed Flux](1D-steady-state-convective-qrad-rod.md) | 1D | Steady | Convection + prescribed flux (BC 303) | Closed-form linear solution | Internal |
| [Multilayer Slab](1D-steady-state-multilayer-slab.md) | 1D | Steady | Three materials in series inside one block, unequal cells | Piecewise-linear exact solution, 1.3 × 10⁻⁷ K | Internal |
| [Oscillating Rod](1D-transient-oscillating-rod.md) | 1D | Transient | Oscillating wall temperature | Peak temperature rise at x = 0.09 m, 94.69 K | NAFEMS T5 |
| [Skewed Plate](2D-steady-state-skewed-plate.md) | 2D | Steady | Walls 301, 302, 303, 305 on a 45° sheared mesh | Linear field, exact to 5.6 × 10⁻⁹ K | Internal |
| [Temperature Plate](2D-steady-state-temperature-plate.md) | 2D | Steady | Prescribed temperatures | Temperature at one point, 260.5 °C | NAFEMS T9 (i) |
| [Volumetric Heat Plate](2D-steady-state-qvol-plate.md) | 2D | Steady | Volumetric heat source | Temperature at one point, 310.1 °C | NAFEMS T9 (ii) |
| [Convective Plate](2D-steady-state-convective-plate.md) | 2D | Steady | Prescribed temperature, convection, insulated edge | Temperature at point E, 18.25 °C | NAFEMS T4 |
| [Multi-Material Plate](2D-transient-multimat-plate.md) | 2D | Transient | Three material layers (two materials) | Fourier-series solution, 6 snapshots | LANL multi-material diffusion |
| [Variable-c<sub>p</sub> Plate](2D-transient-variable-cp-plate.md) | 2D | Transient | Uniform heating, c<sub>p</sub>(T) linear | Closed-form solution | Internal |
| [Residual Convergence](residual-convergence.md) | 2D | Steady | IRS and multigrid against plain time stepping | Iterations to convergence | Internal study |
| [Multigrid Cycles (FAS)](fas-cycles.md) | 2D, 3D | Steady | FAS V / W cycles and IRS on a two-material plate, a four-layer liner, overlapping chimera blocks and anisotropic cells | Iterations and wall time to convergence, field equal to the single-grid one | Internal study |
| [Parallel Scaling](parallel-scaling.md) | 2D | Steady | OpenMP and MPI × OpenMP | Wall-clock time vs. cores | Internal study |

## Running the tests

Each verification case has a `verify.py`. It reads the run's `OUTPUT/`, prints the comparison,
writes its figure (if any), and exits with 0 on PASS and 1 on FAIL.

```bash
cd test/<category>/<case>
./FUSS.sh solve
python3 verify.py
```

CTest runs every case in a copy under the build tree (see [Testing](../development/testing.md)):

```bash
ctest --test-dir build -L fast    # the 9 short verification cases and the quick checks
ctest --test-dir build            # everything, with the smoke runs and the MPI checks
```
