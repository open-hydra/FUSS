# Verification & Validation

This section documents the Verification & Validation (V&V) test suite for FUSS. The purpose is to demonstrate that the solver produces correct results across a range of thermal problems, by comparison against analytical or well-established reference solutions, and that the numerical features perform as expected in terms of convergence and parallel efficiency.

## Test suite

| Test | Dim | Regime | Physics | Verification | Ref |
|---|---|---|---|---|---|
| [Radiative Rod](1D-steady-state-radiative-rod.md) | 1D | Steady-state | Conduction + radiation, nonlinear BC | Analytical (nonlinear algebraic eq.) | NAFEMS Thermal Test 9(iii) |
| [Oscillating Rod](1D-transient-oscillating-rod.md) | 1D | Transient | Oscillating wall temperature | Analytical Fourier series | NAFEMS Thermal Test 5 |
| [Temperature Plate](2D-steady-state-temperature-plate.md) | 2D | Steady-state | Conduction, prescribed-temperature BCs | NAFEMS benchmark value | NAFEMS Thermal Test 9(I) – Study 1 |
| [Volumetric Heat Plate](2D-steady-state-qvol-plate.md) | 2D | Steady-state | Conduction with volumetric heat source | NAFEMS benchmark value | NAFEMS Thermal Test 9(I) – Study 2 |
| [Convective Plate](2D-steady-state-convective-plate.md) | 2D | Steady-state | Mixed Dirichlet / convective BCs | NAFEMS benchmark value (T = 18.25 °C at A) | NAFEMS Thermal Test 9(I) |
| [Multi-Material Plate](2D-transient-multimat-plate.md) | 2D | Transient | Three-material plate, 1-D reduction | Analytical Fourier series, 6 snapshots | LANL multi-material diffusion benchmark |
| [Residual Convergence](residual-convergence.md) | 2D | Steady-state | Acceleration strategy comparison (RK3, IRS, multigrid) | Residual history vs. iteration count | Internal study |
| [Parallel Scaling](time-convergence.md) | 2D | Steady-state | OpenMP and MPI×OpenMP wall-clock scaling | Execution time vs. core count | Internal study |

## Running the tests

Each test case lives under `test/` and contains a `verify.py` script. The script reads the solver output, compares it against the reference solution, and writes a figure to the local `RESULTS/` directory.

```bash
cd test/<category>/<TestName>
./FUSS.sh solve
python3 verify.py          # check errors and export figure
python3 verify.py --plot   # as above, also display figure interactively
```

The numerics-feature cases (residual convergence and parallel scaling) read pre-collected output directly from `SOLUTION/` and do not require a prior solver run to generate intermediate files:

```bash
cd test/numerics-features/residual_convergence
python3 verify.py     # plots residual histories → RESULTS/residual_comparison.svg

cd test/numerics-features/time_convergence
python3 verify.py     # plots execution times → RESULTS/timing_comparison.svg
```

The script exits with code 0 on success and code 1 if any pointwise error exceeds the case-specific acceptance threshold.
