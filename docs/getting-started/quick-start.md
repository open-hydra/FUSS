# Quick Start: 2-D Heated Plate

This page runs the bundled case `test/steady-state/temperature_Plate`: steady conduction in a plate with prescribed wall temperatures, checked against the NAFEMS thermal benchmark T9 (i). The [V&V page](../vv/2D-steady-state-temperature-plate.md) documents the case in full.

!!! note "Prerequisites"
    - `bin/FUSS` built (see [Installation](installation.md)).
    - Python 3 with NumPy and Matplotlib, for `verify.py`.

---

## Problem

A 0.6 m × 0.4 m plate ($k = 52$ W/(m·K), $\rho = 7850$ kg/m³, $c_p = 460$ J/(kg·K)) is held at 1273.15 K on its left edge and at 273.15 K on the other three. The steady temperature satisfies

$$\frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2} = 0$$

| Parameter | Value |
|-----------|-------|
| Grid | 60 × 40 cells, one block |
| Left edge (face 1) | $T = 1273.15$ K (BC 302) |
| Right, bottom, top edges (faces 2, 3, 4) | $T = 273.15$ K (BC 302) |
| Faces 5, 6 | `null` (2-D case) |
| Initial temperature | 273.15 K |
| Reference | $T = 260.5$ °C at point A (x = 0.2 m, y = 0.2 m) |

---

## Case directory

```
test/steady-state/temperature_Plate/
├── input.ini           # FUSS settings (and the pre-processing sections)
├── FUSS.sh             # run script
├── verify.py           # comparison with the NAFEMS value
├── INPUT/
│   ├── ic.tec          # grid and initial field (T, matID)
│   ├── bc.txt          # boundary conditions
│   ├── phase.txt       # material names
│   └── properties.dat  # property table (1 K steps)
├── MESH/mesh.tec       # grid written by the mesh generator (not read by FUSS)
└── SOLUTION/           # a stored field.tec
```

FUSS reads the grid from `INPUT/ic.tec`. The `INPUT/` files were written by ATLAS (and `MESH/` by the GRIB mesh generator) from the `[GRIB-*]`, `[GPB-*]`, `[ICB-*]` and `[BCB-*]` sections of `input.ini`; they ship with the case, so neither tool is needed to run it.

---

## Solver configuration

`input.ini`:

```ini
[FUSS-Parameters]
newrun = true
res-threshold = 1e-8

[FUSS-Numerics]
time-scheme = RK3
vnn = 1.0
time-accurate = false
integration-variables = cons
irs-beta = 0.5

[FUSS-IO]
sol-diter     = 10000
sol-overwrite = true
sol-format    = tecplot ascii
ini-diter = 1000

;

[GRIB-meshgen]
outpath = MESH/
method = fmsh
type = 2Dplane
width = 0.01

[GRIB-Block1]
pivot = 0.0 0.0 0.0
length = 0.6
height = 0.4
nx = 60
ny = 40
strx = U 1.0
stry = U 1.0

;

[GPB-Phase1]
type = solid
material = mat
cp  = 460.0
rho = 7850.0
k   = 52.0
Tmax = 2000

;

[ICB-Block1]
type = homogeneous
material = mat
T = 273.15

;

[BCB-Block1]
face1 = Twall-hot
face2 = Twall-cold
face3 = Twall-cold
face4 = Twall-cold
face5 = null
face6 = null

[Twall-hot]
type = wall
T = 1273.15

[Twall-cold]
type = wall
T = 273.15
```

| Key | Meaning here |
|-----|--------------|
| `newrun = true` | Start from `INPUT/ic.tec` (`false` restarts from `OUTPUT/field.tec`) |
| `res-threshold = 1e-8` | Stop when the residual (L2 norm over all cells of the temperature change of one iteration, K) falls below $10^{-8}$ |
| `time-accurate = false` | Steady run: local time step in each cell |
| `time-scheme = RK3`, `vnn = 1.0` | SSP RK3; `vnn` above the explicit limit (0.628 for RK3), allowed here by IRS |
| `irs-beta = 0.5` | Implicit residual smoothing |
| `integration-variables = cons` | Advance the volumetric enthalpy $h$ |
| `sol-diter`, `sol-overwrite`, `sol-format` | Write `OUTPUT/field.tec` every 10000 iterations and at the end, overwriting it |
| `ini-diter = 1000` | Re-read the runtime keys of `input.ini` every 1000 iterations |

All keys: [Input Parameters](../user/registry.md). Comment lines in `input.ini` start with `!`, `;` or `#`; a comment line may contain `=` (FiNeR 2.1.0 or later).

---

## Running

```bash
cd test/steady-state/temperature_Plate
./FUSS.sh solve
```

`FUSS.sh solve` creates `OUTPUT/` and `bin/`, copies the repository's `bin/FUSS` into `./bin/` when it is newer, sets `ulimit -s unlimited`, `KMP_STACKSIZE=100M` and `OMP_NUM_THREADS`, and runs `./bin/FUSS`. Options go before `solve`:

```bash
./FUSS.sh -p 4 solve         # 4 OpenMP threads
./FUSS.sh -m 2 solve         # 2 MPI ranks (mpirun -np 2 --map-by socket --bind-to socket)
./FUSS.sh -m 2 -p 4 solve    # 2 ranks x 4 threads
./FUSS.sh -b solve           # background: output to logfile and errors_file, PID in .ID
./FUSS.sh kill               # stop the background run
```

`./FUSS.sh compile` runs `install.sh compile` in the repository and copies the new `bin/FUSS`. Running with several threads or ranks needs a build with OpenMP or MPI; see [Parallel Execution](../user/parallel.md).

The run takes 4894 iterations, under a minute on one thread.

---

## Screen output

After a header that reports the threads and ranks of the build, FUSS checks its input files and prints the set-up:

```
 =========================================================================================
 Loading
 =========================================================================================
   Properties                     OK
   Input file                     OK
   Initial conditions             OK
   Boundary conditions            OK
 =========================================================================================

 =========================================================================================
 Set-up
 =========================================================================================
 Domain
   Blocks                         1
   Cells                          2400
   Boundary faces                 5000

 Boundary conditions
   Viscous wall                   200

 Input/Output
   Initial conditions file        INPUT/ic.tec
   Solution format                tecplot ascii
   Probes number                  0
 =========================================================================================

 =========================================================================================
 Physical model
 =========================================================================================
 Solid model
 Equation                         Fourier
 =========================================================================================

 =========================================================================================
 Numerical scheme
 =========================================================================================
 Space
   Diffusive face metric          stencil

 Time
   Scheme                         Third-order Runge-Kutta
   Integration variables          Conservative
   Implicit residual smoothing    Beta set to 0.500000E+00
 =========================================================================================
```

A file that cannot be read or fails its checks is reported as `FAIL` with a message, and the run stops. "Boundary faces" counts every boundary cell face, including the 4800 `null` faces 5 and 6; "Viscous wall" counts the wall faces (BC 301–305).

One line per iteration follows (`shell-diter = 1` by default); an excerpt:

```
FUSS | Iter =        1 | Global iter =        1 | Energy  residual = 0.357298E+04
FUSS | Iter =     1000 | Global iter =     1000 | Energy  residual = 0.320304E+00
FUSS | Iter =     2000 | Global iter =     2000 | Energy  residual = 0.378109E-02
FUSS | Iter =     3000 | Global iter =     3000 | Energy  residual = 0.446604E-04
FUSS | Iter =     4000 | Global iter =     4000 | Energy  residual = 0.527508E-06
FUSS | Iter =     4894 | Global iter =     4894 | Energy  residual = 0.997409E-08
```

The run ends with a timing summary (iterations, solver time, time per iteration, elapsed time).

### Output files

| File | Content |
|------|---------|
| `OUTPUT/field.tec` | Grid and cell fields: `T`, `matID`, `qvol`, `k`, `rho`, `cs`, `h` |
| `OUTPUT/wall.tec` | One zone per wall face (`B1F1` … `B1F4`): wall temperature `Tw` and heat flux into the solid `qw` |
| `OUTPUT/residual-history.dat` | Iteration, time (−1 in steady runs), residual |
| `OUTPUT/diagnostic.tec` | Cell residual `E` and local time step `dt` |

See [Running FUSS](../user/using.md) for all output options.

---

## Verification

```bash
python3 verify.py
```

`verify.py` reads `OUTPUT/field.tec`, interpolates the cell-centre temperatures bilinearly at point A, compares the result with 533.65 K (260.5 °C), and writes the figure `2D-steady-state-temperature-plate-field.svg`. It prints `Result: PASS` and exits with status 0 when the relative error (in kelvin) is below 1 %, and exits with status 1 otherwise. This run gives 260.509 °C, a relative error of 0.0018 %.

---

## Next steps

- Other cases: `test/steady-state/`, `test/transient/`, `test/numerics-features/`; see [Testing](../development/testing.md) and [Verification & Validation](../vv/index.md).
- Setting up a case: [User Guide](../user/index.md) and [Input File](../user/input.md).
- The code: [Code Structure](../development/structure.md) and [Contributing](../development/contributing.md).

---

## References

- NAFEMS Publication P16, *Benchmark Tests for Thermal Analysis*, Test 9 (i), YR3087, Vol. 2, 1986.
