# Running FUSS

A FUSS case is a directory holding `input.ini` (the solver settings) and `INPUT/` (mesh, initial state, materials, boundary conditions). FUSS runs in that directory and writes its results to `OUTPUT/`.

FUSS reads only its own sections of `input.ini`, `[FUSS-*]` and the probe sections (see [Input File](input.md)). The other sections in the test cases (`[GRIB-*]`, `[GPB-*]`, `[ICB-*]`, `[STB-*]`, `[BCB-*]`, the named BC sections) are the input of two pre-processors of hydra, which are not part of FUSS:

- **GRIB** builds a structured multi-block mesh from `[GRIB-meshgen]` and `[GRIB-Block<n>]`.
- **ATLAS** writes the `INPUT/` files from the mesh and the other sections: `phase.txt` and `properties.dat` (GPB), `ic.tec` (ICB), `st.tec` (STB), `bc.txt` (BCB). ATLAS puts them in `fromATLAStoSolver/`; FUSS reads them from `INPUT/`.

Any script can write these files instead; `test/steady-state/skewed_Plate/make_inputs.py` does so for its case. The formats are given in [Initial Conditions](initial-conditions.md) and [Boundary Conditions](boundary-conditions.md).

---

## Case Directory

```
my_case/
├── input.ini            ← solver settings (plus pre-processor sections)
├── FUSS.sh              ← run script of the test cases (optional)
├── Tw_time.dat          ← time tables named in bc.txt (optional)
├── INPUT/
│   ├── ic.tec           ← mesh, initial temperature and material of every cell
│   ├── st.tec           ← volumetric heat source (optional)
│   ├── phase.txt        ← material names
│   ├── properties.dat   ← property tables
│   └── bc.txt           ← boundary conditions (bc2.txt, bc3.txt, ... for coarse multigrid levels)
├── MESH/
│   └── mesh.tec         ← grid written by GRIB in the test cases (not read by FUSS)
└── OUTPUT/              ← must exist before FUSS starts
```

`ic.tec` and `st.tec` take the extension of `ic-format` (`.tec`, `.szplt` or `.vtm`).

---

## Mesh

FUSS takes the grid from the initial-condition file `INPUT/ic.*` (on a restart, from the solution file it restarts from): the node coordinates of every block, one zone per block. Blocks are numbered in the order of the zones, and `bc.txt` refers to them by that number. The file format is described in [Initial Conditions](initial-conditions.md#initial-field-inputic).

The dimension of the problem and the type of the mesh are taken from the nodes of every block:

| Cells | Problem |
|-------|---------|
| $n_k > 1$ in any block | 3-D |
| $n_k = 1$ in every block, $n_j > 1$ in any | 2-D plane: layers one cell thick in $z$; 2-D axisymmetric when the two $k$-planes of every block are a wedge rotated about the $x$ axis |
| $n_j = n_k = 1$ in every block | 1-D |

Until October 2026 block 1 alone decided: a 3-D mesh whose block 1 is one cell thick in $z$ ran as 2-D, with block 1 decoupled from its $k$ neighbours.

On a 2-D mesh the plane or wedge decision looks at every node of every block: the $k = 1$ node plane of a block must be its $k = 0$ plane translated (plane mesh), or turned about the $x$ axis by the same angle in every block (wedge), to within 1 % of the thickness of the block. The coordinates therefore need enough significant digits (a wedge of 5 degrees written with five digits does not tell its type). A block nearer a turn than a translation, or nearer a translation than a turn, but not to within that tolerance, wedges of different angles, and wedge and plane blocks together stop the run at set-up with the block and the misfit. A 3-D mesh has no angle.

On a 2-D mesh faces 5 and 6 of every block must be zero-flux faces (code `0` or `300`, see [Boundary Conditions](boundary-conditions.md#face-numbering)).

---

## Running the Solver

Run FUSS from the case directory:

```bash
mkdir -p OUTPUT
ulimit -s unlimited
export OMP_NUM_THREADS=4
/path/to/FUSS/bin/FUSS
```

FUSS does not create `OUTPUT/`. It may need a large stack: `ulimit -s unlimited`, and `export KMP_STACKSIZE=100M` with the Intel runtime.

The test cases carry `FUSS.sh`; its options go before the task:

```bash
./FUSS.sh solve                  # 1 process, 1 thread
./FUSS.sh -p 4 solve             # 4 OpenMP threads
./FUSS.sh -m 2 solve             # 2 MPI ranks
./FUSS.sh -m 2 -p 10 solve       # 2 ranks x 10 threads
./FUSS.sh -b solve               # in the background: output to logfile and errors_file, PID in .ID
./FUSS.sh kill                   # stop the background run
```

`solve` creates `OUTPUT/` and `bin/`, sets the stack limits and `OMP_NUM_THREADS`, copies the repository's `bin/FUSS` to `./bin/FUSS` when it is newer, and runs it (through `mpirun` with `-m`). The repository is found through `MASTERDIR=../../..`, the depth of the test cases: change it for a case elsewhere. Threads, ranks and the block partition are described in [Parallel Execution](parallel.md).

---

## Log

FUSS prints the loading checks, a summary of the set-up, and one line per iteration (every `shell-diter` iterations). A failed check prints the reason and stops the run. A log of `test/steady-state/temperature_Plate`:

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
   Mesh type                      2-D, plane
   Time-step skew factor (max)    1.000000

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

FUSS | Iter =        1 | Global iter =        1 | Energy  residual = 0.357298E+04
...
FUSS | Iter =     3051 | Global iter =     3051 | Energy  residual = 0.356132E-04

 =========================================================================================
 FUSS timing
 =========================================================================================
   Iterations                     3051
   Solver                          5.86503E+00 s
   Solver per iteration            1.92233E-03 s
   Elapsed                         5.91949E+00 s
 =========================================================================================
```

`Mesh type` is 1-D, 2-D plane, 2-D axisymmetric wedge of an angle, or 3-D; `Time-step skew factor (max)` is the largest factor of the explicit step over the cells (1.000000 unless 3-D cells are skewed in all three pairs of directions). `Boundary faces` counts every record of `bc.txt`, zero-flux faces included; `Viscous wall` counts the wall records (codes 301–305). Time-accurate runs print `Time` and `Delta t` instead of the residual, and multigrid runs prefix the line with the grid level.

---

## Output

| File | Content |
|------|---------|
| `field.*` | Cell values: `T`, `matID`, `qvol`, `k`, `rho`, `cs` ($c_p$), `h` |
| `field<n>.*` | Numbered snapshots, with `sol-overwrite = false` |
| `wall.*` | Wall faces (codes 301–305): wall temperature `Tw` and heat flux into the solid `qw`; written when `bc.txt` has walls |
| `diagnostic.*` | Steady runs: cell residual `E` and local time step `dt` |
| `residual-history.dat` | Iteration, time, residual (every `res-diter` iterations) |
| `<probe>.txt` | Probe time histories |
| `field-level<n>.*`, `wall-level<n>.*` | Multigrid: solution and wall data of coarse level `<n>`, written while it is computed and when the run moves to the next finer level |
| `field-prolongated.*`, `wall-prolongated.*` | Multigrid: the coarse solution prolongated to the finest level, written with the coarse-level output (see `prolongated-diter`) |

The final solution of a run is always written to `field.*`; the solution due at the last iteration (by `sol-diter` or `sol-dtime`) is written before it, with its number and its wall and diagnostic files, so a run that ends on an output iteration leaves the same files as one that goes on. The field and wall files carry the solution time: the physical time of a time-accurate run, minus the global iteration count of a steady one (a negative `SOLUTIONTIME`, no `STRANDID`, in Tecplot; the `TIME` field data in vtk). Every field file also stores the counters a restart needs: `CYCLE` (the global iteration count) and `FUSS_ITERATION`, `FUSS_OUTPUT_NUMBER`, `FUSS_GRID_LEVEL`, `FUSS_LAST_OUTPUT_TIME`, `FUSS_RES_REFERENCE` (the residual of the first iteration, the unit of `res-threshold`), as field data in the vtk block files, as `DATASETAUXDATA` records after the last zone of a Tecplot ASCII file (six lines) and as auxiliary data of a `.szplt`. A script that reads a `.tec` file must stop at the first line that is not a number. The `.vtm` lists its block files relative to its own folder (`vtk/fieldBlock1.vts`). The text formats (`tecplot ascii`, `vtk ascii`) write every real with 17 significant digits and `.szplt` stores 64-bit values, so a restart reproduces the field exactly in every format (see [Restart](initial-conditions.md#restart)).

### Formats

`sol-format` (output and restart) and `ic-format` (initial condition) in `[FUSS-IO]`:

| Value | Files |
|-------|-------|
| `tecplot ascii` | Tecplot ASCII, `.tec` (default) |
| `tecplot binary` | Tecplot binary, `.szplt` (needs TecIO) |
| `vtk ascii` | VTK ASCII: `field.vtm`, block files under `OUTPUT/vtk/` |
| `vtk binary` | VTK base64, same layout |
| `vtk raw` | VTK raw appended binary, same layout |

### Frequency

| Key in `[FUSS-IO]` | Meaning |
|--------------------|---------|
| `sol-diter` | Write the solution every N iterations |
| `sol-dtime` | Write the solution every N seconds of simulated time |
| `sol-overwrite` | `true`: overwrite `field.*`; `false`: write `field1.*`, `field2.*`, ... |

---

## Probes

```ini
[FUSS-Probes]
probe1 = p1
probe2 = p2

[p1]
index-position = 1 90 1 1
variables = T
dtime = 1.0e-1
```

| Key | Meaning |
|-----|---------|
| `index-position` | `block i j k` of the probe cell |
| `position` | `x y z`, used when `index-position` is not given: the probe takes the cell with the nearest centre |
| `variables` | Recorded variables: only `T`; another name stops the run at set-up |
| `dtime` | Write at the multiples of `dtime` seconds of simulated time (at most once per step, from the first multiple after the start or restart time) |
| `diter` | Write every `diter` iterations |

Each probe writes `OUTPUT/<name>.txt` (`OUTPUT/p1.txt` above), one line per record: the time (or the iteration, for `diter` records) followed by the variables. A new run replaces the probe files; a restart appends to them, and creates any that are missing (for example a probe added before the restart); the rows written after the solution a restart reads are dropped first. The time-based rows fall on the multiples of `dtime`, counted with an integer, so a restart on an exact multiple writes the rows of the uninterrupted run. Every row of `residual-history.dat` and of the probe files is flushed when it is written, so a run that is killed leaves every row it wrote; the files are closed before MPI ends. The last step of the run is always sampled when its time or iteration matches.

---

## Restart

Set `newrun = false` in `[FUSS-Parameters]`: FUSS restarts from `OUTPUT/field.*` (or, when it is missing, from the last `OUTPUT/field<n>.*`) in `sol-format`. The thresholds (`res-threshold`, `time-threshold`, `iter-threshold`) can be changed before restarting. What is read back is described in [Initial Conditions](initial-conditions.md#restart): a run split in two and restarted writes the files of the run in one go.

---

## Exit Codes

Every error that stops the run ends it with **exit code 1**, after writing the message and flushing the output; a run that ends well has exit code 0. The line `[RANK r] STOP: <message>` names the stop (`STOP: <message>` in a serial build), and an `[ERROR]` line above it says where: the file, the line or the key and the value. With MPI every rank stops together. The errors are those of the input file, the material tables, `bc.txt` (the file, the line, the record and the type), the time tables, the initial condition or the restart solution (including a value that is not finite, and a solution without a time), the probes, the FAS and multigrid settings, a decomposition map, and the state check of the time step (a temperature that is NaN or outside the properties table, found inside the OpenMP loop: the run ends with exit code 1 on that rank and the launcher ends the others). No solution is written by a run that stops at set-up.

---

## Using FUSS as a Library

The build also produces the static library `libFUSSL.a`. `src/app/main.f90` is the reference driver:

```fortran
use FUSS_Advanced_Types_m, only: FUSS_simulation_type
use FUSS_Config_Types_m,   only: obj_sim_param
use FUSS_Procedures_m,     only: FUSS_type

type(FUSS_type)            :: FUSS
type(FUSS_simulation_type) :: simulation

call FUSS%setup(simulation)
obj_sim_param%TODO = 1
do while (obj_sim_param%TODO <= 2)
  call FUSS%solve(simulation, my_hook)
  if (obj_sim_param%TODO <= 2) call FUSS%postprocess(simulation)
end do
call FUSS%postprocess(simulation)
```

`solve` advances one iteration and calls `my_hook(domain)` at every Runge–Kutta stage, after the boundary fluxes and inside the OpenMP parallel region, with the domain of the level being computed (`main.f90` passes an empty subroutine; hydra adds its coupling terms there). With MPI, `main.f90` also calls `mpi_init_env` and `mpi_finalize_env` (module `FUSS_Mod_MPI`). Linking is described in [Installation](../getting-started/installation.md#library-linking-advanced).
