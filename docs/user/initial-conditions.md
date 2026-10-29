# Initial Conditions

A new run (`newrun = true`) starts from four files in `INPUT/`. FUSS does not read the `[GPB-*]`, `[ICB-*]` or `[STB-*]` sections of `input.ini`: ATLAS turns them into these files.

| File | Content | ATLAS tool |
|------|---------|------------|
| `ic.*` | Grid, initial temperature and material of every cell | ICB, from `[ICB-Block<n>]` |
| `st.*` | Volumetric heat source of every cell (optional) | STB, from `[STB-Block<n>]` |
| `phase.txt` | Material names | GPB, from `[GPB-Phase1]` |
| `properties.dat` | Property table of every material | GPB, from `[GPB-Phase1]` |

A restart replaces `ic.*` by the last solution (see [Restart](#restart)).

---

## Initial Field: `INPUT/ic.*`

The extension follows `ic-format` in `[FUSS-IO]`: `ic.tec` (`tecplot ascii`, default), `ic.szplt` (`tecplot binary`) or `ic.vtm` (`vtk ascii`, `vtk raw`). The file is a structured multi-block field with one zone per block, in block order:

- the node coordinates `x y z`, which are the grid of the run (see [Mesh](using.md#mesh));
- two cell-centred variables, in this order: the temperature `T` (K) and the material index `matID`. FUSS takes them by position, not by name; a block with fewer than two variables stops the run at set-up.

`matID = n` selects the `n`-th material of `phase.txt`, i.e. the `n`-th zone of `properties.dat`. Every initial temperature must lie inside the table, $T_\text{min} \le T < T_\text{max}$.

The header of `ic.tec` in `test/transient/oscillating_Rod` (100 × 2 × 1 cells):

```
 VARIABLES ="x" "y" "z" T matID
 ZONE  T = B1-SP, I=101, J=3, K=2, DATAPACKING=BLOCK, VARLOCATION=([1-3]=NODAL,[4-5]=CELLCENTERED), SOLUTIONTIME=0.100000000000000E+002, STRANDID = 0
```

With ATLAS, `[ICB-Block<n>]` sets block `n`:

```ini
[ICB-Block1]
type     = homogeneous
material = steel
T        = 1500.0
```

`material` is a name of the `[GPB-Phase1]` list; ICB writes its position in that list as `matID`, and writes `matID = 1` for a name that is not in the list. Temperature profiles from a file and interpolation from an old solution are ICB options: see the ATLAS documentation.

---

## Heat Source: `INPUT/st.*`

An optional file in the same format and with the same grid as `ic.*` (FUSS stops if the block dimensions differ). Its first cell-centred variable is the volumetric heat source $\dot q_\text{vol}$ (W/m³) of the cell, added to the energy balance of the cell at every stage and constant in time. Without the file $\dot q_\text{vol} = 0$.

With ATLAS:

```ini
[STB-Block1]
qvol = 1.e6
```

---

## Materials

### `INPUT/phase.txt`

The first line is a header and is skipped. Each further line names one material in its first word, in the order of the zones of `properties.dat`; FUSS ignores the rest of the line. From `test/transient/multimat_Plate`:

```
solid-bulk
insulation 1
conductive 1
```

### `INPUT/properties.dat`

A Tecplot point table with one zone per material. From `test/steady-state/temperature_Plate`:

```
TITLE = "Mass Thermodynamic Properties"
VARIABLES = "Temperature", "Cp", "Density", "Conductivity", "Energy"
ZONE T="mat"
I=2000, F=POINT
1.0 460.000000 7850.000000 52.000000 3611000.000000
2.0 460.000000 7850.000000 52.000000 7222000.000000
...
```

| Column | Unit |
|--------|------|
| Temperature | K |
| Cp ($c_p$) | J/(kg·K) |
| Density ($\rho$) | kg/m³ |
| Conductivity ($k$) | W/(m·K) |
| Energy ($h = \int \rho\,c_p\,\mathrm{d}T$) | J/m³ |

- The rows are 1 K apart. FUSS takes the first temperature of the first zone as $T_\text{min}$ and the number of rows to set $T_\text{max} = T_\text{min} + I - 1$; it does not read the other temperatures. Every zone must cover the same range.
- Properties are interpolated linearly between the two integer temperatures around $T$.
- The Energy column is used up to an additive constant. At load FUSS warns if its increments differ from $\rho\,c_p$ by more than 1 %.
- At load every Cp, Density, Conductivity and Energy must be finite, Cp and Density positive and Conductivity not negative (0 gives a material that conducts no heat, as the insulation of `test/transient/multimat_Plate`); a value outside stops the run at set-up, naming the material, the row and the column, also in a row that no cell reaches.
- A cell temperature outside $[T_\text{min}, T_\text{max})$ (or not finite) stops the run with an error and exit code 1, also with IRS (see [Time Integration](../theory/time-integration.md#admissible-updates)).

The property model, the Energy column and a script for temperature-dependent tables are described in [Material Properties](../theory/thermo.md).

With ATLAS, GPB writes both files:

```ini
[GPB-Phase1]
type     = solid
material = steel aluminium
cp       = 460.0  900.0
rho      = 7850.0 2700.0
k        = 52.0   200.0
Tmax     = 5000
```

With constant `cp`, `rho`, `k` (one value per material, in the order of `material`) the tables run from `Tmin` (default 1 K) to `Tmax` (default 5000 K) with $h = \rho\,c_p\,T$. Leave out the GPB key `name`: it prefixes the file names (`<name>-phase.txt`, ...), while FUSS run on its own reads `INPUT/phase.txt` and `INPUT/properties.dat`. Temperature-dependent tables can also be written by a script, as `test/transient/variable_cp_Plate/make_properties.py` does.

---

## Multi-Material Domains

Give each material its own blocks and join them with block connections (code 101, `connection` in ATLAS). At a connection the interface temperature follows from the continuity of the heat flux, each side with its own conductivity over half a cell: a series resistance, which on a uniform orthogonal mesh is the harmonic mean of the two conductivities. Inside a block a face between two materials takes the series conductivity of the two half cells (see [Face Conductivity](../theory/numerics.md#face-conductivity)), at the mean temperature of the two cells for one material. Test case: [2D Transient Multi-Material Plate](../vv/2D-transient-multimat-plate.md).

---

## Restart

With `newrun = false` FUSS reads `OUTPUT/field.*` in `sol-format`, or, when it is missing, the last of `OUTPUT/field1.*`, `field2.*`, .... The final solution of a run is always written to `field.*`.

- The grid, `T` and `matID` come from that file; $\dot q_\text{vol}$ is read again from `INPUT/st.*`, so it can be changed between runs. `INPUT/st.*` is an input file: it has the extension of `ic-format` (`st.tec`, `st.szplt` or `st.vtm`), whatever `sol-format` is (it was read with `sol-format`, and a case with `st.tec` restarted in `vtk raw` ran with $\dot q_\text{vol} = 0$ in silence).
- **Kind of run.** A solution is the one of a steady run when its time is negative (the writers stamp a steady solution with minus the global iteration count), or, in a Tecplot ASCII file written before, when its zones say `STRANDID = 0`; otherwise it is time-accurate. A solution without a time (vtk without `TIME`, Tecplot ASCII without `SOLUTIONTIME`) cannot tell where to go on: the run stops at set-up and says to use the file as `INPUT/ic.*` with `newrun = true`. So does a time or a counter that is not finite.
- **Counters.** A solution of the same kind of run as the restart stores the counters of its run (see [Output](using.md#output)), and the restart restores them: the global iteration count, the iteration count of the finest level (the VNN ramp), the number of the last numbered output, the time of the last time-based output and the grid level. The run split in two then writes the files of the run in one go, byte for byte, in every format (`tecplot ascii`, `tecplot binary`, `vtk ascii`, `vtk binary`, `vtk raw`), with probes, time-based outputs, multigrid and FAS cycles (a run cut inside a FAS cycle restarts from the cut state: split on a cycle boundary).
- Time-accurate runs continue from the solution time of the file. A time-accurate run restarted from a steady solution starts at $t = 0$ with new counters (its first row is that of the run started from the same field as `ic.*`); a steady run restarted from a time-accurate solution counts its iterations from the integer part of the time.
- A solution without the counters (written by a version before them, or of the other kind of run) restores T, the time and, for a steady run, the iteration count: the counters start as in a new run, the multigrid starts from the coarsest level and the numbered outputs are numbered by scanning `OUTPUT/` (the legacy path below). Move the earlier snapshots away first if you do not want them kept.
- `time-threshold` is the simulated time, while `iter-threshold` and `level<m>-iter` count the iterations of the restarted run only.
- `residual-history.dat` and the probe files are appended to; the rows written after the solution read (a run that went on past a numbered solution) are dropped first, a row cut in the middle too, and the files end with a newline.
- With `sol-overwrite = false` the numbering continues from `FUSS_OUTPUT_NUMBER` of the solution: the snapshots of the previous run are kept. For a solution without counters it continues after the last `field<n>.*` found in `OUTPUT/` (legacy path).
- With multigrid (`levels > 1`) a solution with counters, written at level 1 (`field.*`, `field<n>.*`), goes on at level 1: the coarse levels are restricted from it and the run does not start again from the coarsest level. A solution without counters starts from the coarsest level, as before.

**Limits.** The residual threshold `res-threshold` is absolute, so no reference value is stored. With `res-diter > 1` the first iteration of a restart computes a fresh residual where the uninterrupted run tests the one of its last multiple, so an end by `res-threshold` can come earlier within that window. With `vnn-rise-threshold > 0` and FAS cycles the counters of the coarse levels' ramp are not stored. A `.szplt` written before the signed time (ORION before a978930) reads as time-accurate at $t = |{\rm iteration}|$: a steady restart from it counts right, a time-accurate restart from it starts at that time.
