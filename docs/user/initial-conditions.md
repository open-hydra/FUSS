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
- two cell-centred variables, in this order: the temperature `T` (K) and the material index `matID`. FUSS takes them by position, not by name.

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
- A cell temperature outside $[T_\text{min}, T_\text{max})$ stops the run with an error; with IRS the smoothed update is replaced first (see [Time Integration](../theory/time-integration.md#admissible-updates)).

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

Give each material its own blocks and join them with block connections (code 101, `connection` in ATLAS). At a connection the interface temperature follows from the continuity of the heat flux, each side with its own conductivity over half a cell: a series resistance, which on a uniform orthogonal mesh is the harmonic mean of the two conductivities. Inside a block the flux through a face uses the conductivity of the lower-index cell at the mean temperature of the two cells, so a material change inside a block is not treated as a series resistance. Test case: [2D Transient Multi-Material Plate](../vv/2D-transient-multimat-plate.md).

---

## Restart

With `newrun = false` FUSS reads `OUTPUT/field.*` in `sol-format`, or, when it is missing, the last of `OUTPUT/field1.*`, `field2.*`, .... The final solution of a run is always written to `field.*`.

- The grid, `T` and `matID` come from that file; $\dot q_\text{vol}$ is read again from `INPUT/st.*` (with the extension of `sol-format`), so it can be changed between runs.
- Time-accurate runs continue from the solution time of the file, steady runs from its iteration count (a `.szplt` file gives 0, see [Output](using.md#output)).
- `time-threshold` is the simulated time, while `iter-threshold` counts the iterations of the restarted run only.
- `residual-history.dat` and the probe files are appended to.
- With `sol-overwrite = false` the numbering continues after the last `field<n>.*` found in `OUTPUT/`: the snapshots of the previous run are kept.
- With multigrid (`levels > 1`) the run starts again from the coarsest level, restricted from the restart file, and the prolongation to the finer levels replaces the restart solution. To continue on the finest grid, restart with `levels = 1`.
