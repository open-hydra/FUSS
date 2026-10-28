# Material Properties

FUSS models heat conduction in solid media.  Every cell carries a **material index** that selects the thermophysical properties required by the heat equation: density $\rho$, specific heat capacity $c_p$, and thermal conductivity $k$.  All three can vary with temperature; they are read from a table, one per material, and interpolated linearly in temperature.

---

## Thermophysical Properties

The three properties that appear in the heat equation are:

| Property | Symbol | Unit | Role |
|----------|--------|------|------|
| Thermal conductivity | $k(T)$ | W/(m·K) | Controls the rate of spatial heat diffusion |
| Density | $\rho(T)$ | kg/m³ | Inertia of the thermal response |
| Specific heat capacity | $c_p(T)$ | J/(kg·K) | Energy stored per unit mass per degree |

Together they define the **local thermal diffusivity**

$$
\alpha(T) = \frac{k(T)}{\rho(T)\,c_p(T)} \quad \text{[m²/s]}
$$

The time step of each cell uses $\alpha$ at the cell temperature; time-accurate runs advance all cells with the smallest of these steps (see [Time Integration](time-integration.md#diffusive-stability-condition-vnn)).

---

## Property Table: `properties.dat`

The solver reads the properties only from the **Tecplot ASCII table**

```
INPUT/properties.dat
```

and the material names from `INPUT/phase.txt`.  The table holds one ZONE per material, each a list of rows at consecutive integer temperatures in kelvin.

### File format

```
TITLE = "Mass Thermodynamic Properties"
VARIABLES = "Temperature", "Cp", "Density", "Conductivity", "Energy"
ZONE T="<material_name>"
I=<N>, F=POINT
<T_min>.0  <cp_1>  <rho_1>  <k_1>  <h_1>
<T_min+1>.0  <cp_2>  <rho_2>  <k_2>  <h_2>
...
```

| Column | Symbol | Unit | Notes |
|--------|--------|------|-------|
| Temperature | $T$ | K | $T_\text{min}, T_\text{min}+1, \ldots, T_\text{min}+N-1$ |
| Cp | $c_p(T)$ | J/(kg·K) | Specific heat at constant pressure |
| Density | $\rho(T)$ | kg/m³ | Material density |
| Conductivity | $k(T)$ | W/(m·K) | Thermal conductivity |
| Energy | $h(T)$ | J/m³ | Volumetric enthalpy $\int \rho\,c_p\,\mathrm{d}T$ (see below) |

The columns are read by position.  Of the Temperature column only the first value is used: FUSS sets $T_\text{min}$ to it (rounded to an integer) and $T_\text{max} = T_\text{min} + N - 1$ from the row count $N$ of the first zone, and assumes the rows to be 1 K apart.  Every zone must therefore cover the same temperatures as the first.  ATLAS GPB writes the rows from `Tmin` (default 1) to `Tmax` (default 5000) of the `[GPB-Phase1]` section.

### Multiple materials

When the simulation uses more than one material, each material occupies its own ZONE.  The material index of a cell is the position of its ZONE in the file:

```
ZONE T="insulation"
I=200, F=POINT
1.0  1.0  1.0  1.e-12  1.0
...
ZONE T="conductive"
I=200, F=POINT
1.0  1.0  1.0  1.0  1.0
...
```

`INPUT/phase.txt` lists the materials in the same order.  FUSS skips its first line and reads the first word of each following line as a material name, used in messages only; the number of these lines is the number of materials.  ATLAS writes the phase type on the first line and a group number after each name:

```
solid-bulk
insulation 1
conductive 1
```

---

## Property Interpolation

The routines `co_Kappa`, `co_Rho`, `co_Cp` and `co_H` interpolate linearly between the two rows that bracket the temperature:

$$
k(T) = k(\lfloor T \rfloor)
       + \bigl(T - \lfloor T \rfloor\bigr)\,
         \bigl[k(\lfloor T \rfloor + 1) - k(\lfloor T \rfloor)\bigr]
$$

and analogously for $\rho$, $c_p$ and $h$.  The rows are addressed directly by the integer part of $T$, so a lookup costs $O(1)$.  The inverse, $T(h)$ in `co_T`, finds the first row above $h$ and interpolates linearly in that interval; below the first row or above the last it extrapolates the end interval.  The row is found by bisection, $O(\log n)$ in the number of rows, when the material's Energy column does not decrease anywhere (always, for a physical table), and by a scan from $T_\text{min}$ otherwise; on a non-decreasing column the two find the same row.

Where the properties are evaluated:

| Property | Used for | Temperature |
|---|---|---|
| $k$ | interior faces | mean of the two cell temperatures; between two materials the series conductivity of the two half cells |
| $k$ | block-connection, chimera and wall faces | temperature of each cell |
| $\rho\,c_p$ | `prim` update | cell temperature |
| $\alpha = k/(\rho\,c_p)$ | time step | cell temperature |
| $h$, $T(h)$ | `cons` update, chimera ghost cells, multigrid restriction | cell temperature |

!!! warning "Temperature range"

    A lookup reads the rows $\lfloor T \rfloor$ and $\lfloor T \rfloor + 1$, so it is valid for $T_\text{min} \le T < T_\text{max}$.  After every Runge–Kutta stage each cell is checked, with or without IRS: a temperature outside this range, or not finite, stops the run with an error naming the cell and exit code 1 (see [Admissible updates](time-integration.md#admissible-updates); until October 2026 a smoothed update was halved instead).  Temperatures are not clamped.  The table must therefore cover every temperature the solution can reach, prescribed wall temperatures included.

---

## Energy Column and the Enthalpy Inversion

The **Energy** column is the volumetric enthalpy

$$
h(T) = \int^{T} \rho(T')\,c_p(T')\,\mathrm{d}T',
$$

up to an additive constant, which plays no role.  With `integration-variables = cons` (the default) every Runge–Kutta stage advances $h$ and recovers $T$ by inverting the table; with `prim` the stage advances $T$ dividing by $\mathrm{d}h/\mathrm{d}T = \rho\,c_p$, the same equation (see [Governing Equations](governing-equations.md#state-vector)).  When it loads the table FUSS compares each increment $h(T+1) - h(T)$ with the mean of $\rho\,c_p$ at $T$ and $T+1$ and warns if they differ by more than 1 % anywhere: a column holding $\rho\,c_p\,T$ for a temperature-dependent material would give the conservative path the heat capacity $\rho\,c_p + T\,\mathrm{d}(\rho\,c_p)/\mathrm{d}T$.  For constant properties the column reduces to $h(T) = \rho\,c_p\,T$, as in the single-material test cases, whose first row reads

$$
h = 7850\;\text{kg/m}^3 \times 460\;\text{J/(kg·K)} \times 1\;\text{K} = 3\,611\,000\;\text{J/m}^3
$$

For **variable properties** each row must accumulate the product $\rho(T)\,c_p(T)$ up to that temperature, integrating the actual variation.

---

## Variable Properties: How to Build the Table

To use temperature-dependent properties, provide a `properties.dat` file in which the Cp, Density, and Conductivity columns vary from row to row.  A Python snippet that writes the table for a material with known $k(T)$, $\rho(T)$, $c_p(T)$, accumulating $h$ with the trapezoidal rule that the load-time check uses:

```python
import numpy as np

T_min = 1         # first integer temperature [K]
T_max = 2000      # last  integer temperature [K]
T_arr = np.arange(T_min, T_max + 1, dtype=float)

# Define property functions (replace with real data or polynomial fits)
def k_func(T):   return 52.0 - 0.01 * (T - 300)          # W/(m·K)
def rho_func(T): return 7850.0 - 0.05 * (T - 300)        # kg/m³
def cp_func(T):  return 460.0  + 0.20 * (T - 300)        # J/(kg·K)

k_arr   = k_func(T_arr)
rho_arr = rho_func(T_arr)
cp_arr  = cp_func(T_arr)

# Volumetric enthalpy: h(T_min) = rho*cp*T_min, then trapezoidal increments over 1 K
c     = rho_arr * cp_arr
h_arr = c[0] * T_min + np.concatenate(([0.0], np.cumsum(0.5 * (c[1:] + c[:-1]))))

with open("INPUT/properties.dat", "w") as f:
    f.write('TITLE = "Mass Thermodynamic Properties"\n')
    f.write('VARIABLES = "Temperature", "Cp", "Density", "Conductivity", "Energy"\n')
    f.write('ZONE T="mymat"\n')
    f.write(f'I={len(T_arr)}, F=POINT\n')
    for i, T in enumerate(T_arr):
        f.write(f"{T:.1f} {cp_arr[i]:.6f} {rho_arr[i]:.6f} "
                f"{k_arr[i]:.6f} {h_arr[i]:.6f}\n")
```

The companion `INPUT/phase.txt` needs a first line (skipped) and one line per material:

```
solid-bulk
mymat 1
```

!!! note "Constant-property shortcut"

    When the `[GPB-Phase1]` section of `input.ini` gives scalar values for `cp`, `rho`, and `k`, ATLAS GPB writes a table with the same values in every row from `Tmin` to `Tmax` and the Energy column $\rho\,c_p\,T$.  A table with temperature-dependent columns can be written directly, as above.

---

## Phase Definition (`input.ini`)

The `[GPB-Phase1]` section, read by ATLAS, declares the materials used in the simulation:

```ini
[GPB-Phase1]
type     = solid
material = steel aluminium
cp       = 460.0  900.0
rho      = 7850.0 2700.0
k        = 52.0   200.0
Tmax     = 5000
```

The `material` list and the scalar property arrays (`cp`, `rho`, `k`) must be ordered consistently; ATLAS writes the ZONEs of `properties.dat` and the lines of `phase.txt` in the order of `material`.  `Tmin` and `Tmax` set the temperature range of the table, and with it the admissible range of the solution (see [Temperature range](#property-interpolation)).

---

## Multi-Material Block Assignment

Each block is assigned a material through its initial-condition section, which ATLAS turns into the per-cell material index of the initial-condition file `INPUT/ic.*`:

```ini
[ICB-Block1]
type     = homogeneous
material = steel
T        = 300.0

[ICB-Block2]
type     = homogeneous
material = aluminium
T        = 300.0
```

The `material` keyword must match one of the names declared in `[GPB-Phase1]`.  A material may also change inside a block: the face between two materials takes the series conductivity of the two half cells (see [Face Conductivity](numerics.md#face-conductivity)).

---

## Multi-Material Interface Treatment

At a block connection between two materials each side keeps its own conductivity, $k_P(T_P)$ and $k_N(T_N)$, over its half cell, and the interface temperature follows from the continuity of the heat flux.  In one dimension the flux is that of two resistances in series, $(T_N - T_P)/\bigl(\Delta x_P/(2k_P) + \Delta x_N/(2k_N)\bigr)$, the harmonic mean of $k_P$ and $k_N$ for equal cell sizes.  Since $k$ depends on $T$, the interface conductances change as the solution evolves.  The discretisation is given in [Spatial Discretisation](numerics.md#block-connections-101-and-chimera-faces-102); the [multi-material plate](../vv/2D-transient-multimat-plate.md) tests it.

---

## Volumetric Heat Generation

Each cell carries a heat source $\dot{q}_\text{vol}$ (W/m³), read from `INPUT/st.*` when that file exists (zero otherwise).  ATLAS (STB) writes it from the `[STB-Block<n>]` sections, one value per block:

```ini
[STB-Block1]
qvol = 5.0e5
```

The source enters every stage as $\Delta t\,\dot q_\text{vol}$ on $h$ (divided by $\rho\,c_p$ on $T$ with `prim`).  FUSS itself makes it depend neither on time nor on temperature; on coarse multigrid levels it is the volume average of the fine cells.

---

## References

1. F. P. Incropera, D. P. DeWitt, T. L. Bergman, A. S. Lavine, *Fundamentals of Heat and Mass Transfer*, 7th ed., Wiley, 2011.
2. H. S. Carslaw, J. C. Jaeger, *Conduction of Heat in Solids*, 2nd ed., Oxford University Press, 1959.
3. C. J. Greenshields, H. G. Weller, *Notes on Computational Fluid Dynamics: General Principles*, CFD Direct Ltd., 2022.
