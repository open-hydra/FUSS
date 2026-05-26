# Material Properties

FUSS models heat conduction in solid media.  Each computational block is assigned a **material phase** that carries the three thermophysical properties required by the heat equation: density $\rho$, specific heat capacity $c_p$, and thermal conductivity $k$.  All three properties can vary with temperature; they are evaluated at every cell and every residual call by linear interpolation of a precomputed table.

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

When properties vary with temperature, $\alpha(T)$ changes cell-by-cell and the governing equation becomes non-linear.  The explicit time-step limit is set conservatively using the maximum $\alpha$ across the domain (see [Time Integration](time-integration.md)).

---

## Property Table: `properties.dat`

Material properties are not read directly from the scalar values in `input.ini`.  Instead, the solver always reads them from a precomputed **Tecplot ASCII table** located at

```
INPUT/properties.dat
```

This file contains one Tecplot ZONE per material, each zone providing a look-up table indexed by **integer temperature in Kelvin**.

### File format

```
TITLE = "Mass Thermodynamic Properties"
VARIABLES = "Temperature", "Cp", "Density", "Conductivity", "Energy"
ZONE T="<material_name>"
I=<N>, F=POINT
<T_start>.0  <cp_1>  <rho_1>  <k_1>  <h_1>
<T_start+1>.0  <cp_2>  <rho_2>  <k_2>  <h_2>
...
```

| Column | Symbol | Unit | Notes |
|--------|--------|------|-------|
| Temperature | $T$ | K | Integer values: $T_\text{start}, T_\text{start}+1, \ldots, T_\text{start}+N-1$ |
| Cp | $c_p(T)$ | J/(kg·K) | Specific heat at constant pressure |
| Density | $\rho(T)$ | kg/m³ | Material density |
| Conductivity | $k(T)$ | W/(m·K) | Thermal conductivity |
| Energy | $h(T)$ | J/m³ | Cumulative volumetric enthalpy (see below) |

The table must start at an integer temperature and contain **one row per integer Kelvin**.  A common choice is $T_\text{start} = 1$ K running up to `Tmax` K specified in the `[GPB-Phase1]` section of the input file.

### Multiple materials

When the simulation uses more than one material, each material occupies its own ZONE, ordered in the same sequence as the names listed in `[GPB-Phase1]`:

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

---

## Property Interpolation

At every residual evaluation, FUSS calls dedicated routines (`co_Kappa`, `co_Rho`, `co_Cp`) that perform **linear interpolation** between the two table rows that bracket the current cell temperature:

$$
k(T) = k(\lfloor T \rfloor)
       + \bigl(T - \lfloor T \rfloor\bigr)\,
         \bigl[k(\lfloor T \rfloor + 1) - k(\lfloor T \rfloor)\bigr]
$$

and analogously for $\rho$ and $c_p$.  Because the table is indexed by integer Kelvin, the floor and ceiling indices are simply $T_m = \lfloor T \rfloor$ and $T_p = T_m + 1$, so no binary search is needed — interpolation is $O(1)$.

!!! warning "Temperature range"

    If the cell temperature falls outside the range covered by the table, the indexing will read out-of-bounds memory.  The `Tmax` parameter in `[GPB-Phase1]` should be set comfortably above the maximum expected temperature to prevent this.  Typical practice is to add 20–30 % margin above the highest anticipated temperature.

---

## Energy Column and the Enthalpy Inversion

The **Energy** column stores the cumulative volumetric enthalpy

$$
h(T_n) = \sum_{i=T_\text{start}}^{T_n} \rho(i)\,c_p(i)\;\Delta T, \qquad \Delta T = 1\;\text{K}
$$

which is the discrete approximation of $\int_{T_\text{start}}^{T} \rho\,c_p\,\mathrm{d}T'$.

This column is used by the `co_H` and `co_T` routines to convert between temperature and stored energy — a path that arises when the integration variable is switched to conservative (energy) form.  For constant properties the formula reduces to $h(T) = \rho\,c_p\,T$, as seen in the single-material test cases:

$$
h = 7850\;\text{kg/m}^3 \times 460\;\text{J/(kg·K)} \times 1\;\text{K} = 3\,611\,000\;\text{J/m}^3
$$

For **variable properties** each row must accumulate the product $\rho(T)\,c_p(T)$ up to that temperature, integrating the actual variation.

---

## Variable Properties: How to Build the Table

To exploit the temperature-dependent capability, provide a `properties.dat` file in which the Cp, Density, and Conductivity columns vary from row to row.  A Python snippet to generate the file for a material with known $k(T)$, $\rho(T)$, $c_p(T)$ functions:

```python
import numpy as np

T_start = 1       # first integer temperature [K]
T_end   = 2000    # last  integer temperature [K]
T_arr   = np.arange(T_start, T_end + 1, dtype=float)

# Define property functions (replace with real data or polynomial fits)
def k_func(T):   return 52.0 - 0.01 * (T - 300)          # W/(m·K)
def rho_func(T): return 7850.0 - 0.05 * (T - 300)        # kg/m³
def cp_func(T):  return 460.0  + 0.20 * (T - 300)        # J/(kg·K)

k_arr   = k_func(T_arr)
rho_arr = rho_func(T_arr)
cp_arr  = cp_func(T_arr)

# Cumulative volumetric enthalpy (trapezoid rule, ΔT = 1 K)
h_arr = np.cumsum(rho_arr * cp_arr * 1.0)

with open("INPUT/properties.dat", "w") as f:
    f.write('TITLE = "Mass Thermodynamic Properties"\n')
    f.write('VARIABLES = "Temperature", "Cp", "Density", "Conductivity", "Energy"\n')
    f.write(f'ZONE T="mymat"\n')
    f.write(f'I={len(T_arr)}, F=POINT\n')
    for i, T in enumerate(T_arr):
        f.write(f"{T:.1f} {cp_arr[i]:.6f} {rho_arr[i]:.6f} "
                f"{k_arr[i]:.6f} {h_arr[i]:.6f}\n")
```

The companion `INPUT/phase.txt` file must list the material name with its sequential index:

```
solid phase
mymat 1
```

!!! note "Constant-property shortcut"

    When the `[GPB-Phase1]` section in `input.ini` provides scalar values for `cp`, `rho`, and `k`, the preprocessor expands them into a constant-valued table that spans `[1, Tmax]` K.  The resulting `properties.dat` has the same format but identical values in every row — effectively making $k$, $\rho$, and $c_p$ constant.  Users who need temperature-varying properties must supply their own `properties.dat` rather than relying on the scalar shortcut.

---

## Phase Definition (`input.ini`)

The `[GPB-Phase1]` section declares the materials used in the simulation:

```ini
[GPB-Phase1]
type     = solid
material = steel aluminium
cp       = 460.0  900.0
rho      = 7850.0 2700.0
k        = 52.0   200.0
Tmax     = 5000
```

The `material` list and the scalar property arrays (`cp`, `rho`, `k`) must be ordered consistently and match the ZONE order in `properties.dat`.  `Tmax` sets the upper temperature bound of the table; the solver clamps cell temperatures to `[1, Tmax]` to prevent out-of-range table access.

---

## Multi-Material Block Assignment

Each block is independently assigned a material through its initial condition section:

```ini
[ICB-Block1]
type     = homogeneous
material = steel
T        = 300.0
qvol     = 0.0

[ICB-Block2]
type     = homogeneous
material = aluminium
T        = 300.0
qvol     = 0.0
```

The `material` keyword must match one of the names declared in `[GPB-Phase1]` (and present as a ZONE in `properties.dat`).  The solver stores an integer **material ID** at every cell, which is used as the first index into the property table when evaluating $k$, $\rho$, and $c_p$ during the residual loop.

---

## Multi-Material Interface Treatment

When two adjacent blocks carry different materials, the diffusive flux at their shared face uses the **harmonic mean conductivity**

$$
k_f = \frac{2\,k_L(T_L)\,k_R(T_R)}{k_L(T_L) + k_R(T_R)}
$$

where $k_L$ and $k_R$ are evaluated at the temperatures of the two cells on either side of the face.  Because $k$ can now vary with $T$, the face conductivity changes each time step as the solution evolves.  The harmonic mean remains the correct choice for a series resistance regardless of how the individual conductivities depend on temperature.

---

## Volumetric Heat Generation

Each block can include a uniform internal heat source $\dot{q}_\text{vol}$ (W/m³) specified in the initial condition section:

```ini
[ICB-Block1]
type     = homogeneous
material = mat
T        = 273.15
qvol     = 5.0e5
```

The source term is added directly to the cell residual and has no temperature dependence in the current implementation.

---

## References

1. F. P. Incropera, D. P. DeWitt, T. L. Bergman, A. S. Lavine, *Fundamentals of Heat and Mass Transfer*, 7th ed., Wiley, 2011.
2. H. S. Carslaw, J. C. Jaeger, *Conduction of Heat in Solids*, 2nd ed., Oxford University Press, 1959.
3. C. J. Greenshields, H. G. Weller, *Notes on Computational Fluid Dynamics: General Principles*, CFD Direct Ltd., 2022.
