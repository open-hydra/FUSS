# Governing Equations

FUSS solves the **heat equation for solid continua** on structured multi-block grids.  The single unknown is the temperature field $T(\mathbf{x}, t)$.  Both transient and steady-state problems are supported.

---

## Transient Heat Equation

The governing equation for heat conduction in a solid with volumetric heat generation is

$$
\rho(T)\, c_p(T)\, \frac{\partial T}{\partial t}
= \nabla\!\cdot\!\bigl(k(T)\,\nabla T\bigr) + \dot{q}_\text{vol}
$$

where

| Symbol | Quantity | Unit |
|--------|----------|------|
| $\rho(T)$ | Density | kg/m³ |
| $c_p(T)$ | Specific heat capacity | J/(kg·K) |
| $T$ | Temperature | K |
| $k(T)$ | Thermal conductivity | W/(m·K) |
| $\dot{q}_\text{vol}$ | Volumetric heat generation rate | W/m³ |

All three material properties — $k$, $\rho$, and $c_p$ — can depend on temperature.  When they do, the equation is **non-linear**: the diffusivity $\alpha(T) = k(T)/[\rho(T)\,c_p(T)]$ changes with the local solution, and $\nabla\!\cdot\!(k(T)\,\nabla T)$ contains the term $\mathrm{d}k/\mathrm{d}T\,\lvert\nabla T\rvert^2$.  The explicit time step of each cell uses the $\alpha$ of its current temperature (see [Time Integration](time-integration.md#diffusive-stability-condition-vnn)).

Properties are tabulated per material against temperature and interpolated linearly (see [Material Properties](thermo.md)).  $\dot q_\text{vol}$ is given per cell in `INPUT/st.*` (zero without that file); FUSS itself keeps it constant.

---

## Steady-State Limit

When `time-accurate = false`, FUSS marches the transient equation to convergence with a local time step per cell, optionally accelerated by implicit residual smoothing and grid sequencing (see [Time Integration](time-integration.md)).  The converged solution satisfies the steady-state heat equation

$$
\nabla\!\cdot\!\bigl(k\,\nabla T\bigr) + \dot{q}_\text{vol} = 0
$$

The convergence measure is the Euclidean norm, over all cells of all blocks, of the temperature change in one iteration,

$$
r = \Bigl[\sum_{i} \bigl(T_i^{n+1} - T_i^{n}\bigr)^2\Bigr]^{1/2},
$$

evaluated every `res-diter` iterations.  A steady run stops when $r$ has fallen to `res-threshold` times its value at the first iteration (of the grid level), or when the iteration limit is reached; a time-accurate run stops at the time or the iteration limit.

---

## Finite Volume Integral Form

In conservative form the stored quantity is the volumetric enthalpy

$$
h(T) = \int^{T}\!\rho(T')\,c_p(T')\,\mathrm{d}T', \qquad \frac{\mathrm{d}h}{\mathrm{d}T} = \rho(T)\,c_p(T),
$$

the **Energy** column of `properties.dat` (see [Material Properties](thermo.md)), so that $\partial h/\partial t = \rho\,c_p\,\partial T/\partial t$ also with temperature-dependent properties.  Integrating over a cell $\Omega_i$ and applying the divergence theorem gives

$$
\frac{\mathrm{d}}{\mathrm{d}t}\!\int_{\Omega_i}\!h(T)\,\mathrm{d}V
= \oint_{\partial\Omega_i}\! k(T)\,\nabla T\!\cdot\!\hat{\mathbf{n}}\,\mathrm{d}A
\;+\;\dot{q}_\text{vol}\,V_i
$$

Note that $h$ is not $\rho\,c_p\,T$ unless the properties are constant: with $\rho\,c_p$ varying, $\mathrm{d}(\rho\,c_p\,T)/\mathrm{d}T = \rho\,c_p + T\,\mathrm{d}(\rho\,c_p)/\mathrm{d}T$ would add a spurious term to the heat capacity.

Approximating all integrals with midpoint quadrature, and evaluating material properties at the cell-centre temperature $T_i$, yields the **semi-discrete update**

$$
\rho(T_i)\,c_p(T_i)\,\frac{\mathrm{d}T_i}{\mathrm{d}t}
= \frac{1}{V_i}\sum_{f=1}^{N_f}\!\bigl(k_f\,\nabla T\!\cdot\!\hat{\mathbf{n}}\bigr)_f\,A_f
\;+\;\dot{q}_{\text{vol},i}
$$

where $V_i$ is the cell volume, $f$ runs over the (up to six) faces of the hexahedral cell, $A_f$ is the face area, $\hat{\mathbf{n}}$ is the outward unit normal, and $k_f$ is the face conductivity.  How $\nabla T$ and $k_f$ are evaluated at interior, block-connection and boundary faces is described in [Spatial Discretisation](numerics.md).

---

## State Vector

FUSS advances a single cell-centred variable, the temperature:

$$
\mathbf{P} = [T]
$$

Each cell also carries a material index, read with $T$ from the initial-condition file, and $\dot q_\text{vol}$ (see [Material Properties](thermo.md#volumetric-heat-generation)).  There is no velocity, pressure, or species field.  The key `integration-variables` selects what the Runge–Kutta stages advance; both options integrate the same equation:

| `integration-variables` | Advanced quantity | Update of a stage |
|---|---|---|
| `cons` (default) | volumetric enthalpy $h_i = h(T_i)$ | $h_i$ by $\Delta t\,[\tfrac{1}{V_i}\sum_f (k_f\nabla T\cdot\hat{\mathbf n})_f A_f + \dot q_{\text{vol},i}]$, then $T_i$ by inverting the $h(T)$ table |
| `prim` | temperature $T_i$ | the same increment divided by $\mathrm{d}h/\mathrm{d}T = \rho(T_i)\,c_p(T_i)$ |

The two agree as long as the Energy column is the integral of $\rho\,c_p$, as ATLAS GPB writes it; FUSS checks the table when it loads it and warns when $\mathrm{d}h/\mathrm{d}T$ differs from $\rho\,c_p$ by more than 1 % (see [Material Properties](thermo.md#energy-column-and-the-enthalpy-inversion)).  Before 2026-10-09 the `prim` update divided by $\mathrm{d}(\rho\,c_p\,T)/\mathrm{d}T$, which added $T\,\mathrm{d}(\rho\,c_p)/\mathrm{d}T$ to the heat capacity of temperature-dependent materials (see the [variable-$c_p$ plate](../vv/2D-transient-variable-cp-plate.md)).

---

## Boundary Conditions

The solver reads one boundary condition per boundary cell face from `INPUT/bc.txt` (`INPUT/bc<m>.txt` for the coarse multigrid level $m$), so the data may vary from cell to cell along a block face.  Each entry carries a numeric type code.

**Sign convention.** With $\hat{\mathbf n}$ the outward normal of the solid,

$$
q_\text{in} = k\,\nabla T\cdot\hat{\mathbf n}
$$

is the heat flux **into the solid** (W/m²): positive when the solid is heated through the face.  $T_w$ denotes the temperature at the face centre.

| Code | Condition | Face relation | `bc.txt` data line |
|:----:|-----------|---------------|--------------------|
| `0` | null | $q_\text{in} = 0$ | — |
| `300` | symmetry | $q_\text{in} = 0$ | — |
| `101` | block connection | $T$ and $k\,\nabla T\cdot\hat{\mathbf n}$ continuous | `bs is js ks fs d11 d12 d21 d22` |
| `102` | chimera | ghost cells from overlapping donor cells | `n1 n2`, then `n1+n2` lines `bD iD jD kD volume_fraction` |
| `103`, `104` | coupling with the hydra fluid solver | $q_\text{in}$ set by hydra | as `101` / as `102` |
| `301` | heat flux | $q_\text{in} = q_w$ | $q_w$, or a table file |
| `302` | temperature | $T_w$ prescribed | $T_w$, or a table file |
| `303` | convection (+ $q_\text{rad}$) | $q_\text{in} = h_\text{conv}\,(T_\text{ref} - T_w) + q_\text{rad}$ | `hconv qrad Tref` |
| `304` | radiation | $q_\text{in} = \varepsilon\,\sigma\,(T_\text{ref}^4 - T_w^4)$ | `eps Tref` |
| `305` | convection + radiation | $q_\text{in} = h_\text{conv}\,(T_\text{ref} - T_w) + \varepsilon\,\sigma\,(T_\text{ref}^4 - T_w^4)$ | `hconv eps Tref` |

with $\sigma = 5.67\times10^{-8}$ W/(m²·K⁴).  Null and symmetry faces are adiabatic; both mirror the interior cells into the ghost cells.  A heat-flux wall with $q_w = 0$ is adiabatic as well.  For 301 and 303–305 the wall temperature $T_w$ is an unknown: it is the face temperature at which the conduction flux of the solid matches the prescribed $q_\text{in}$, and it is written to the wall output with $q_\text{in}$.  The discrete form of every condition, including the wall ghost cells, is given in [Spatial Discretisation](numerics.md#boundary-faces).  Codes 103 and 104 are only meaningful when FUSS runs inside hydra.

**Time tables.**  For 301 and 302 the data line may hold a file name instead of a number.  The file has two columns, time and value; if its first line contains the word `periodic`, the table repeats with period $t_\text{last} - t_\text{first}$.  Values are interpolated linearly, held at the first value before the first time and at the last value after the last time.  The table is evaluated at every ghost-cell fill; in a time-accurate run stage $s$ of the Runge–Kutta scheme reads it at $t^n + c_s\,\Delta t$ (the stage time), which preserves the temporal order of the scheme.  A run with a time table does not stop on `res-threshold` while `time-accurate = true`: the table keeps forcing the solution, so only `time-threshold` or `iter-threshold` end it.

**ATLAS keys.**  ATLAS writes `bc.txt` from the `type = wall` sections of `input.ini` (see [Boundary Conditions](../user/boundary-conditions.md)): `q` or `q-time-file` → 301, `T` or `T-time-file` → 302, `hconv` + `qrad` + `Tref` → 303 (pure convection needs `qrad = 0` written), `eps` + `Tref` → 304, `hconv` + `eps` + `Tref` → 305.  ATLAS 7d3c72d tests the 304 combination first, so a section with `hconv`, `eps` and `Tref` is written as 304 and `hconv` is dropped: code 305 has to be written into `bc.txt` by hand.  The face keywords `null` and `symmetry` give 0 and 300; `connection` and `chimera` give 101 or 102.

---

## Multi-Block Structure

FUSS supports an arbitrary number of structured hexahedral blocks.  Each cell carries a material index, an initial temperature and $\dot q_\text{vol}$; each boundary cell face carries a boundary condition.  Blocks exchange heat through block connections (101) and chimera overlaps (102), and the solver treats the assembled multi-block domain as a single computational grid.  A material change must coincide with a block boundary (see [Spatial Discretisation](numerics.md#face-conductivity)).

---

## References

1. O. C. Zienkiewicz, R. L. Taylor, *The Finite Element Method*, Vol. 1, 6th ed., Butterworth-Heinemann, 2005.
2. E. Kreyszig, *Advanced Engineering Mathematics*, 10th ed., Wiley, 2011.
3. F. P. Incropera, D. P. DeWitt, T. L. Bergman, A. S. Lavine, *Fundamentals of Heat and Mass Transfer*, 7th ed., Wiley, 2011.
