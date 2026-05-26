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

All three material properties — $k$, $\rho$, and $c_p$ — can depend on temperature.  When they do, the equation is **non-linear**: the effective diffusivity $\alpha(T) = k(T)/[\rho(T)\,c_p(T)]$ changes with the local solution, and the diffusion term $\nabla\!\cdot\!(k(T)\,\nabla T)$ introduces a non-linearity through both the varying coefficient and the $\nabla k \cdot \nabla T$ cross-term that emerges when $k$ is not spatially uniform.  The explicit time-step limit is evaluated using the cell-local $\alpha$ and the global minimum is used (see [Time Integration](time-integration.md)).

Properties are looked up at every residual evaluation by linear interpolation of a precomputed table indexed by temperature (see [Material Properties](thermo.md)).

---

## Steady-State Limit

When `time-accurate = false`, FUSS marches the transient equation to convergence using local time stepping or implicit residual smoothing.  The converged solution satisfies the steady-state heat equation

$$
\nabla\!\cdot\!\bigl(k\,\nabla T\bigr) + \dot{q}_\text{vol} = 0
$$

Convergence is monitored through the residual of this equation.  The run terminates when the maximum residual falls below the user-specified threshold `res-threshold`.

---

## Finite Volume Integral Form

Integrating the transient heat equation over a cell $\Omega_i$ and applying the divergence theorem gives

$$
\frac{\mathrm{d}}{\mathrm{d}t}\!\int_{\Omega_i}\!\rho(T)\,c_p(T)\,T\,\mathrm{d}V
= \oint_{\partial\Omega_i}\! k(T)\,\nabla T\!\cdot\!\hat{\mathbf{n}}\,\mathrm{d}A
\;+\;\dot{q}_\text{vol}\,V_i
$$

Approximating all integrals with midpoint quadrature, and evaluating material properties at the cell-centre temperature $T_i$, yields the **semi-discrete update**

$$
\rho(T_i)\,c_p(T_i)\,\frac{\mathrm{d}T_i}{\mathrm{d}t}
= \frac{1}{V_i}\sum_{f=1}^{N_f}\!\bigl(k_f\,\nabla T\!\cdot\!\hat{\mathbf{n}}\bigr)_f\,A_f
\;+\;\dot{q}_{\text{vol},i}
$$

where $V_i$ is the cell volume, $f$ runs over the (up to six) faces of the hexahedral cell, $A_f$ is the face area, $\hat{\mathbf{n}}$ is the outward unit normal, and $k_f$ is the face conductivity evaluated from the temperatures of the two adjacent cells (see [Spatial Discretisation](numerics.md) and [Material Properties](thermo.md)).

---

## State Vector

FUSS stores a single cell-centred variable per cell:

$$
\mathbf{P} = [T]
$$

Unlike CFD solvers, there is no velocity, pressure, or species field.  The integration variable is always temperature (selected via `integration-variables = prim` in the input file).

---

## Boundary Conditions

Boundary conditions are enforced through ghost-cell values that are set before each residual evaluation. At the INI level all wall BCs use `type = wall` and the **physical behaviour is selected by the parameter combination** in the named section. The pre-processor translates each parameter combination into an internal numeric code consumed by the solver. The catalogue of physical behaviours is:

| Behaviour | Mathematical form | Parameter combination | Internal code |
|-----------|------------------|-----------------------|:-------------:|
| Block connection | continuity of $T$ and $k\,\nabla T\!\cdot\!\hat{\mathbf n}$ | `face* = connection` | `101` |
| Symmetry | $k\,\nabla T\!\cdot\!\hat{\mathbf n} = 0$ | (set by ATLAS at symmetry planes) | `300` |
| Prescribed heat flux | $-k\,\nabla T\!\cdot\!\hat{\mathbf n} = q_w$ | `q` (or `q-time`) | `301` |
| Prescribed temperature | $T_\text{wall} = T_w$ (constant or $T_w(t)$) | `T` (or `T-time`) | `302` |
| Convection (+ optional $q_\text{rad}$) | $q = h_\text{conv}(T_\text{ref} - T_\text{wall}) + q_\text{rad}$ | `hconv`, `Tref`, `qrad` | `303` |
| Radiation | $q = \varepsilon\,\sigma(T_\text{ref}^4 - T_\text{wall}^4)$ | `eps`, `Tref` | `304` |
| Convection + radiation | $q = h_\text{conv}(T_\text{ref} - T_\text{wall}) + \varepsilon\,\sigma(T_\text{ref}^4 - T_\text{wall}^4)$ | `hconv`, `eps`, `Tref` | `305` |

The `null` keyword (used on degenerate 2-D z-faces) deactivates the face entirely. Adiabatic walls are expressed as a prescribed-flux BC with `q = 0`.

For a step-by-step description of the INI assignment and full BC examples see [Boundary Conditions](../user/boundary-conditions.md).

### Prescribed temperature

A constant wall temperature is specified with

```ini
[MyBC]
type = wall
T    = 1273.15
```

A time-varying wall temperature is specified by pointing to a two-column ASCII table:

```ini
[MyBC]
type   = wall
T-time = Tw_time.dat
```

The file can optionally begin with the keyword `periodic`, in which case the table is repeated cyclically for the entire simulation. Linear interpolation is used between table entries.

### Prescribed flux

A constant heat flux is specified with

```ini
[MyBC]
type = wall
q    = 1.0e5     ; W/m² (positive = into the solid)
```

Setting `q = 0` yields an adiabatic wall. A time-varying flux is specified via `q-time = Q_time.dat`.

### Convection

$$
q_\text{conv} = h_\text{conv}\,(T_\text{ref} - T_\text{wall}) + q_\text{rad}
$$

Parameters `hconv` (W/(m²·K)), `Tref` (K), and the optional radiative pre-load `qrad` (W/m²) are set in the input file. The sign convention is positive inward (heat into the solid when $T_\text{ref} > T_\text{wall}$).

### Radiation

$$
q_\text{rad} = \varepsilon\,\sigma\,(T_\text{ref}^4 - T_\text{wall}^4)
$$

with the Stefan–Boltzmann constant $\sigma = 5.67 \times 10^{-8}$ W/(m²·K⁴), surface emissivity $\varepsilon \in (0, 1]$, and radiation sink temperature $T_\text{ref}$.

!!! note "Non-linearity of radiation"

    Pure-radiation and combined convection+radiation BCs introduce a strong non-linearity in $T$ (quartic dependence). When these BCs are active, smaller VNN values (more conservative time stepping) are recommended to maintain stability.

### Combined convection + radiation

When both heat-transfer modes are present at the same face, the section carries `hconv`, `eps`, and `Tref` together; the solver applies the sum of both fluxes.

### Block interface

At internal block boundaries (`face* = connection`), FUSS enforces continuity of both temperature and normal heat flux. The face conductivity is computed as the harmonic mean of the two adjacent cell conductivities (see [Spatial Discretisation](numerics.md)), which correctly represents a series resistance across the interface even when the two materials differ strongly in conductivity.

---

## Multi-Block Structure

FUSS supports an arbitrary number of structured hexahedral blocks.  Each block carries an independent material assignment (phase), initial condition, and set of boundary conditions on its six faces.  Blocks communicate through `connection` interfaces (internal code 101), and the solver treats the assembled multi-block domain as a single computational grid.

---

## References

1. O. C. Zienkiewicz, R. L. Taylor, *The Finite Element Method*, Vol. 1, 6th ed., Butterworth-Heinemann, 2005.
2. E. Kreyszig, *Advanced Engineering Mathematics*, 10th ed., Wiley, 2011.
3. F. P. Incropera, D. P. DeWitt, T. L. Bergman, A. S. Lavine, *Fundamentals of Heat and Mass Transfer*, 7th ed., Wiley, 2011.
