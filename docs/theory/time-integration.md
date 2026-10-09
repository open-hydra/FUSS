# Time Integration

FUSS advances the semi-discrete heat equation in time with an explicit **strong-stability-preserving Runge–Kutta** scheme.  This page describes the time-stepping algorithm, the diffusive stability condition that limits the time step, the two operating modes (steady-state and time-accurate), and the optional implicit residual smoothing that can accelerate convergence for steady problems.

---

## SSP Runge–Kutta Scheme

The default (and only) time integrator is the **3-stage SSP RK3** scheme (Shu–Osher):

$$
\begin{aligned}
T^{(1)} &= T^n + \Delta t\;\mathcal{L}(T^n)\\[4pt]
T^{(2)} &= \tfrac{3}{4}\,T^n
  + \tfrac{1}{4}\bigl(T^{(1)} + \Delta t\;\mathcal{L}(T^{(1)})\bigr)\\[4pt]
T^{n+1} &= \tfrac{1}{3}\,T^n
  + \tfrac{2}{3}\bigl(T^{(2)} + \Delta t\;\mathcal{L}(T^{(2)})\bigr)
\end{aligned}
$$

where $\mathcal{L}(T)$ is the spatial residual of the heat equation:

$$
\mathcal{L}(T_i) =
\frac{1}{\rho_i\,c_{p,i}}\!\left[
\frac{1}{V_i}\sum_{f} Q_f\,A_f
\;+\;\dot{q}_{\text{vol},i}
\right]
$$

This is the `prim` form.  With `integration-variables = cons` (the default) the stages advance the volumetric enthalpy $h_i$ with the bracket alone, and $T_i$ is recovered from the $h(T)$ table after each stage; since $\mathrm{d}h/\mathrm{d}T = \rho\,c_p$ the two forms integrate the same equation (see [Governing Equations](governing-equations.md#state-vector)).

The scheme is:

- **Third-order accurate** in time for smooth problems
- **SSP (strong-stability-preserving)**: each stage is a convex combination of forward-Euler steps, preserving positivity and monotonicity inherited from a TVD spatial operator
- **TVD** under the same stability limit as forward Euler

---

## Diffusive Stability Condition (VNN)

Pure diffusion is governed by the **von Neumann number** (VNN) condition.  For cell $i$ in coordinate direction $d$, the stable time-step limit is

$$
\Delta t_{\text{VNN},i}^{(d)} =
\text{VNN} \times \frac{\rho(T_i)\,c_p(T_i)\,(\Delta x_{i}^{(d)})^2}{k(T_i)}
= \text{VNN} \times \frac{(\Delta x_{i}^{(d)})^2}{\alpha(T_i)}
$$

where $\alpha(T_i) = k(T_i) / [\rho(T_i)\,c_p(T_i)]$ is the **local** thermal diffusivity evaluated at the current cell temperature, and $\Delta x_{i}^{(d)}$ is the cell width in direction $d$.  When material properties vary with temperature, each cell computes its own diffusivity from the property table at every time step, so the global limit tightens automatically in regions where $\alpha$ is large.

The global time step is the minimum over all cells and all three coordinate directions:

$$
\Delta t = \min_{i,\,d}\;\Delta t_{\text{VNN},i}^{(d)}
$$

!!! tip "Choosing VNN"

    For the 3-stage SSP RK3 scheme the stability limit for the 1-D diffusion equation is $\text{VNN} \leq 1$.  In practice, values between 0.5 and 0.8 are used to maintain a safety margin on non-uniform grids.  For steady-state runs with IRS enabled, VNN can be raised to 2 or higher without loss of stability (the smoothing effectively extends the stable region).

---

## Operating Modes

### Time-accurate mode (`time-accurate = true`)

The global minimum time step is used at every iteration.  The solver advances to the physical time `time-threshold` specified in the input file, writing field snapshots at intervals controlled by `sol-dtime` or `sol-diter`.

This mode is required for all transient problems: oscillating boundary conditions, evolving temperature profiles, probe monitoring.

### Steady-state mode (`time-accurate = false`)

Each cell uses its **local** time step (the cell's own VNN limit without taking the global minimum).  This is equivalent to pseudo-time marching and accelerates convergence to steady state by large factors on grids with widely varying cell sizes.

Convergence is declared when the maximum residual (absolute or relative, user-selectable) falls below `res-threshold`.

!!! warning "Local time stepping and transient accuracy"

    Local time stepping destroys temporal accuracy.  It must not be used when the transient history is physically meaningful.

---

## Implicit Residual Smoothing (IRS)

IRS increases the effective stability limit of the explicit RK scheme by applying a Laplacian smoothing operator to the residual before the state update.  It is particularly effective for steady-state problems on stretched meshes.

### Smoothing equation

In each coordinate direction $d$ the smoothed residual $\mathcal{L}^\ast$ satisfies the implicit tridiagonal system:

$$
\mathcal{L}^\ast_i =
\frac{\mathcal{L}_i
  + \varepsilon\,\bigl(\mathcal{L}^\ast_{i-1} + \mathcal{L}^\ast_{i+1}\bigr)}
{1 + 2\varepsilon}
$$

with smoothing parameter $\varepsilon$ specified by `irs-beta` in the input file.  The system is solved approximately with **2 Jacobi iterations**.

### Effect on stability

The smoothing amplifies the stable VNN limit by approximately $(1 + 2\varepsilon)$ per direction.  With $\varepsilon = 0.5$ (a common choice) the effective limit is $(1 + 1)^3 \approx 8$ times the un-smoothed value in 3-D, allowing VNN values of 2 or more to be used without instability.

### Admissible updates

The smoothed update of a cell mixes the updates of its neighbours. Above the explicit limit it can overshoot the temperatures imposed at the boundaries and leave the range of the properties table, where the property lookups would read past its ends. Each smoothed update is therefore checked: the new temperature must be finite and lie in $[T_\text{min}, T_\text{max})$, the first and last temperature of `properties.dat`. A cell that fails takes its own, unsmoothed update instead, halved until the new state is admissible (at most 10 times), and otherwise keeps no update in that stage: every stage is then a convex combination of admissible states. The number of such cells, summed over the MPI ranks, is printed every `shell-diter` iterations. An admissible smoothed update is applied unchanged.

Without IRS nothing is replaced: a temperature outside the table, or NaN, stops the run with an error message naming the cell, and a non-zero exit status.

### Input settings

```ini
[FUSS-Numerics]
irs      = true
irs-beta = 0.5
```

IRS is switched on at setup, by `irs = true` or by `irs-beta > 0`; only then are its arrays allocated. During the run `irs-beta` and `irs = false` take effect at the next re-read of `input.ini` (every `ini-diter` iterations), while `irs = true` added to a run started without IRS is ignored with a warning. The default is no smoothing.

---

## Multigrid (grid sequencing)

With `levels = n` (n > 1) in `[FUSS-Multigrid]`, FUSS builds n grids, each coarser level merging 2 × 2 × 2 cells (2 × 2 in 2-D) of the finer one, so the cell counts of every block must be divisible by $2^{n-1}$. The run starts on the coarsest level, from the initial field restricted to it (volume average of $h$). Level $m$ runs `level<m>-iter` iterations, or until the residual reaches `res-threshold`; its solution is then prolongated (trilinear interpolation) to the next finer level, where the run continues. The finest level ends after `level1-iter` iterations or at `res-threshold` / `time-threshold`: with multigrid, `iter-threshold` is not used. A coarse level is not visited again (no multigrid cycles).

`level<m>-iter`, and `iter-threshold` without multigrid, are re-read with `input.ini` every `ini-diter` iterations: lowering the value of the level being computed below the iterations already done ends that level at the next re-read.

Every field write of a coarse level also prolongates its solution to the finest level and writes it (`field-prolongated`, `wall-prolongated`). On a large mesh that write can cost more than the coarse iterations between two writes: `prolongated-diter` (re-read at runtime) spaces it to at least that many iterations of the level. The end of a level always writes it, so the one written at the end of level 2 is the exact start of level 1.

---

## Solution Output

### Steady-state runs

A single field file is written (with `sol-overwrite = true`) every `sol-diter` iterations.  This overwrites the previous snapshot, retaining only the latest state.

### Transient runs

With `sol-overwrite = false`, a new numbered file is written every `sol-diter` iterations or every `sol-dtime` physical seconds (whichever is specified).  Files are named `field1.tec`, `field2.tec`, …, `field{N}.tec`.

### Probes

Point probes record the temperature at a user-specified cell index with a time step of `dtime` seconds (independent of the field-write interval).  Output is a two-column ASCII file `(time, T)` in the `OUTPUT/` directory.

---
