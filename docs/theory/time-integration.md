# Time Integration

FUSS advances the semi-discrete heat equation with an explicit time scheme.  This page describes the schemes, the diffusive stability condition that limits the time step, the two operating modes (steady-state and time-accurate), the optional implicit residual smoothing, and multigrid: grid sequencing and FAS cycles.

---

## Time Schemes

`time-scheme` selects **euler** (default, forward Euler), **RK2** (Heun) or **RK3** (SSP, Shu–Osher):

$$
\begin{aligned}
\textsf{euler}:\quad T^{n+1} &= T^n + \Delta t\;\mathcal{L}(T^n)\\[6pt]
\textsf{RK2}:\quad T^{(1)} &= T^n + \Delta t\;\mathcal{L}(T^n)\\
T^{n+1} &= \tfrac{1}{2}\,T^n + \tfrac{1}{2}\bigl(T^{(1)} + \Delta t\;\mathcal{L}(T^{(1)})\bigr)\\[6pt]
\textsf{RK3}:\quad T^{(1)} &= T^n + \Delta t\;\mathcal{L}(T^n)\\
T^{(2)} &= \tfrac{3}{4}\,T^n
  + \tfrac{1}{4}\bigl(T^{(1)} + \Delta t\;\mathcal{L}(T^{(1)})\bigr)\\
T^{n+1} &= \tfrac{1}{3}\,T^n
  + \tfrac{2}{3}\bigl(T^{(2)} + \Delta t\;\mathcal{L}(T^{(2)})\bigr)
\end{aligned}
$$

$\mathcal{L}(T)$ is the spatial residual of the heat equation:

$$
\mathcal{L}(T_i) =
\frac{1}{\rho_i\,c_{p,i}}\!\left[
\frac{1}{V_i}\sum_{f} Q_f\,A_f
\;+\;\dot{q}_{\text{vol},i}
\right]
$$

This is the `prim` form.  With `integration-variables = cons` (the default) the stages advance the volumetric enthalpy $h_i$ with the bracket alone, and $T_i$ is recovered from the $h(T)$ table after each stage; since $\mathrm{d}h/\mathrm{d}T = \rho\,c_p$ the two forms integrate the same equation (see [Governing Equations](governing-equations.md#state-vector)).

The schemes are first, second and third-order accurate in time; the RK stages are convex combinations of forward-Euler steps (SSP).  Boundary-condition time tables are read at the stage time $t^n + c_s\,\Delta t$ ($c = 0$; $0, 1$; $0, 1, \tfrac12$), which preserves the temporal order.  Ghost cells are refilled before every stage.

---

## Diffusive Stability Condition (VNN)

The explicit update is stable while the diffusion number of every cell stays below the limit of the scheme. The time step of cell $i$ is

$$
\Delta t_i = \frac{\text{VNN}}{2\,\alpha(T_i)\,\sum_{d}\lvert\nabla\xi_d\rvert_i^{2}},
\qquad \alpha = \frac{k}{\rho\,c_p}
$$

where $\xi_d$ is the cell index along direction $d$ (the sum runs over the two or three directions of the mesh) and $\nabla\xi_d$ is row $d$ of the cell metric tensor $M$ (see [Grid Metrics](numerics.md#grid-metrics)). On an orthogonal cell $\lvert\nabla\xi_d\rvert = 1/\Delta x_d$ and the condition is the familiar $\alpha\,\Delta t \sum_d 1/\Delta x_d^2 \le \text{VNN}/2$; on a skewed cell $\lvert\nabla\xi_d\rvert$ is larger, $1/(\Delta x_d \sin\varphi)$ for a parallelogram whose edges meet at the angle $\varphi$. The diffusivity is taken at the current cell temperature, so with temperature-dependent properties the step follows the local $\alpha$.

On a uniform mesh, orthogonal or skewed, the largest eigenvalue of the discrete diffusion operator is $4\,\alpha \sum_d \lvert\nabla\xi_d\rvert^2$ (the mixed-derivative terms vanish in that mode). VNN is therefore the fraction of the stability limit of forward Euler and RK2, whose stability interval on the negative real axis is $[-2, 0]$: they are stable for $\text{VNN} \le 1$. RK3 reaches $-2.51$, i.e. $\text{VNN} \le 1.25$.

Until October 2026 the sum used the mean cell edge lengths, $\sum_d 1/\Delta l_d^2$, which on a skewed cell overstates the stable step by $1/\sin^2\varphi$: on a mesh sheared by 45° the explicit runs diverged above $\text{VNN} = 0.5$, and now above 1.

**Material jumps inside a block (`dt-face-k`).** The step above uses the conductivity of the cell, but an interior face carries the conductivity of its flux (see [Face Conductivity](numerics.md#face-conductivity)): with `face-k = mean-T` that of the lower-index cell's material at the mean temperature of the two cells, with `cell-mean` the mean of the two cells' conductivities. Where a block holds a material jump, the first cells of the low-$k$ material therefore exchange heat through up to $k_\text{high}/k_\text{low}$ times their own conductivity (`mean-T`, cells above a high → low jump along the index direction) or $(1 + k_\text{high}/k_\text{low})/2$ times (`cell-mean`, both sides), and the cell step overstates their stable step by up to that ratio. On the four-layer liner of the FAS stress campaign ($k$ = 5 / 300 / 15 / 150 W/(m K) through the thickness in one 64 × 128 block, ratio 20 at the 300 → 15 jump; the `multimat_liner` test case is the same wall split into two blocks along it) the run at VNN = 0.5 leaves the properties table at the first jacket cell ($T = -9.9$ K) after 5 iterations; at VNN = 0.2 it runs, slowly.

With `dt-face-k = true` (default `false`) direction $d$ of the sum takes the larger conductivity of the cell's two faces in that direction, each evaluated as its flux evaluates it:

$$
\Delta t_i = \frac{\text{VNN}\,\rho c_p}{2\sum_{d} k_d\,\lvert\nabla\xi_d\rvert_i^{2}},
\qquad k_d = \max\bigl(k_{d-1/2},\,k_{d+1/2}\bigr).
$$

A block-connection or chimera face takes the harmonic mean $2k_Pk_N/(k_P+k_N)$ of the two sides, the one-dimensional conductivity of its flux (the neighbour's from the first ghost cell); walls, symmetry and coupling faces take the cell's own $k$, as the wall fluxes do. A cell whose $k_d$ all equal its own $k$ keeps the cell step unchanged, so on one material at constant $k$ the key changes no bit of a run (CTest `DtFaceKInertWalls`, `DtFaceKInertChimera`). The key adds a second pass over the cells, with up to nine property lookups per cell and iteration (seven in 2-D): +4 % per iteration on one material, +8 % on the liner. It is re-read at runtime.

On the liner the key gives back the whole RK3 range: the run is stable up to VNN = 1.25 (at 1.3 it diverges inside the 300 W/(m K) layer: the RK3 limit of the uniform mesh) and reaches $10^{-8}$ in 239 377 iterations at VNN = 1.2, against 1 259 270 at VNN = 0.2 without the key; the converged field agrees with the W-cycle reference to $4\times10^{-6}$ K. The same mismatch held the IRS and W + IRS runs of the [FAS stress campaign](../vv/fas-cycles.md) in their residual orbits: with the key, `cycle = w` reaches $10^{-8}$ in 1 746 iterations at VNN = 1.0 (7 222 at VNN = 0.2 without it), W + IRS in 2 882 with no admissibility fallback, and IRS alone runs at VNN = 2 without fallbacks. With `face-k = cell-mean` the last coat cell sees $(5 + 300)/2 = 30$ times its own conductivity: without the key the run leaves the table there within 5 iterations already at VNN = 0.2 (it needs 0.1), with it the liner converges in 233 645 iterations at VNN = 1.2, to the `cell-mean` W-cycle field within $4\times10^{-6}$ K.

Taking the larger of the two face conductivities is conservative. The row sum of the explicit operator (Gershgorin's bound, which the cell formula meets with equality on a uniform mesh) needs only their mean, which at the 300 → 15 jump would allow a 1.9 times longer step to the first jacket cell; on the liner the mean saves 3 % of the iterations (231 619 against 239 377), so the safer bound is kept.

With `vnn-rise-threshold = N` the step is scaled by $n/N$ during the first $N$ iterations $n$.

!!! tip "Choosing VNN"

    The default is 0.3. Without IRS stay below the limit of the scheme (1 for euler and RK2, 1.25 for RK3), with a margin where the mesh varies quickly. Blocks that hold a material jump need `dt-face-k = true` for these limits to apply. With IRS, steady runs can use larger values: the test cases run at 2 with `irs-beta = 0.5` (3 did not converge on a test plate); the [admissible-update check](#admissible-updates) catches overshoots out of the properties table.

---

## Operating Modes

### Time-accurate mode (`time-accurate = true`)

The global minimum time step is used at every iteration.  The solver advances to the physical time `time-threshold`, writing field snapshots at intervals controlled by `sol-dtime` or `sol-diter`.  With a boundary-condition time table the residual threshold does not end the run (the table keeps forcing the solution): only `time-threshold` or `iter-threshold` do.

This mode is required for all transient problems: oscillating boundary conditions, evolving temperature profiles, probe monitoring.

### Steady-state mode (`time-accurate = false`)

Each cell uses its **local** time step (the cell's own VNN limit without taking the global minimum).  This is equivalent to pseudo-time marching and accelerates convergence to steady state by large factors on grids with widely varying cell sizes.

Convergence is declared when the residual — the Euclidean norm over all cells of the temperature change of one iteration, evaluated every `res-diter` iterations — falls below `res-threshold` (see [Governing Equations](governing-equations.md#steady-state-limit)).

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

The smoothed update of a cell mixes the updates of its neighbours. Above the explicit limit it can overshoot the temperatures imposed at the boundaries and leave the range of the properties table, where the property lookups would read past its ends. Each smoothed update is therefore checked: the new temperature must be finite and lie in $[T_\text{min}, T_\text{max})$, the first and last temperature of `properties.dat`. A cell that fails halves its smoothed update until the new state is admissible (at most 10 times), and otherwise keeps no update in that stage. The halved update keeps the implicit damping: the earlier fallback to the cell's own *unsmoothed* update reintroduced exactly the explicit update IRS exists to stabilise, and on a stiff multimaterial mesh (large conductivity jumps, squeezed cells) enough cells fell back at once to lock the run into a bounded residual orbit. The number of halved or dropped cells, summed over the MPI ranks, is printed every `shell-diter` iterations. An admissible smoothed update is applied unchanged.

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

With `levels = n` (n > 1) in `[FUSS-Multigrid]`, FUSS builds n grids, each coarser level merging 2 × 2 × 2 cells (2 × 2 in 2-D) of the finer one, so the cell counts of every block must be divisible by $2^{n-1}$, and the coarsest level must keep at least 2 cells per direction (a 1-cell-wide coarse block has degenerate ghost metrics). The run starts on the coarsest level, from the initial field restricted to it (volume average of $h$). Level $m$ runs `level<m>-iter` iterations, or until the residual reaches `res-threshold`; its solution is then prolongated (trilinear interpolation) to the next finer level, where the run continues. The finest level ends after `level1-iter` iterations or at `res-threshold` / `time-threshold`: with multigrid, `iter-threshold` is not used. With `cycle = none` (the default) a coarse level is not visited again: the descent is plain grid sequencing.

`level<m>-iter`, and `iter-threshold` without multigrid, are re-read with `input.ini` every `ini-diter` iterations: lowering the value of the level being computed below the iterations already done ends that level at the next re-read.

A restart (`newrun = false`) reads only the finest-level field and **redoes the descent**: the sequencing prolongation then overwrites the restart field on the finest level (exactly at interior cells, with one-sided-stencil errors in the boundary cell ring), and the finest level re-converges from there. Restarting is therefore safe but not free; a converged field is not preserved through a restart with `levels > 1`.

### FAS cycles

`cycle = v | w` (default `none`, which is bit-identical to the pre-FAS behaviour) turns the levels below the running one into a FAS multigrid (Full Approximation Scheme), for steady runs with `integration-variables = cons`: each solver step becomes one cycle rooted at the running level, so the grid-sequencing descent above acts as the FMG shell. A cycle smooths `pre-iter` iterations (default 2), restricts the state (volume average of $h$) and the volume-integrated steady residuals $N = -R + \dot q V$ (conservative sum of the 2 × 2 (× 2) children), builds the coarse forcing $F = \sum N_\text{fine} - N_\text{coarse}(w_0)$ and adds it to every Runge–Kutta stage, recurses, prolongates the correction in $h$ with the trilinear weights, and smooths `post-iter` iterations (default 1). The coarsest level of a cycle runs `coarse-iter` iterations (default 10); when the *root* of the descent is itself the coarsest level, that stint smooths max(`pre-iter` + `post-iter`, 1) iterations per step instead.

**W convention.** `cycle = w` applies the recursion count γ = 2 at *every* level below the running one, including at the root: with $L$ levels the coarsest is visited $2^{L-1}$ times per cycle, each intermediate visit re-restricting from its (already corrected) finer level. The textbook W applies γ only below the root, so this convention runs twice its coarse sub-tree. Measured consequences: on 2 levels it is *bit-identical* to a V-cycle with a doubled coarsest budget (the convention only relabels the work: W with `coarse-iter` c reproduces V with 2c to the byte), and on deep hierarchies the doubled sub-tree is what buys the measured robustness — stable at any depth without IRS and insensitive to `coarse-iter`, where V is not — at a per-cycle cost that stays small because each coarser level carries 4–8× fewer cells.

**Stability.** `w` converges with textbook, grid-independent rates on every verified isotropic hierarchy and is insensitive to `coarse-iter` there (on a 2-level hierarchy the coarse solve is the knob). `v` *without IRS* is reliable up to 3 levels; on 4 or more it can lock into a bounded limit cycle (the residual oscillates at a large value and the log fills with `FAS: ... updates halved` reports) whatever `coarse-iter` and even from a converged FMG descent — the explicit pre/post smoothing is too weak for the accumulated coarse corrections. Use `w`, or enable IRS, which stabilises `v` at any depth.

**Chimera.** Cycles work unrestricted across chimera (102) boundaries. The FAS transfers never
read a receiver cell: restriction, forcing and correction sums cover interior cells only, the
receivers are ghost cells refreshed by each level's own chimera fill (so the coarse donors come
from the restricted donor field, which the forcing absorbs like any other coarse-operator
inconsistency), and with IRS the residual ghost at a 102 face is zero-gradient, the same policy
IRS applies to every non-connection ghost. Both `v` and `w` are verified on overlapping blocks
([chimera numbers](../vv/fas-cycles.md#chimera-overset-blocks)): convergence to the `none` field,
1-vs-2-rank bit-identity with donors crossing ranks on every level, no admissibility fallbacks.

**Guards.** Setup aborts a cycle request with `levels = 1`, `time-accurate = true`, `integration-variables = prim` or `pre-iter + post-iter = 0`. At runtime (`ini-diter` re-reads) switching cycles on from `none` is ignored with a warning (the coarse-level arrays are allocated at setup only), and an edit that breaks the setup guards restores the setup values; `v` ↔ `w` and valid `pre/post/coarse-iter` edits take effect. Two admissibility guards protect the properties table and report at the `shell-diter` cadence: a prolongated correction that leaves the table is halved up to 10 times, else dropped (`FAS: n prolongated corrections halved or dropped`), and a coarse smoothing update pushed outside the table by the forcing is halved the same way instead of stopping the run (`FAS: n coarse smoothing updates halved`). Occasional reports during a transient are harmless; a steady stream of them means the cycle is not converging (see Stability above).

Under cycles the root-level iteration counter advances by `pre-iter` + `post-iter` per solver step, and the `res-diter` / `shell-diter` / `sol-diter` conditions are only tested at the end of each cycle: an interval that is not a multiple of the per-cycle count is effectively stretched (up to that factor). With the defaults (`res-diter = shell-diter = 1`) the residual history and the shell get one entry per cycle, carrying the residual of the cycle's last fine iteration.

The cycles converge to the fine discrete solution: the forcing absorbs any inconsistency of the coarse operators and boundary values, so the converged field (and, on an exactly representable solution, the bit-exact fixed point) is the same with and without cycles. Measured on the verification cases ([V&V](../vv/fas-cycles.md)): with the default `coarse-iter`, 2-level cycles on the `residual_convergence` plate need 14× (`v`) to 103× (`w` + IRS) fewer fine-grid iterations than the single-grid run and 7.5–24× less wall time (the mesh caps the hierarchy at 2 levels); with full hierarchies (3–5 levels) `w` + IRS converges in a roughly $N$-independent number of cycles — up to three orders of magnitude fewer iterations than single-grid at $128^2$ — and nonlinear $k(T)$, $c_p(T)$ or layered-material cases gain ~200× in iterations. Strongly anisotropic cells (AR ≫ 1) remain the known weak spot of full coarsening with a pointwise smoother: at aspect ratio 100 the cycles are no faster than IRS alone, until a line-implicit smoother is available.

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
