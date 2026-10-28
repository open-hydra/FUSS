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
\Delta t_i = \frac{\text{VNN}}{\alpha(T_i)\,\sum_{d}\lvert\nabla\xi_d\rvert_i^{2}},
\qquad \alpha = \frac{k}{\rho\,c_p}
$$

where $\xi_d$ is the cell index along direction $d$ (the sum runs over the two or three directions of the mesh) and $\nabla\xi_d$ is row $d$ of the cell metric tensor $M$ (see [Grid Metrics](numerics.md#grid-metrics)). On an orthogonal cell $\lvert\nabla\xi_d\rvert = 1/\Delta x_d$ and the condition is the familiar $\alpha\,\Delta t \sum_d 1/\Delta x_d^2 \le \text{VNN}$; on a skewed cell $\lvert\nabla\xi_d\rvert$ is larger, $1/(\Delta x_d \sin\varphi)$ for a parallelogram whose edges meet at the angle $\varphi$. The diffusivity is taken at the current cell temperature, so with temperature-dependent properties the step follows the local $\alpha$.

On a uniform 2-D mesh, orthogonal or skewed, the largest eigenvalue of the discrete diffusion operator is $4\,\alpha \sum_d \lvert\nabla\xi_d\rvert^2$ (the mixed-derivative terms vanish in that mode). In 3-D it is the same on every cell but those skewed in all three pairs of directions: with $c_{12}, c_{13}, c_{23}$ the cosines between the rows $\nabla\xi_d$ of $M$, when $c_{12}^2 + c_{13}^2 + c_{23}^2 + 2c_{12}c_{13}c_{23} > 1$ (a rhombohedron with edges 115° apart), it exceeds that value by a factor below 9/8. The step takes the bound $f = (2 + m + 1/m)/4$ of that factor, $m + 1$ the largest eigenvalue of the matrix of the cosines ($m$ the largest root of $m^3 - s\,m - 2c_{12}c_{13}c_{23} = 0$, $s$ the sum of their squares), which is exact on rhombohedra; $f = 1$ on every other cell, and on every 2-D mesh. The sum of the step is multiplied by $f$ in both passes (below); in the face pass, where the $k_d$ differ by direction, $f$ multiplies the whole sum, exact when the $k_d$ agree and a frozen-coefficient bound otherwise. VNN is therefore $\Delta t\,\lambda_{\max}/4$, and the stability limit of forward Euler and RK2, whose stability interval on the negative real axis is $[-2, 0]$, is $\text{VNN} \le 1/2$. RK3 reaches $-2.51$, i.e. $\text{VNN} \le 0.628$.

**Convention.** VNN is defined as in MOSE, $\Delta t = \text{VNN}\,\rho c_p / \sum_d k_d \lvert\nabla\xi_d\rvert^2$, so that one value of `vnn` means one step in a coupled deck; in 3-D, on cells skewed in all three pairs of directions, the sum is multiplied by the factor $f \le 9/8$ above. Until October 2026 FUSS carried a factor 2 in the denominator (limits 1 for Euler and RK2, 1.25 for RK3): **VNN$_{\text{new}}$ = VNN$_{\text{old}}$ / 2**. A deck written before this change runs at twice its step unless its `vnn` is halved; every deck of the repository was halved in the same change, with no change in any result (halving is exact in binary). The default is 0.15 (it was 0.3).

Until October 2026 the sum used the mean cell edge lengths, $\sum_d 1/\Delta l_d^2$, which on a skewed cell overstates the stable step by $1/\sin^2\varphi$: on a mesh sheared by 45° the explicit runs diverged above $\text{VNN} = 0.25$, and now above 0.5.

**Face conductivities.** The formula above takes the conductivity of the cell, but the flux of each face carries the conductivity of that face (see [Face Conductivity](numerics.md#face-conductivity)). The step therefore takes, in each direction $d$, the larger of the conductivities of the cell's two faces in that direction, each evaluated as its flux evaluates it:

$$
\Delta t_i = \frac{\text{VNN}\,\rho c_p}{\sum_{d} k_d\,\lvert\nabla\xi_d\rvert_i^{2}},
\qquad k_d = \max\bigl(k_{d-1/2},\,k_{d+1/2}\bigr).
$$

An interior face takes $k$ at the mean temperature of its two cells, or between two materials the series conductivity of its two half cells, which can reach twice the conductivity of the less conductive cell; a block-connection or chimera face takes the harmonic mean $2k_Pk_N/(k_P+k_N)$ of the two sides (the neighbour's from the first ghost cell); walls, symmetry and coupling faces take the cell's own $k$, as the wall fluxes do. Where every face carries the cell's own $k$ the formula is the cell formula above, and the cell keeps that step to the last bit; a block of one material with a constant $k$ column, connected only to that material, is not even visited, so such meshes cost nothing more. Elsewhere the second pass costs up to nine property lookups per cell and iteration (seven in 2-D).

Taking the larger of the two face conductivities is conservative: the row sum of the explicit operator (Gershgorin's bound, which the cell formula meets with equality on a uniform mesh) needs only their mean.

Until October 2026 the step used the cell conductivity everywhere, and the face of a material jump inside a block took the conductivity of the lower-index cell: at a 300 → 15 W/(m K) jump the first cell of the low-conductivity material exchanged heat through 20 times its own $k$, and the four-layer liner of the FAS stress campaign left the properties table at VNN = 0.25 within 5 iterations.

With `vnn-rise-threshold = N` the step is scaled by $n/N$ during the first $N$ iterations $n$.

!!! tip "Choosing VNN"

    The default is 0.15. Without IRS stay below the limit of the scheme (0.5 for euler and RK2, 0.628 for RK3), with a margin where the mesh varies quickly. With IRS, steady runs can use larger values: the test cases run at 1 with `irs-beta = 0.5`; the [admissible-update check](#admissible-updates) stops the run on an overshoot out of the properties table.

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

The smoothing amplifies the stable VNN limit by approximately $(1 + 2\varepsilon)$ per direction.  With $\varepsilon = 0.5$ (a common choice) the effective limit is $(1 + 1)^3 \approx 8$ times the un-smoothed value in 3-D. This is the estimate behind the decks of the suite that run IRS at VNN = 1; the limit of a given case is not computed, and the admissible-update check below catches an overshoot.

### Admissible updates

The smoothed update of a cell mixes the updates of its neighbours. Above the explicit limit it can overshoot the temperatures imposed at the boundaries and leave the range of the properties table, where the property lookups would read past its ends. Every update, smoothed or explicit, is therefore checked after each stage: the new temperature must be finite and lie in $[T_\text{min}, T_\text{max})$, the first and last temperature of `properties.dat`. A cell that fails stops the run: the rank that finds it writes an error message naming the block, the cell, the temperature and the table range, and ends with exit code 1 (the check runs inside an OpenMP loop, where no MPI call is allowed: the launcher ends the other ranks); nothing is clamped or replaced. This is MOSE's rule (`Check_And_Fix_State`). Until October 2026 a smoothed update that failed was halved until admissible (at most 10 times) and otherwise dropped, with an `IRS:` count printed every `shell-diter` iterations; the fallback had been introduced for a stiff multimaterial mesh whose runs it kept in a bounded residual orbit, a symptom of the face-conductivity defect fixed on 2026-10-14 (see [FAS cycles](../vv/fas-cycles.md#stress-campaign)), and no case of the test suite or of the FAS acceptance set ever used it (zero fallbacks on every row of `fuss-stress/ref/237c9c2/fas_vv.md`), so it was removed with the sync of October 2026 (decision D2). The two guards of the FAS cycles stay: a prolongated correction that leaves the table is halved or dropped, and a coarse-level smoothing update under a forcing is halved (`FAS:` reports), because a forced coarse problem is not a physical time evolution ([FAS cycles](#fas-cycles)).

### Input settings

```ini
[FUSS-Numerics]
irs-beta = 0.5
```

`irs-beta > 0` switches IRS on at setup; only then are its arrays allocated. During the run a new `irs-beta` takes effect at the next re-read of `input.ini` (every `ini-diter` iterations), and `irs-beta = 0` switches the smoothing off; a positive value in a run started without IRS is ignored with a warning. The default, 0, is no smoothing. Until October 2026 a separate key `irs = true` switched IRS on as well, with `irs-beta = 0` a smoothing that changed nothing; it was removed.

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
