# Time Integration

FUSS advances the semi-discrete heat equation with an explicit time scheme. This page describes the
schemes, the stability limit of the time step, the steady and time-accurate modes, implicit residual
smoothing, and multigrid: grid sequencing and FAS cycles.

---

## Time Schemes

`time-scheme` selects **euler** (default, forward Euler), **RK2** (Heun) or **RK3** (SSP,
Shu–Osher):

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

This is the `prim` form. With `integration-variables = cons` (the default) the stages advance the
volumetric enthalpy $h_i$ with the bracket alone, and $T_i$ is recovered from the $h(T)$ table after
each stage. Since $\mathrm{d}h/\mathrm{d}T = \rho\,c_p$, both forms integrate the same equation (see
[Governing Equations](governing-equations.md#state-vector)).

- The schemes are first, second and third order in time. The RK stages are convex combinations of
  forward-Euler steps (SSP).
- Boundary time tables are read at the stage time $t^n + c_s\,\Delta t$ ($c = 0$; $0, 1$;
  $0, 1, \tfrac12$), which keeps the temporal order.
- Ghost cells are refilled before every stage.

---

## Diffusive Stability Condition (VNN)

The explicit update is stable while the diffusion number of every cell stays below the limit of the
scheme. The time step of cell $i$ is

$$
\Delta t_i = \frac{\text{VNN}\,\rho c_p}{\sum_{d} k_d\,\lvert\nabla\xi_d\rvert_i^{2}},
$$

where $\nabla\xi_d$ is row $d$ of the cell metric (see [Grid Metrics](numerics.md#grid-metrics)) and
the sum runs over the directions of the mesh. $\rho$ and $c_p$ are taken at the cell temperature.

- **Orthogonal cells**: $\lvert\nabla\xi_d\rvert = 1/\Delta x_d$, and the condition is the familiar
  $\alpha\,\Delta t \sum_d 1/\Delta x_d^2 \le \text{VNN}$, $\alpha = k/(\rho c_p)$.
- **Skewed cells**: $\lvert\nabla\xi_d\rvert = 1/(\Delta x_d \sin\varphi)$ on a parallelogram whose
  edges meet at the angle $\varphi$, so the step shrinks with the skew.
- **Temperature-dependent properties**: the step follows the local values.

**Limits.** VNN is $\Delta t\,\lambda_{\max}/4$, with $\lambda_{\max}$ the largest eigenvalue of the
discrete diffusion operator. The stability limit is therefore

| Scheme | Limit |
|---|---|
| euler, RK2 | VNN ≤ 0.5 |
| RK3 | VNN ≤ 0.628 |

This is the convention of MOSE, so one value of `vnn` means one time step in a coupled deck. The
default is 0.15.

**Cells skewed in three directions.** On a uniform mesh $\lambda_{\max} = 4\,\alpha \sum_d
\lvert\nabla\xi_d\rvert^2$, except on 3-D cells skewed in all three pairs of directions. Let
$c_{12}, c_{13}, c_{23}$ be the cosines between the rows of $M$. When
$c_{12}^2 + c_{13}^2 + c_{23}^2 + 2c_{12}c_{13}c_{23} > 1$ (for example a rhombohedron with edges
115° apart), $\lambda_{\max}$ exceeds that value by a factor below 9/8. The step then multiplies its
sum by the bound $f = (2 + m + 1/m)/4$ of that factor, with $m + 1$ the largest eigenvalue of the
cosine matrix. $f$ is exact on rhombohedra and equal to 1 on every other cell and on every 2-D mesh.
The log prints its largest value as `Time-step skew factor (max)`.

**Face conductivities.** Each face flux carries the conductivity of its face (see
[Face Conductivity](numerics.md#face-conductivity)), so $k_d$ is the larger of the two face
conductivities of the cell in direction $d$, each evaluated as its flux evaluates it:

$$
k_d = \max\bigl(k_{d-1/2},\,k_{d+1/2}\bigr).
$$

| Face | Conductivity |
|---|---|
| interior, one material | $k$ at the mean temperature of the two cells |
| interior, two materials | the series conductivity of the two half cells (up to twice the lower $k$) |
| block connection, chimera | $2k_Pk_N/(k_P+k_N)$ |
| wall, symmetry, coupling | the cell's own $k$ |

Under pyrolysis this face pass is not made: the step takes $k_s/(\rho_s\,c_{p,s})$ of the cell (see [Time Step with Pyrolysis](#time-step-with-pyrolysis)).

Where every face carries the cell's own $k$ the step is the plain cell formula. A block of one
material with a constant $k$ column is not even visited, so such meshes cost nothing more. Elsewhere
the face pass costs up to nine property lookups per cell and iteration (seven in 2-D). Taking the
larger of the two faces is conservative: Gershgorin's bound needs only their mean.

**VNN ramp.** With `vnn-rise-threshold = N` the step is scaled by $n/N$ during the first $N$
iterations.

!!! tip "Choosing VNN"
    Without IRS stay below the limit of the scheme (0.5 for euler and RK2, 0.628 for RK3), with a
    margin where the mesh varies quickly. The test suite checks both limits on Cartesian and skewed
    3-D cells (`VnnStability`, [Testing](../development/testing.md)). With IRS, steady runs can go
    higher: the test cases run at 1 with `irs-beta = 0.5`. An overshoot out of the property table
    stops the run ([Admissible updates](#admissible-updates)).

---

## Operating Modes

### Time-accurate mode (`time-accurate = true`)

Every cell advances with the global minimum time step. The run goes on to `time-threshold` (or
`iter-threshold`), writing snapshots every `sol-dtime` or `sol-diter`. `res-threshold` does not end
it. Use this mode for every transient problem.

The residual is still printed and written to `residual-history.dat` every `res-diter` iterations: the Euclidean norm over all cells of the change of one iteration, $r = \sqrt{\sum (T - T^{\mathrm{old}})^2 + \sum \sum_r (\chi_r - \chi_r^{\mathrm{old}})^2 + \sum \left((p - p^{\mathrm{old}})/p_{\mathrm{ref}}\right)^2}$, where the reaction advancements `chi_r` count with pyrolysis and the gas pressure with the type-2 gas closure ($p_{\mathrm{ref}}$ = `[FUSS-Darcy] p-ref`, default 101325 Pa); without them it is the temperature change alone.

### Steady mode (`time-accurate = false`)

Each cell uses its own **local** time step. This pseudo-time marching reaches steady state much
faster on grids with widely varying cell sizes.

!!! warning "Local time stepping and transient accuracy"
    Local time stepping destroys temporal accuracy. Do not use it when the transient history
    matters.

**Convergence test.** The residual is the Euclidean norm, over all cells, of the temperature change
of one iteration, evaluated every `res-diter` iterations. A steady run ends when the residual has
fallen by the factor `res-threshold` from its value at the first iteration (of the grid level, with
multigrid). The default, 1e-10, is a fall of ten orders.

- **The test is relative, as in MOSE.** The norm grows with the temperature level and with the
  square root of the cell count, so an absolute threshold would stop the same problem at different
  depths on different meshes.
- **A residual that rises before it falls** ends the run when it is `res-threshold` times the first
  one: later than a test relative to its peak, never earlier.
- **Restart.** The first residual is stored with every solution (`FUSS_RES_REFERENCE`), so a restart
  ends where the uninterrupted run ends.
- `residual-history.dat` and the log report the residual as it is, not divided by the first one.

---

## Time Step with Pyrolysis

With `[FUSS-Pyrolysis]` the run is time-accurate, with `integration-variables = cons`, without IRS and on one grid level (the set-up refuses anything else). The Runge–Kutta stages advance the stored energy $E$ and the reaction advancements $\chi$ together; the type-2 gas closure adds a pressure solve before them. The model is that of [Pyrolysis](pyrolysis.md).

### Order of operations

One step from $t^n$ to $t^{n+1} = t^n + \Delta t$:

1. **Time step.** $\Delta t$ is the global minimum of the cell steps (below).
2. **Wall values.** The env tables of the ablating walls (306), and the exchange values of the coupled ones (307), are evaluated once, at $t^{n+1}$.
3. **Start-of-step state.** $T^n$, $\chi^n$ and, under type 2, $p^n$ and the stored gas $G^n$ are saved.
4. **Pressure (type 2).** The pressure equation is solved once for $p^{n+1}$, with coefficients from the start-of-step state and the production predicted at the mid-step temperature. The first step of a new run is run twice, the second time with the history the first pass supplies.
5. **Stages** $s = 1 \dots n_\text{RK}$, at the stage times of the scheme:
    1. ghost cells of $T$ and $\chi$;
    2. conduction fluxes with $k_s$, then the boundary fluxes; an ablating wall solves its surface energy balance here;
    3. **kinetics**: for every reaction, the closed-form increment $\delta\chi_r^{(s)} = \Phi(\chi_r^n;\ \kappa_r(T^{(s)}),\ \Delta t) - \chi_r^n$ from the start-of-step advancement over the full step, with $\kappa$ at the stage temperature;
    4. **gas enthalpy**: type 1 integrates the stage production along the outflow lines and adds the enthalpy flow; type 2 adds the enthalpy flow of the Darcy mass flux of $p^{n+1}$;
    5. **update**: $E$ and $\chi$ are advanced with the weights of the scheme, then $T$ is found by inverting $E$ at the new $\chi$, and the cached $\rho_s$, $c_{p,s}$, $k_s$ are refreshed.
6. **Gas-mass closure (type 2).** The difference between the gas the solid released in the stages and the gas the pressure solve expelled goes into $G^{n+1}$, and $T$ is inverted again at constant $E$.

**Kinetics in the stages.** Each stage applies its increment as a forward-Euler stage applies $\Delta t\,\mathcal L$, so with the coefficients of [Time Schemes](#time-schemes)

$$
\chi^{n+1} = \chi^n + \sum_s b_s\,\delta\chi^{(s)},
\qquad
b = (1);\ (\tfrac12, \tfrac12);\ (\tfrac16, \tfrac16, \tfrac23)
$$

for euler, RK2 and RK3. Each increment is exact for the rate frozen at its stage temperature, and the weights are positive and sum to 1, so a stiff reaction ($\kappa\Delta t \gg 1$) cannot overshoot $\chi = 1$. The stored energy is advanced by the stage formula of the scheme itself, with the stage increment $\Delta t\,(-R/V + \dot q_\text{vol})$.

**Pressure in the stages (type 2).** Every stage reads $p^{n+1}$ as frozen data. The state at which a stage inverts $E$ carries the share of the step the stage has reached: $p^{(s)} = p^n + c_s(p^{n+1} - p^n)$ and $G^{(s)} = G^n + c_s(G^{n+1} - G^n)$, with $c_s$ the time fraction of the stage's output (SSP RK3: 1, ½, 1). The first stage starts from $E^n$ itself.

### What limits the time step

The cell step keeps the VNN form of [Diffusive Stability Condition](#diffusive-stability-condition-vnn) with the cached properties of the cell:

$$
\Delta t_i = \frac{\text{VNN}}{\alpha_i\,f\sum_d \lvert\nabla\xi_d\rvert_i^2},
\qquad
\alpha_i = \frac{k_{s,i}}{\rho_{s,i}\,c_{p,s,i}} .
$$

A face between two materials needs no face pass. It takes the two half cells in series, so its term of the energy
norm of the discrete operator is

$$
C_f\,(T_L - T_R)^2 = \min_{T_f}\Bigl[D_L\,(T_L - T_f)^2 + D_R\,(T_f - T_R)^2\Bigr] \le D_L\,T_L^2 + D_R\,T_R^2,
\qquad D = \frac{2k_sA}{h},
$$

$D$ being the conductance of a half cell. The largest eigenvalue of the operator is therefore at most that of the two
materials apart, each with a fixed temperature on the interface, and the cell step of each side covers that problem:
the step is stable up to the nominal `vnn` (0.5 for euler and RK2, 0.628 for RK3) at a material jump as in one
material. At the nominal value itself the shortest mode is not damped, for one material or two (its amplification
factor has modulus 1), so a sharp start leaves an oscillation from cell to cell that does not decay: accuracy needs
a `vnn` below the limit. Within one material the face takes the mean of the two $k_s$, which has no such bound: across
a front sharp to one cell, a $k_s$ ratio above about 1.2 can exceed the nominal limit.

| Enters $\Delta t$ | Does not enter $\Delta t$ |
|--------------------|----------------------------|
| $k_s$, $\rho_s$, $c_{p,s}$ of the cell (current $T$ and $\chi$) | the face conductivities (no face pass; a face between two materials is bounded above) |
| the cell metric and the skew factor $f$ | the heat capacity of the stored gas (leaving it out only shortens $\Delta t$) |
| `vnn` and its ramp | the reaction rates: the kinetics are integrated in closed form |
| | the gas pressure: it is solved implicitly |
| | the gas advection: no Péclet or Courant limit is applied |
| | the surface energy balance of an ablating wall |

The time step does not resolve the reactions. A step that is long compared with $1/\kappa$ stays stable but loses accuracy, and under type 2 a mid-step prediction that crosses a threshold temperature can force the [end-of-step correction](pyrolysis.md#end-of-step-mass-correction) to draw on the solid.

### Residual

The residual of a pyrolysis run adds the change of the advancements, and under type 2 that of the pressure scaled by `p-ref` (see [Time-accurate mode](#time-accurate-mode-time-accurate-true)). It is reported only: a time-accurate run ends at `time-threshold` or `iter-threshold`.

---

## Implicit Residual Smoothing (IRS)

IRS smooths the residual with a Laplacian operator before the update. This raises the stability
limit of the explicit scheme, which pays off on steady problems on stretched meshes.

### Smoothing equation

In each direction $d$ the smoothed residual $\mathcal{L}^\ast$ solves

$$
\mathcal{L}^\ast_i =
\frac{\mathcal{L}_i
  + \varepsilon\,\bigl(\mathcal{L}^\ast_{i-1} + \mathcal{L}^\ast_{i+1}\bigr)}
{1 + 2\varepsilon}
$$

with $\varepsilon$ = `irs-beta`, approximated by **two Jacobi iterations**.

### Effect on stability

The smoothing raises the stable VNN by about $(1 + 2\varepsilon)$ per direction. With
$\varepsilon = 0.5$ that is about 8 times the plain limit in 3-D. This estimate is behind the decks
of the suite that run IRS at VNN = 1. The limit of a given case is not computed; the admissibility
check below catches an overshoot.

### Admissible updates

Above the explicit limit a smoothed update can overshoot the temperatures imposed at the boundaries
and leave the property table. Every update, smoothed or not, is therefore checked after each stage:
the new temperature must be finite and lie in $[T_\text{min}, T_\text{max})$, the range of
`properties.dat`.

A cell that fails stops the run. The rank that finds it names the block, the cell, the temperature
and the table range, and ends with exit code 1; the launcher ends the other ranks. Nothing is
clamped or replaced. This is MOSE's rule.

The FAS cycles keep two guards of their own, because a forced coarse problem is not a physical time
evolution ([FAS cycles](#fas-cycles)).

### Input

```ini
[FUSS-Numerics]
irs-beta = 0.5
```

- `irs-beta > 0` switches IRS on at set-up; the default, 0, is no smoothing.
- During the run `irs-beta` keeps its set-up value: it is not one of the keys a reading of `input.ini` may change (see [Changes During a Run](../user/input.md#changes-during-a-run)), and a changed value is listed in `OUTPUT/<prefix>runtime-ini-ignored.txt`, so IRS can be neither switched on nor off during a run.

---

## Multigrid (grid sequencing)

With `levels = n` (n > 1) in `[FUSS-Multigrid]`, FUSS builds n grids. Each coarser level merges
2 × 2 × 2 cells (2 × 2 in 2-D) of the finer one.

- The cell counts of every block must be divisible by $2^{n-1}$, and the coarsest level must keep at
  least 2 cells per direction.
- The run starts on the coarsest level, from the initial field restricted to it (volume average of
  $h$).
- Level $m$ runs `level<m>-iter` iterations, or until its residual has fallen by `res-threshold` from
  its first iteration. Its solution is then prolongated (trilinear interpolation) to the next finer
  level.
- The finest level ends at `level1-iter`, `res-threshold` or `time-threshold`. `iter-threshold` is
  not used.
- With `cycle = none` (the default) a coarse level is not visited again: this is plain grid
  sequencing.

`iter-threshold` (without multigrid) is re-read with `input.ini`: lowering it to or below the
iterations already done ends the run at the next reading. `level<m>-iter` keeps its set-up value; a
changed value is listed in `OUTPUT/<prefix>runtime-ini-ignored.txt`.

**Restart.** A restart reads a solution written at level 1 and goes on at level 1: the coarse levels
are restricted from it. A solution without the restart counters starts again from the coarsest level
(see [Restart](../user/initial-conditions.md#restart)).

### FAS cycles

`cycle = v | w` (default `none`) turns the levels below the running one into a FAS multigrid (Full
Approximation Scheme). It needs a steady run with `integration-variables = cons`. Each solver step
becomes one cycle rooted at the running level, and the grid-sequencing descent acts as the full
multigrid (FMG) shell.

One cycle:

1. smooths `pre-iter` iterations (default 2);
2. restricts the state (volume average of $h$) and the volume-integrated steady residuals
   $N = -R + \dot q V$ (sum of the children);
3. builds the coarse forcing $F = \sum N_\text{fine} - N_\text{coarse}(w_0)$ and adds it to every
   Runge–Kutta stage of the coarse level;
4. recurses to the next level;
5. prolongates the correction in $h$ with trilinear weights;
6. smooths `post-iter` iterations (default 1).

The coarsest level of a cycle runs `coarse-iter` iterations (default 10). When the root of the
descent is itself the coarsest level, it smooths max(`pre-iter` + `post-iter`, 1) iterations per
step instead.

**W convention.** `cycle = w` recurses twice at *every* level below the running one, the root
included: with $L$ levels the coarsest is visited $2^{L-1}$ times per cycle. The textbook W recurses
twice only below the root, so FUSS's W runs its coarse sub-tree twice as often.

- On two levels W is bit-identical to V with twice the coarse budget.
- On deep hierarchies the doubled sub-tree makes W stable at any depth without IRS and insensitive to
  `coarse-iter`. The extra cost is small, because each coarser level has 4–8× fewer cells.

**Stability.** W converges at grid-independent rates on every isotropic hierarchy tested. V
*without IRS* is reliable up to 3 levels. On 4 or more it can stall: the residual oscillates at a
large value and the log fills with `FAS: ... updates halved` reports. Use W, or enable IRS, which
stabilises V at any depth.

**Chimera.** Cycles work across chimera (102) boundaries without restriction. The FAS transfers never
read a receiver cell: the restriction, forcing and correction sums cover interior cells only. The
receivers are ghost cells, refreshed by each level's own chimera fill, and the forcing absorbs the
difference between coarse and fine donors. With IRS the residual ghost at a 102 face is
zero-gradient. V and W both converge on overlapping blocks, bit-identically on 1 and 2 ranks
([Chimera results](../vv/fas-cycles.md#chimera-overset-blocks)).

**Guards.**

- Set-up stops a cycle request with `levels = 1`, `time-accurate = true`,
  `integration-variables = prim` or `pre-iter + post-iter = 0`.
- At runtime the multigrid keys keep their set-up values: none of them is a key a reading of `input.ini` may change (see [Changes During a Run](../user/input.md#changes-during-a-run)), and a changed `cycle`, `pre-iter`, `post-iter` or `coarse-iter` is listed in `OUTPUT/<prefix>runtime-ini-ignored.txt`.
- A prolongated correction that leaves the property table is halved up to 10 times, else dropped
  (`FAS: n prolongated corrections halved or dropped`).
- A coarse smoothing update pushed out of the table by the forcing is halved the same way
  (`FAS: n coarse smoothing updates halved`).

Occasional reports during a transient are harmless; a steady stream means the cycle is not
converging.

**Counters.** Under cycles the fine-level iteration counter advances by `pre-iter` + `post-iter` per
step. `res-diter`, `shell-diter` and `sol-diter` are tested at the end of each cycle only, so an
interval that is not a multiple of the per-cycle count stretches. With the defaults the residual
history gets one entry per cycle. A run that reaches `level1-iter` inside a cycle finishes that cycle
and ends on its boundary, the state a restart goes on from. `res-threshold` and `time-threshold` end
the run at once.

**Results.** The cycles converge to the fine-grid solution: on non-smooth 3-level plates the field
equals the single-grid one to $2.2\times10^{-7}$ K (`MultigridSolution`,
[Testing](../development/testing.md)). Measured on the verification cases
([Multigrid Cycles](../vv/fas-cycles.md)):

- on a two-level plate, 14× (V) to 91× (W + IRS) fewer fine iterations than the single grid;
- on full hierarchies (3–5 levels), W + IRS converges in a mesh-independent number of cycles, up to
  three orders of magnitude fewer iterations at 128²;
- nonlinear $k(T)$, $c_p(T)$ and layered cases gain about 200× in iterations;
- on cells of aspect ratio 100 the cycles are no faster than IRS alone
  ([Anisotropic cells](../vv/fas-cycles.md#anisotropic-cells)).

**Coarse-level output.** Every field write of a coarse level also writes its solution prolongated to
the finest level (`field-prolongated`, `wall-prolongated`). On a large mesh that write can cost more
than the coarse iterations between two writes: `prolongated-diter` (its set-up value) spaces it out. The end of a level
always writes it.

---

## Solution Output

- **Steady runs**: with `sol-overwrite = true` one field file is written every `sol-diter`
  iterations, overwriting the previous one.
- **Transient runs**: with `sol-overwrite = false` a numbered file (`field1.tec`, `field2.tec`, …) is
  written every `sol-diter` iterations or `sol-dtime` seconds.
- **Probes**: the temperature at a cell, every `dtime` seconds, independent of the field writes. See
  [Probes](../user/using.md#probes).
