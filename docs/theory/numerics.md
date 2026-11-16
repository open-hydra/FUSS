# Spatial Discretisation

FUSS uses a **cell-centred finite volume method** on structured multi-block hexahedral grids. This
page describes the heat flux at interior faces, block connections and boundary faces, the ghost
cells, and the grid metrics.

---

## Finite Volume Framework

All variables are stored at **cell centres**. The semi-discrete update of cell $P$ (derived in
[Governing Equations](governing-equations.md)) is

$$
\rho_P\,c_{p,P}\,\frac{\mathrm{d}T_P}{\mathrm{d}t}
= \frac{1}{V_P}\sum_{f=1}^{6} Q_f\,A_f
\;+\;\dot{q}_{\text{vol},P}
$$

where $Q_f = \bigl(k\,\nabla T\!\cdot\!\hat{\mathbf{n}}_\text{out}\bigr)_f$ is the heat flux entering
the cell through face $f$ and $A_f$ the face area. The code accumulates the residual
$R_P = -\sum_f Q_f\,A_f$, the heat leaving the cell. A forward-Euler stage adds
$\Delta t_P\,(-R_P/V_P + \dot q_{\text{vol},P})$ to $h_P$, or that divided by $\rho_P\,c_{p,P}$ to
$T_P$ (see [Time Integration](time-integration.md)).

Each stage computes the fluxes in three steps: ghost-cell fill, interior faces of each block,
boundary faces. On a 2-D grid (one cell in $k$) the gradient component along $k$ is dropped, and
faces 5 and 6 must be zero-flux faces (code `0` or `300`).

---

## Diffusive Flux at Interior Faces

Let $L$ and $R$ be the two cells of a face normal to the index direction $\xi_n$, and $\xi_t$ the two
other index directions. The face gradient is built from index differences: across the face for the
normal part, and as the mean of the central differences of the two cells for the tangential parts:

$$
\left.\frac{\partial T}{\partial \xi_n}\right|_f = T_R - T_L,
\qquad
\left.\frac{\partial T}{\partial \xi_t}\right|_f = \tfrac14\Bigl[\bigl(T_{L,t+1} - T_{L,t-1}\bigr) + \bigl(T_{R,t+1} - T_{R,t-1}\bigr)\Bigr].
$$

The face metric $M_f$, whose rows are $(\nabla\xi_d)_f$ (see [Grid Metrics](#grid-metrics)), turns
them into a Cartesian gradient:

$$
\nabla T_f = \sum_{d} \left.\frac{\partial T}{\partial \xi_d}\right|_f (\nabla\xi_d)_f,
\qquad
Q_f\,A_f = k_f\,A_f\,\nabla T_f\cdot\hat{\mathbf n}_f .
$$

$\hat{\mathbf n}_f$ points from $L$ to $R$: $Q_f A_f$ enters $L$ and leaves $R$. Next to a block
boundary the tangential differences read the ghost cells.

### Face metric (`diffusive-metric`)

The key `diffusive-metric` of `[FUSS-Numerics]` selects how $M_f$ is built. It is read at set-up
only.

| `diffusive-metric` | Interior faces | Block connections (101) | Walls (301–305) | Chimera (102) |
|:---:|---|---|---|---|
| `stencil` (default) | $M_f = J_f^{-1}$, $J_f$ from differences of the cell centres | one gradient, as at an interior face | $J$ from the same differences as the wall gradient | two-sided balance with the cell metrics |
| `mean` | $M_f = \tfrac12(M_L + M_R)$ | two-sided balance with the cell metrics | metric of the cell | two-sided balance |

**`stencil`.** The metric takes of the cell-centre positions $\mathbf{x}$ the same differences that
the gradient takes of $T$. They are the columns of the face Jacobian

$$
J_f = \bigl[D_a\mathbf{x},\; D_b\mathbf{x},\; D_c\mathbf{x}\bigr],
\qquad
M_f = J_f^{-1},
$$

with $D_a\mathbf{x} = \mathbf{x}_R - \mathbf{x}_L$ and $D_{b,c}\mathbf{x}$ the tangential differences
built like those of $T$. The centre of a cell is the mean of its eight nodes.

- **A linear field is exact on any mesh.** For $T = T_0 + \mathbf{g}\cdot\mathbf{x}$ every difference
  is $D T = J_f^{T}\mathbf{g}$, so the face gradient is $\mathbf{g}$: on skewed, stretched or folded
  meshes and across jumps of cell size (Vinokur, NASA CR-177416, 1986). The patch test of the suite
  checks this on 22 mesh families to $10^{-13}$ of $k|\mathbf{g}|A$
  ([Testing](../development/testing.md)).
- **The scheme is second order** on random, size-jump and sheared meshes (`ConvergenceOrder`).
- The column across the face, and the columns of index directions with one cell (the $k$ direction
  of a 2-D mesh), come from the mean of the two cell edges.
- A stencil that does not span the space (a degenerate face or cell) stops the run at set-up. The
  message names the block, the cell, the face and the three differences.

**`mean`.** The mean of the two cell metrics. In one dimension $M = 1/h$, so it divides $T_R - T_L$ by
the harmonic mean of $h_L$ and $h_R$ instead of the centre distance $\tfrac12(h_L + h_R)$. It
therefore **overstates the heat flux by $(1+r)^2/(4r) - 1$ at a cell-size ratio $r = h_R/h_L$**:
+4.2 % at $r = 1.5$, +8.9 % at $r = 1.8$. The bias does not decrease with refinement; the metric is
exact only where the two cells are identical, and first order or less on non-smooth meshes.

**Ghost-cell positions.** For the metric, each ghost cell sits where its value is right for a linear
field:

| Face | Ghost centre |
|---|---|
| block connection | the source cell |
| wall, chimera | $2\mathbf{x}_f - \mathbf{x}_1$ ($\mathbf{x}_f$ the face centre, $\mathbf{x}_1$ the cell) |
| symmetry, null | the mirror image of $\mathbf{x}_1$ |
| coupling (103, 104) | $\mathbf{x}_1$ |

The face metrics are stored at set-up: 18 values per cell in 2-D, 27 in 3-D.

---

## Face Conductivity

### Inside a block

Between two cells of the same material, $k$ is taken at the mean temperature of the two cells:

$$
k_f = k\bigl(\tfrac12(T_L + T_R)\bigr).
$$

Between two materials the two half cells are in series, each with its own material at its own
temperature:

$$
k_f = \frac{h_L + h_R}{\dfrac{h_L}{k_L(T_L)} + \dfrac{h_R}{k_R(T_R)}},
$$

where $h_L$, $h_R$ are the widths of the two cells normal to the face (the edge of each cell across
the face, projected on $\hat{\mathbf n}_f$).

- The face gradient divides $T_R - T_L$ by $\tfrac12(h_L + h_R)$, so the flux is that of two layers in
  series: the same one a block connection balances (below).
- For equal cells $k_f$ is the harmonic mean $2k_Lk_R/(k_L + k_R)$.
- In 1-D, with $k$ constant in each material, a field linear in each layer is exact
  ([Multilayer Slab](../vv/1D-steady-state-multilayer-slab.md)).
- A side with $k = 0$ makes the face adiabatic.

### Block connections (101) and chimera faces (102)

With `diffusive-metric = stencil` a block connection is an interior face: one gradient of the two
cells, the stencil metric, and the interior conductivity rule. The tangential neighbours on the other
side are the cells of the other block.

With `mean`, and at every chimera face, the two sides are two half cells in series. The normal flux
is written on each side with that side's metric, conductivity and tangential differences, as a
function of the interface temperature $T_\text{int}$. On the side of cell $P$:

$$
q_P(T_\text{int}) = k_P\Bigl[2\,(T_\text{int} - T_P)\,\nabla\xi_n^P\cdot\hat{\mathbf n}_f\,s + \sum_{t}\tfrac12\bigl(T_{P,t+1} - T_{P,t-1}\bigr)\,\nabla\xi_t^P\cdot\hat{\mathbf n}_f\Bigr],
$$

with $s = +1$ if the face is on the high-index side of $P$ and $-1$ otherwise; $q_N(T_\text{int})$ is
built the same way for the neighbour $N$. Both are linear in $T_\text{int}$, so $q_P = q_N$ gives
$T_\text{int}$ in closed form, and the face flux is $q_P(T_\text{int})$. In 1-D this is the series
resistance

$$
Q = \frac{T_N - T_P}{\dfrac{\Delta x_P}{2\,k_P} + \dfrac{\Delta x_N}{2\,k_N}}.
$$

For a block connection the neighbour is the cell of the connected block, copied into the ghost
cells, with its metric. For a chimera face it is the ghost cell, filled from the donors, with its own
extrapolated metric.

**One flux per connection face.** The two records of a connection face (one per block) are paired at
set-up. One of them owns the face, computes the flux once, and gives it with opposite signs to the
two cells; between MPI ranks one value per face travels from the owner. The two cells therefore
receive exactly opposite shares. A mesh split into blocks gives the fluxes of the unsplit mesh when
the index directions agree and the face nodes are identical (test `ConnectionSplit`).

---

## Boundary Faces

### Null and symmetry (0, 300)

No flux crosses the face. Ghost layer $g$ takes the temperature and material of interior cell $g$
counted from the face ($g = 1, 2$).

### Walls (301–305)

All wall types share one discretisation (`Lib_BC_Fluxes_Wall.f90`). At a wall face of cell $P$ the
face gradient takes its normal part from $T_w$ across half a cell, and its tangential parts from the
central differences of $P$ along the face. The heat flux into the solid is then linear in $T_w$:

$$
q_\text{in} = G\,(T_w - T_P) + C,
\qquad
G = 2\,k_P\,\nabla\xi_n\cdot\hat{\mathbf n}_f,
\qquad
C = k_P\sum_{t}\tfrac12\bigl(T_{P,t+1} - T_{P,t-1}\bigr)\,\nabla\xi_t\cdot\hat{\mathbf n}_\text{out}.
$$

- $k_P = k(T_P)$; $\nabla\xi_d$ are the rows of the wall-face metric $M_w$; $\hat{\mathbf n}_f$ points
  towards increasing $\xi_n$ and $\hat{\mathbf n}_\text{out}$ out of the solid.
- $G$ is the conductance of the half cell: $2k_P/h_n$ on a parallelepiped, $h_n$ the cell height
  normal to the face.
- $C$ is the cross-diffusion term. It vanishes on orthogonal cells and makes a linear field exact on
  skewed ones.
- With `stencil`, $M_w = J_w^{-1}$ with $J_w = [\,2s(\mathbf{x}_f - \mathbf{x}_P),\; (\mathbf{x}_{hi} - \mathbf{x}_{lo})/(hi - lo)\,]$:
  the same position differences the gradient takes of the temperature. A linear field then gives
  $T_w = T(\mathbf{x}_f)$ and the exact flux on any mesh. With `mean`, $M_w$ is the cell metric.

| Code | Prescribed | $T_w$ | $q_\text{in}$ |
|:----:|------------|-------|---------------|
| `302` | $T_w$ | given | $G\,(T_w - T_P) + C$ |
| `301` | $q_w$ | $T_P + (q_w - C)/G$ | $q_w$ |
| `303` | $h_\text{conv}$, $q_\text{rad}$, $T_\text{ref}$ | root of $h_\text{conv}(T_\text{ref} - T_w) + q_\text{rad} = G\,(T_w - T_P) + C$ | $h_\text{conv}(T_\text{ref} - T_w) + q_\text{rad}$ |
| `304` | $\varepsilon$, $T_\text{ref}$ | root of $\varepsilon\sigma(T_\text{ref}^4 - T_w^4) = G\,(T_w - T_P) + C$ | $\varepsilon\sigma(T_\text{ref}^4 - T_w^4)$ |
| `305` | $h_\text{conv}$, $\varepsilon$, $T_\text{ref}$ | root of $h_\text{conv}(T_\text{ref} - T_w) + \varepsilon\sigma(T_\text{ref}^4 - T_w^4) = G\,(T_w - T_P) + C$ | $h_\text{conv}(T_\text{ref} - T_w) + \varepsilon\sigma(T_\text{ref}^4 - T_w^4)$ |

For 301 and 303–305, Newton's method finds $T_w$ starting from $T_P$. One step is exact without
radiation; with radiation the iteration stops at $\lvert\Delta T_w\rvert \le 10^{-12}\,T_w$ (at most
50 steps). The cell receives $A_f\,q_\text{in}$.

**Ghost cells.** Ghost layer $g$ holds $2\,T_w - T_g$, with $T_g$ the interior cell $g$ counted from
the face. The wall temperature is then the mean of each mirrored pair, so for a linear field on a
uniform mesh the ghosts hold the exact values. The ghosts are filled first, with tangential
differences that stay inside the block (one-sided at the face edges). The flux evaluation then
recomputes $T_w$ with central differences that may read the ghosts of neighbouring faces. The two
values differ only at the face edges.

**Output.** `OUTPUT/wall.*` holds $T_w$ and $q_w = q_\text{in}$ for every wall type: the values the
residual applies to the written field, face edges included. A converged run therefore closes its
energy balance, $\sum q_w A + \sum \dot q_\text{vol} V$, to round-off. The output does not change the
run: the run is the same whatever the output cadence.

**Stiffness.** $T_w$ is solved at every evaluation, so
$\partial q_\text{in}/\partial T_P = -G\,\beta/(G + \beta)$ with
$\beta = h_\text{conv} + 4\varepsilon\sigma T_w^3$. This is smaller in magnitude than the $-G$ of a
prescribed temperature: convective and radiative walls do not make a boundary cell stiffer than a
temperature wall.

The [Skewed Plate](../vv/2D-steady-state-skewed-plate.md) checks the wall model with walls 302, 303,
301 and 305 on a sheared mesh.

### Coupling with hydra (103, 104)

The face flux is set by the hydra fluid solver; FUSS adds it to the residual unchanged. The ghost
cells are mirrored as for symmetry.

---

## Ghost Cells and Multi-Block Connectivity

Each block carries two layers of ghost cells, filled before the fluxes of every Runge–Kutta stage:

- **Block connection (101)**: copies of the temperature and material of the connected block, with
  the index permutation of the face. Between MPI ranks the values travel in non-blocking messages,
  posted before the local boundary work.
- **Walls, null, symmetry and coupling faces**: as described above.
- **Chimera (102)**: each ghost layer has its own donor list (`n1`, `n2` lines in `bc.txt`). The
  ghost value is the volume-weighted enthalpy of the donors, $h_g = \sum_c f_c\,h(T_c)$, converted
  back to $T_g$. All donors of a ghost must share one material.
- The tangential neighbours of connection and chimera faces, which may be ghost cells of a
  neighbouring block, are collected last.

The flux loops then run on each block independently. OpenMP threads share the cells of a block, and
MPI distributes whole blocks (see [Parallel Execution](../user/parallel.md)).

---

## Grid Metrics

### Cell metric

The metric of a cell comes from the matrix $A$ whose rows are its mean edge vectors
$\partial\mathbf{x}/\partial\xi$, $\partial\mathbf{x}/\partial\eta$, $\partial\mathbf{x}/\partial\zeta$
(each the mean of the four parallel edges):

$$
\begin{pmatrix}\Delta x \\ \Delta y \\ \Delta z\end{pmatrix}
= A^{T}
\begin{pmatrix}\Delta\xi \\ \Delta\eta \\ \Delta\zeta\end{pmatrix}
$$

Each cell stores $M = A^{-T}$. Its rows are $\nabla\xi$, $\nabla\eta$, $\nabla\zeta$: it maps index
differences of $T$ to the Cartesian gradient. The time step uses the cell metric (see
[Time Integration](time-integration.md#diffusive-stability-condition-vnn)); the faces use the metric
of [Diffusive Flux at Interior Faces](#diffusive-flux-at-interior-faces).

### Face areas and cell volumes

The face area vector $A_f\,\hat{\mathbf{n}}_f$ is half the cross product of the two face diagonals,
oriented towards increasing index. The cell volume is the exact volume of the trilinear cell:

$$
V = \frac{1}{6}\sum_{d=1}^{3} \mathbf{e}_d \cdot \left(\mathbf{S}_{d-} + \mathbf{S}_{d+}\right),
$$

with $\mathbf{e}_d$ the mean edge of direction $d$ and $\mathbf{S}_{d\mp}$ the area vectors of its two
faces, both along $+d$. It is exact on every cell, twisted faces included, and the cell volumes of a
block sum to the volume enclosed by its boundary faces. The formula uses only differences of nodes,
so its rounding is that of the cell size, not of the coordinates.

### Mesh checks at set-up

FUSS stops at set-up, with exit code 1 on every rank, on a mesh it cannot solve correctly. Each
message names the block and the cell.

- **Degenerate cell**: $\lvert\det A\rvert \le 10^{-12}\,\lvert A_1\rvert\lvert A_2\rvert\lvert A_3\rvert$
  ($\lvert A_d\rvert$ the length of row $d$), for every cell and every ghost cell.
- **Cell turned over**: a volume whose sign differs from that of $\det A$ at the cell centre.
- **Cell folded at a corner**, on every grid level: at a corner,
  $s\,\mathbf{e}_i\cdot(\mathbf{e}_j\times\mathbf{e}_k) \le 10^{-12}\,|\mathbf{e}_i||\mathbf{e}_j||\mathbf{e}_k|$,
  with $\mathbf{e}_{i,j,k}$ the three edges from that corner and $s$ the orientation of the block.
  Corners with an edge of zero length (on an axis) are skipped. In 2-D the corners decide for the
  whole cell; in 3-D they are a necessary condition.
- **Block one cell across**, on every grid level:
    - the block across a block connection must have at least two cells across its face, because the
      two ghost layers of the connection are its first two cells;
    - a block with one cell in a direction $d$ (other than $k$ of a 1-D or 2-D mesh) must have the
      lines of $d$ orthogonal to the others ($|\cos| \le 10^{-3}$). The face gradient takes no
      difference along a direction with one cell, so a sheared block would lose part of its flux.
      An orthogonal block one cell across runs.
- **Ghost layer turned over behind a chimera face**: the first ghost layer must keep the sign and
  size of the first cell ($\det A \cdot \det A_1 > 10^{-12}\det A_1^2$). The message gives the ratio
  of the first two cells across the face. The same mesh with a wall instead runs.

### Ghost metrics

- Block connection: copied from the connected block.
- Symmetry: mirrored, one layer at a time.
- Every other face, chimera included: computed from nodes extrapolated **linearly** along the index
  lines, $N_{-g} = 2N_{-g+1} - N_{-g+2}$.

Every ghost layer thus keeps the edge across the face of the first cell. Its edges along the face
are extrapolated too: where they change sharply across the first cell they can turn a layer over (on
a plane 2-D mesh, growth by 3 or more turns the first layer over, 5/3 the second). Only the first
ghost layer behind a chimera face is read by a flux; that is why only it is checked for sign.

---

## References

1. J. H. Ferziger, M. Perić, R. L. Street, *Computational Methods for Fluid Dynamics*, 4th ed., Springer, 2020.
2. F. Moukalled, L. Mangani, M. Darwish, *The Finite Volume Method in Computational Fluid Dynamics*, Springer, 2016.
3. R. Eymard, T. Gallouët, R. Herbin, "Finite volume methods," in *Handbook of Numerical Analysis*, Vol. VII, 2000.
4. M. Vinokur, *An analysis of finite-difference and finite-volume formulations of conservation laws*, NASA CR-177416, 1986.
