# Spatial Discretisation

FUSS uses a **cell-centred finite volume method** (FVM) on structured multi-block hexahedral grids.  This page describes how the heat flux is discretised at interior faces, block connections and boundary faces, how the ghost cells are filled, and how the grid metrics are computed.

---

## Finite Volume Framework

All variables are stored at **cell centres**.  The semi-discrete update for cell $P$ (derived in [Governing Equations](governing-equations.md)) is

$$
\rho_P\,c_{p,P}\,\frac{\mathrm{d}T_P}{\mathrm{d}t}
= \frac{1}{V_P}\sum_{f=1}^{6} Q_f\,A_f
\;+\;\dot{q}_{\text{vol},P}
$$

where $Q_f = \bigl(k\,\nabla T\!\cdot\!\hat{\mathbf{n}}_\text{out}\bigr)_f$ is the heat flux entering the cell through face $f$ and $A_f$ is the face area.  The code accumulates the residual $R_P = -\sum_f Q_f\,A_f$ (the heat leaving the cell); the forward-Euler increment of a stage is $\Delta t_P\,(-R_P/V_P + \dot q_{\text{vol},P})$ on $h_P$, or that divided by $\rho_P\,c_{p,P}$ on $T_P$ (see [Time Integration](time-integration.md)).

The fluxes of a stage are computed in three steps: ghost-cell fill, interior faces of each block, boundary faces.  On a 2-D grid (one cell in $k$, detected on the first block) the gradient component along $k$ is dropped; faces 5 and 6 must then be zero-flux faces (code `0` or `300`).

---

## Diffusive Flux at Interior Faces

Let $L$ and $R$ be the two cells of an interior face normal to the index direction $\xi_n$, and $\xi_t$ the two other index directions.  The face gradient is built from index differences: the normal one across the face, the tangential ones as the mean of the central differences of the two cells,

$$
\left.\frac{\partial T}{\partial \xi_n}\right|_f = T_R - T_L,
\qquad
\left.\frac{\partial T}{\partial \xi_t}\right|_f = \tfrac14\Bigl[\bigl(T_{L,t+1} - T_{L,t-1}\bigr) + \bigl(T_{R,t+1} - T_{R,t-1}\bigr)\Bigr],
$$

and converted to a Cartesian gradient with the face metric $M_f$, whose rows are $(\nabla\xi_d)_f$ (see [Grid Metrics](#grid-metrics)):

$$
\nabla T_f = \sum_{d} \left.\frac{\partial T}{\partial \xi_d}\right|_f (\nabla\xi_d)_f,
\qquad
Q_f\,A_f = k_f\,A_f\,\nabla T_f\cdot\hat{\mathbf n}_f .
$$

$\hat{\mathbf n}_f$ points from $L$ to $R$: $Q_f A_f$ enters $L$ and leaves $R$.  The tangential differences of the cells next to a block boundary read the ghost cells.

!!! info "Face metric"

    How $M_f$ is built from the metrics of the two cells is selected by `diffusive-metric` in `[FUSS-Numerics]`:

    | `diffusive-metric` | Face metric | |
    |---|---|---|
    | `inverse-mean` (default) | $M_f = \left[\tfrac12\left(A_L + A_R\right)\right]^{-T}$ | the normal difference is divided by the centre-to-centre distance: a linear field is reproduced exactly on stretched orthogonal meshes |
    | `mean` | $M_f = \tfrac12\left(M_L + M_R\right)$ | faster; biased on stretched meshes |

    In one dimension $M = 1/h$, so `mean` multiplies $T_R - T_L$ by the mean of $1/h_L$ and $1/h_R$, i.e. divides it by their harmonic mean instead of the centre-to-centre distance $\tfrac12(h_L + h_R)$. With a cell-size ratio $r = h_R/h_L$ it overstates the heat flux by $(1+r)^2/(4r)$: $+0.2\,\%$ at $r = 1.1$, $+0.8\,\%$ at $1.2$, $+4.2\,\%$ at $1.5$, $+8.9\,\%$ at $1.8$. The bias does not decrease with grid refinement at a fixed growth ratio and is present in the steady solution; `mean` is adequate where adjacent cells grow by about 1.2 or less. With `inverse-mean` FUSS stores $A^{T}$ per cell (9 values), so each face needs one $3\times 3$ inversion.

    The choice applies to interior faces only. Block-connection and wall faces use the metric of each cell over half a cell (below).

---

## Face Conductivity

### Inside a block

At an interior face the conductivity is evaluated at the mean temperature of the two cells,

$$
k_f = k\bigl(\tfrac12(T_L + T_R)\bigr),
$$

with the material of the lower-index cell $L$.  The material of $R$ is not used, so a material change inside a block would not be represented: each material must fill whole blocks.

### Block connections (101) and chimera faces (102)

At a block connection the two sides are treated as two half cells in series.  The normal flux $k\,\nabla T\cdot\hat{\mathbf n}_f$ is written on each side with that side's metric, conductivity (at its cell temperature, for its material) and tangential differences, as a function of the interface temperature $T_\text{int}$; on the side of cell $P$

$$
q_P(T_\text{int}) = k_P\Bigl[2\,(T_\text{int} - T_P)\,\nabla\xi_n^P\cdot\hat{\mathbf n}_f\,s + \sum_{t}\tfrac12\bigl(T_{P,t+1} - T_{P,t-1}\bigr)\,\nabla\xi_t^P\cdot\hat{\mathbf n}_f\Bigr],
$$

where $s = +1$ if the face is on the high-index side of $P$ and $-1$ otherwise; likewise $q_N(T_\text{int})$ for the neighbour cell $N$ beyond the face.  $T_\text{int}$ is the value at which $q_P = q_N$; since both are linear in $T_\text{int}$, it follows in closed form, and the face flux is $q_P(T_\text{int})$.  In one dimension this gives the series resistance

$$
Q = \frac{T_N - T_P}{\dfrac{\Delta x_P}{2\,k_P} + \dfrac{\Delta x_N}{2\,k_N}},
$$

which for $\Delta x_P = \Delta x_N = \Delta x$ is the harmonic mean $2k_Pk_N/(k_P + k_N)$ over $\Delta x$.  As $k_N \to \infty$ the flux tends to $2k_P(T_N - T_P)/\Delta x_P$; as $k_N \to 0$ it vanishes.

For a block connection the neighbour values are the cells of the connected block (copied into the ghost cells) and its metric; for a chimera face they are the ghost cell, filled from the donor cells, and the ghost cell's own (extrapolated) metric.

---

## Boundary Faces

### Null and symmetry (0, 300)

No flux crosses the face.  The ghost cells mirror the interior cells: ghost layer $g$ takes the temperature and material of interior cell $g$ counted from the face ($g = 1, 2$).

### Walls (301–305)

All wall types share one discretisation (`src/lib/numerics/fluxes/bc/Lib_BC_Fluxes_Wall.f90`).  At a wall face of cell $P$, normal to the index direction $\xi_n$, the face gradient takes its normal part from the wall temperature $T_w$ across half a cell and its tangential parts from the central differences of $P$ along the face.  The heat flux into the solid is then linear in $T_w$:

$$
q_\text{in} = G\,(T_w - T_P) + C,
\qquad
G = 2\,k_P\,\nabla\xi_n\cdot\hat{\mathbf n}_f,
\qquad
C = k_P\sum_{t}\tfrac12\bigl(T_{P,t+1} - T_{P,t-1}\bigr)\,\nabla\xi_t\cdot\hat{\mathbf n}_\text{out},
$$

where $k_P = k(T_P)$, $\nabla\xi_d$ are the rows of the cell metric $M_P$, $\hat{\mathbf n}_f$ is the face normal oriented towards increasing $\xi_n$ and $\hat{\mathbf n}_\text{out}$ the outward normal.  $G$ is the conductance of the half cell ($2k_P/h_n$, with $h_n$ the cell height normal to the face, on a parallelepiped cell).  $C$ is the cross-diffusion term: it vanishes on orthogonal cells and is needed for a linear field to be exact on skewed ones.

| Code | Prescribed | $T_w$ | $q_\text{in}$ |
|:----:|------------|-------|---------------|
| `302` | $T_w$ | given | $G\,(T_w - T_P) + C$ |
| `301` | $q_w$ | $T_P + (q_w - C)/G$ | $q_w$ |
| `303` | $h_\text{conv}$, $q_\text{rad}$, $T_\text{ref}$ | root of $h_\text{conv}(T_\text{ref} - T_w) + q_\text{rad} = G\,(T_w - T_P) + C$ | $h_\text{conv}(T_\text{ref} - T_w) + q_\text{rad}$ |
| `304` | $\varepsilon$, $T_\text{ref}$ | root of $\varepsilon\sigma(T_\text{ref}^4 - T_w^4) = G\,(T_w - T_P) + C$ | $\varepsilon\sigma(T_\text{ref}^4 - T_w^4)$ |
| `305` | $h_\text{conv}$, $\varepsilon$, $T_\text{ref}$ | root of $h_\text{conv}(T_\text{ref} - T_w) + \varepsilon\sigma(T_\text{ref}^4 - T_w^4) = G\,(T_w - T_P) + C$ | $h_\text{conv}(T_\text{ref} - T_w) + \varepsilon\sigma(T_\text{ref}^4 - T_w^4)$ |

For 301 and 303–305, $T_w$ is found by Newton iteration from $T_w = T_P$: one step is exact without radiation; with radiation the iteration stops when $\lvert\Delta T_w\rvert \le 10^{-12}\,T_w$ (at most 50 steps).  The cell then receives $A_f\,q_\text{in}$ ($R_P \mathrel{-}= A_f\,q_\text{in}$).

**Ghost cells.**  Ghost layer $g$ is set to $2\,T_w - T_g$, with $T_g$ the interior cell $g$ counted from the face: the wall temperature is the mean of each mirrored pair, so for a linear field on a uniform mesh the ghosts hold the exact values.  The ghosts are filled first, with $T_w$ computed from tangential differences that stay inside the block (one-sided at the face edges); the flux evaluation then recomputes $T_w$ with central differences that may read the ghost cells of the neighbouring faces.  The two values differ only at the face edges.

**Output.**  `OUTPUT/wall.tec` holds $T_w$ and $q_w = q_\text{in}$, the heat flux into the solid, for every wall type, computed with the same routine as the flux but from interior values only (the one-sided form above, with the time tables at the output time): the ghost cells hold the last stage of the step.

**Stiffness.**  Since $T_w$ is solved at every evaluation, $\partial q_\text{in}/\partial T_P = -G\,\beta/(G + \beta)$ with $\beta = h_\text{conv} + 4\varepsilon\sigma T_w^3$, smaller in magnitude than the $-G$ of a prescribed temperature (302): convective and radiative walls do not make the boundary cell stiffer than a temperature wall.

The linear patch test on a plate sheared by 45°, with walls 302, 303, 301 and 305 on its four sides, reproduces the exact field to $2\times10^{-11}$ K (see [Skewed Plate](../vv/2D-steady-state-skewed-plate.md)).

### Coupling with hydra (103, 104)

The face flux is the value set by the hydra fluid solver; FUSS adds it to the residual unchanged.  The ghost cells are mirrored as for symmetry.

---

## Ghost Cells and Multi-Block Connectivity

Each block carries two layers of ghost cells.  They are filled before the fluxes of every Runge–Kutta stage:

- **Block connection (101)**: the two layers are copies of the temperature and material of the connected block, with the index permutation of the face.  Between MPI ranks the values travel in non-blocking messages, posted before the local boundary work.
- **Walls (301–305)**, **null / symmetry** and the coupling faces: as described above.
- **Chimera (102)**: each of the two ghost layers has its own donor list (`n1`, `n2` lines in `bc.txt`); its value is the volume-fraction-weighted enthalpy of the donor cells, $h_g = \sum_c f_c\,h(T_c)$, converted back to $T_g$ with the donors' material (all donors of a ghost must share one material).
- The tangential neighbours of the connection and chimera faces, which may themselves be ghost cells of the neighbouring block, are collected after the steps above.

The flux loops then run on each block independently: interior faces read only the block's own cells and ghost cells; boundary faces add the connection, chimera, wall and coupling fluxes.  OpenMP threads share the cells of a block and MPI distributes whole blocks (see [Parallel Execution](../user/parallel.md)).

---

## Grid Metrics

The metrics of a cell are built from the matrix $A$ whose rows are its mean edge vectors $\partial\mathbf{x}/\partial\xi$, $\partial\mathbf{x}/\partial\eta$, $\partial\mathbf{x}/\partial\zeta$ (each the mean of the four parallel edges of the hexahedron), so that $A^{T}$ maps index increments to physical increments:

$$
\begin{pmatrix}\Delta x \\ \Delta y \\ \Delta z\end{pmatrix}
= A^{T}
\begin{pmatrix}\Delta\xi \\ \Delta\eta \\ \Delta\zeta\end{pmatrix}
$$

Each cell stores the metric tensor $M = A^{-T}$, computed from the cofactors of $A$ (and, with `inverse-mean`, $A^{T}$); the rows of $M$ are $\nabla\xi$, $\nabla\eta$, $\nabla\zeta$, and it maps index differences of the temperature to its Cartesian gradient.  At an interior face the metrics of the two cells are combined as described in [Diffusive Flux at Interior Faces](#diffusive-flux-at-interior-faces); the time step uses the cell metric (see [Time Integration](time-integration.md#diffusive-stability-condition-vnn)).

The face area vector $A_f\,\hat{\mathbf{n}}_f$ is half the cross product of the two diagonals of the face, oriented towards increasing index.  The cell volume is the sum of the volumes of five tetrahedra into which the hexahedron is split.  A degenerate cell (zero determinant of $A$) stops the run.

The metrics of the ghost cells are copied from the connected block for block connections, mirrored for symmetry faces, and computed from node coordinates extrapolated outwards for all other faces (chimera included).

---

## References

1. J. H. Ferziger, M. Perić, R. L. Street, *Computational Methods for Fluid Dynamics*, 4th ed., Springer, 2020 (chapters on FVM for diffusion).
2. F. Moukalled, L. Mangani, M. Darwish, *The Finite Volume Method in Computational Fluid Dynamics*, Springer, 2016.
3. R. Eymard, T. Gallouët, R. Herbin, "Finite volume methods," in *Handbook of Numerical Analysis*, Vol. VII, 2000.
