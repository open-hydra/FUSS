# Spatial Discretisation

FUSS uses a **cell-centred finite volume method** (FVM) on structured multi-block hexahedral grids.  This page describes how the diffusive heat flux is discretised at cell faces, how conductivity is handled at material interfaces, and how the multi-block connectivity is assembled.

---

## Finite Volume Framework

All variables are stored at **cell centres**.  The semi-discrete update for cell $i$ (derived in [Governing Equations](governing-equations.md)) is

$$
\rho_i\,c_{p,i}\,\frac{\mathrm{d}T_i}{\mathrm{d}t}
= \frac{1}{V_i}\sum_{f=1}^{6} Q_f\,A_f
\;+\;\dot{q}_{\text{vol},i}
$$

where $Q_f = \bigl(k\,\nabla T\!\cdot\!\hat{\mathbf{n}}\bigr)_f$ is the face-normal heat flux and $A_f$ is the face area.  The sum runs over all six faces of the hexahedral cell.

---

## Diffusive Flux at Interior Faces

At a face shared by two cells $L$ (left) and $R$ (right) the normal heat flux is approximated with a **second-order central difference**:

$$
Q_f = k_f\,\frac{T_R - T_L}{\delta}
$$

where $\delta$ is the distance between the two cell centres projected onto the face normal $\hat{\mathbf{n}}_f$, and $k_f$ is the face-averaged conductivity.

!!! info "Physical coordinates"

    On a non-uniform mesh the distance $\delta$ and the face normal $\hat{\mathbf{n}}_f$ are computed from the cell metric tensor, which maps computational indices $(\xi, \eta, \zeta)$ to physical coordinates $(x, y, z)$.  The flux therefore remains second-order accurate on smoothly stretched grids.

---

## Face Conductivity

### Single-material faces

When both cells share the same material, the face conductivity is the **arithmetic mean**:

$$
k_f = \tfrac{1}{2}(k_L + k_R)
$$

For uniform $k$ within a block this reduces to $k_f = k$.

### Multi-material interfaces (block connections)

At a block-to-block interface where $k_L \neq k_R$, the **harmonic mean** is used instead:

$$
k_f = \frac{2\,k_L\,k_R}{k_L + k_R}
$$

This choice is physically motivated: in steady one-dimensional conduction the heat flux through two resistors in series is

$$
Q = \frac{k_L k_R}{\tfrac{1}{2}(k_L \Delta x_L + k_R \Delta x_R)}\,(T_L - T_R)
$$

For equal cell spacings ($\Delta x_L = \Delta x_R$) this exactly reproduces the harmonic-mean conductivity, ensuring that the discretisation is consistent with the analytical interface condition (continuity of $T$ and $k\,\partial T / \partial n$).

!!! warning "High conductivity ratios"

    When $k_R / k_L \gg 1$ the harmonic mean is dominated by the low-conductivity side ($k_f \approx 2 k_L$ as $k_R \to \infty$).  This is the correct physical behaviour: a perfect conductor adds no thermal resistance.  For near-insulating materials ($k \to 0$), the face flux approaches zero regardless of the temperature difference, which is also correct.

---

## Boundary Face Fluxes

The flux at a boundary face depends on which wall behaviour the face carries (see [Governing Equations](governing-equations.md) for the parameter combinations that select each behaviour and the internal numeric codes):

| Behaviour | Internal code | Face flux $Q_f$ |
|-----------|:-------------:|------------------|
| Symmetry / `null` | `300` / inactive | $0$ |
| Prescribed flux | `301` | $q_w$ (set directly; `q = 0` gives an adiabatic wall) |
| Prescribed temperature | `302` | $k_\text{boundary}\,(T_\text{wall} - T_i) / \delta_b$, where $\delta_b$ is half the cell width |
| Convection (+ $q_\text{rad}$) | `303` | $h_\text{conv}\,(T_\text{ref} - T_\text{wall}) + q_\text{rad}$ |
| Radiation | `304` | $\varepsilon\,\sigma\,(T_\text{ref}^4 - T_\text{wall}^4)$ |
| Convection + radiation | `305` | $h_\text{conv}\,(T_\text{ref} - T_\text{wall}) + \varepsilon\,\sigma\,(T_\text{ref}^4 - T_\text{wall}^4)$ |

For the Dirichlet BC (code 302) the ghost-cell value is set to $T_\text{ghost} = 2\,T_\text{wall} - T_i$ so that the central-difference formula evaluated at the face centre returns exactly $T_\text{wall}$. For the nonlinear convective and radiative BCs the wall temperature $T_\text{wall}$ is recovered from the local energy balance between conduction in the solid and the boundary flux, then used in the ghost-cell update.

---

## Multi-Block Connectivity

FUSS treats the full computational domain as a collection of structured blocks that share faces.  The assembly follows these steps:

1. **Ghost-cell layer**: each block is surrounded by one layer of ghost cells that store the values needed for face-flux computation across block boundaries.
2. **Block-connection interface update**: before each residual evaluation, ghost cells at connection interfaces are filled with the interior values of the adjacent block (with appropriate index permutation for non-aligned faces).
3. **Independent residual evaluation**: the flux loop runs independently on each block; the ghost-cell fill ensures the correct cross-block flux without explicit coupling in the residual assembly.

This approach maps naturally onto shared-memory and distributed-memory parallel architectures, since block residuals are independent given the ghost-cell data.

---

## Grid Metrics

Each cell face stores a precomputed **metric tensor** $M_{3\times 3}$ that maps computational coordinate increments $(\Delta\xi, \Delta\eta, \Delta\zeta)$ to physical increments $(\Delta x, \Delta y, \Delta z)$:

$$
\begin{pmatrix}\Delta x \\ \Delta y \\ \Delta z\end{pmatrix}
= M
\begin{pmatrix}\Delta\xi \\ \Delta\eta \\ \Delta\zeta\end{pmatrix}
$$

The face area vector $A_f\,\hat{\mathbf{n}}_f$ is obtained from the cross product of two edge vectors of the face, evaluated at the face centre from the nodal coordinates.  Cell volumes are computed by the divergence theorem applied to the six face area vectors.

---

## References

1. J. H. Ferziger, M. Perić, R. L. Street, *Computational Methods for Fluid Dynamics*, 4th ed., Springer, 2020 (chapters on FVM for diffusion).
2. F. Moukalled, L. Mangani, M. Darwish, *The Finite Volume Method in Computational Fluid Dynamics*, Springer, 2016.
3. R. Eymard, T. Gallouët, R. Herbin, "Finite volume methods," in *Handbook of Numerical Analysis*, Vol. VII, 2000.
