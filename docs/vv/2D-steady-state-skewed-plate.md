# 2-D Skewed Plate (Linear Patch Test)

A finite-volume scheme should reproduce a linear field exactly on any mesh, provided the walls are discretised consistently with the interior. This test checks the wall model on a
skewed mesh.

Setup and run: `test/steady-state/skewed_Plate/README.md`.

---

## Problem

One block of 20 × 20 cells, sheared by 45°: the nodes are at $x = i/20 + y$, $y = j/20$. The left
and right faces are inclined; the bottom and top are horizontal. The exact solution is

$$T(x, y) = 1000 + 200\,x + 300\,y \quad \text{K},$$

with k = 52 W/(m·K). Each wall gets the data that make this field exact. The heat flux into the
solid is $q = k\,\nabla T\cdot\hat{\mathbf n}_\text{out}$:

| Face | Condition | Data per cell | q into the solid (W/m²) |
|---|---|---|---|
| 1 (left, inclined) | temperature (302) | T<sub>w</sub> = exact T at the face centre | +3677.0 |
| 2 (right, inclined) | convection (303) | h = 100 W/(m²·K), q<sub>rad</sub> = 0, T<sub>ref</sub> = T<sub>w</sub> + q/h | −3677.0 |
| 3 (bottom) | heat flux (301) | q<sub>w</sub> = −k ∂T/∂y | −15 600 |
| 4 (top) | convection + radiation (305) | h = 100 W/(m²·K), ε = 0.8, T<sub>ref</sub> solving the wall balance | +15 600 |

## Result

| Quantity | Error | Tolerance |
|---|---|---|
| Field, max \|T − T<sub>exact</sub>\| | 5.6 × 10⁻⁹ K | 10⁻⁶ K |
| Faces 1–4, max \|T<sub>w</sub> − T<sub>exact</sub>\| | ≤ 5.6 × 10⁻⁹ K | 10⁻⁶ K |
| Faces 1–4, max \|q<sub>w</sub> − q\| / \|q\| | ≤ 2.2 × 10⁻¹⁰ | 10⁻⁶ |

The errors are those of the iteration, not of the scheme: they fall with the stopping threshold.

The result rests on three parts of the wall formulation (see
[Boundary Conditions](../user/boundary-conditions.md)):

- the tangential part of the face gradient in every wall flux;
- the wall temperature solved with that tangential part;
- ghost cells set to $2T_w - T$.

The face and wall metrics are built from the same differences of the cell centres as the gradient
(`diffusive-metric = stencil`). A linear field is therefore exact on any mesh, not only on
parallelograms. The `PatchTest` of the test suite checks this on random, folded, symmetric and 3-D
meshes ([Testing](../development/testing.md)).
