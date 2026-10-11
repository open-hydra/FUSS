# 2-D Steady State – Linear Patch Test on a Skewed Plate

A linear temperature field is reproduced exactly by the finite-volume scheme on any mesh of parallelograms, provided the boundary conditions are discretised consistently with the interior. The test imposes a linear field on a plate sheared by 45° through four different wall conditions and checks the solution and the wall output to round-off.

---

## Problem setup

One block of 20 × 20 cells, x sheared by y: the nodes are at $x = i/20 + y$, $y = j/20$ (0.01 m thick). The left and right faces are inclined at 45°, the bottom and top faces are horizontal. The exact solution is

$$T(x, y) = 1000 + 200\,x + 300\,y \quad \text{K},$$

with $k = 52$ W/mK, ρ = 7850 kg/m³, c<sub>p</sub> = 460 J/kgK. On each face the heat flux into the solid is $q = k\,\nabla T\cdot\hat{\mathbf n}_\text{out}$:

| Face | Condition | Data per cell | $q$ into the solid (W/m²) |
|---|---|---|---|
| 1 (left, inclined) | Temperature (302) | $T_w$ = exact T at the face centre | +3677.0 |
| 2 (right, inclined) | Convection (303) | h = 100 W/m²K, q<sub>rad</sub> = 0, $T_\text{ref} = T_w + q/h$ | −3677.0 |
| 3 (bottom) | Heat flux (301) | $q_w = -k\,\partial T/\partial y$ | −15 600 |
| 4 (top) | Convection + radiation (305) | h = 100 W/m²K, ε = 0.8, $T_\text{ref}$ from $h(T_\text{ref}-T_w) + \varepsilon\sigma(T_\text{ref}^4 - T_w^4) = q$ | +15 600 |

The wall data vary from cell to cell, so `make_inputs.py` writes `MESH/mesh.tec` and `INPUT/` directly instead of ATLAS.

## Numerical setup

| Parameter | Value |
|---|---|
| Time scheme | RK3, VNN 2.0, local time step |
| Integration variables | prim |
| Implicit residual smoothing | β = 0.5 |
| Initial condition | 1300 K |
| Residual threshold | 1 × 10⁻¹² (reached at iteration 3154) |

## Results

| Quantity | Error | Tolerance |
|---|---|---|
| Field, max \|T − T<sub>exact</sub>\| | 2.0 × 10⁻¹¹ K | 10⁻⁶ K |
| Faces 1–4, max \|T<sub>w</sub> − T<sub>exact</sub>\| | ≤ 2.0 × 10⁻¹¹ K | 10⁻⁶ K |
| Faces 1–4, max \|q<sub>w</sub> − q\| / \|q\| | ≤ 9.4 × 10⁻¹³ | 10⁻⁶ |

The test needs the three ingredients of the wall formulation (see [Boundary Conditions](../user/boundary-conditions.md)): the tangential part of the face gradient in every wall flux, the wall temperature solved with it, and ghost cells $2T_w - T$. Before October 2026 the walls ignored the tangential term (303–305), mirrored the ghost cells and had sign errors in BC 301 (odd faces) and in the tangential term of BC 302: the same case ended with a maximum error of 190 K. With only the sign errors fixed the error was still 4.6 K, at the cells next to the walls.

```bash
cd test/steady-state/skewed_Plate
./FUSS.sh solve && python3 verify.py
```
