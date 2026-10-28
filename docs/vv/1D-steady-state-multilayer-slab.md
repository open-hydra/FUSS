# 1-D Steady State – Multilayer Slab

Three material layers in series inside one block. With $k$ constant in each material the steady field is linear in each layer and carries the same heat flux through all of them; a face between two materials takes the series conductivity of its two half cells (see [Face Conductivity](../theory/numerics.md#face-conductivity)), so the finite-volume scheme reproduces that field exactly, whatever the cell widths on the two sides of an interface.

---

## Problem setup

One block, 24 × 2 cells, through-thickness $x \in [0, 0.015]$ m:

| Layer | Thickness (mm) | Cells | Cell width (mm) | $k$ (W/mK) | ρ (kg/m³) | c<sub>p</sub> (J/kgK) |
|---|---|---|---|---|---|---|
| coat | 2.5 | 8 | 0.3125 | 5 | 6000 | 500 |
| liner | 7.5 | 6 | 1.25 | 300 | 8900 | 390 |
| jacket | 5.0 | 10 | 0.5 | 15 | 7900 | 480 |

Every interface joins two materials and two cell widths (ratios 4 and 2.5), in both orientations of the jump (5 → 300 and 300 → 15 along $x$). Walls 302 at $x = 0$ (1000 K) and $x = 0.015$ m (400 K); the $y$ and $z$ faces are adiabatic. The exact solution has the heat flux

$$
q = \frac{T_\text{left} - T_\text{right}}{\sum_i L_i / k_i} = 699\,029 \ \text{W/m}^2
$$

and is linear in each layer. `make_inputs.py` writes `MESH/mesh.tec` and `INPUT/`.

## Numerical setup

| Parameter | Value |
|---|---|
| Time scheme | RK3, VNN 0.25, local time step |
| Integration variables | cons |
| Initial condition | 700 K |
| Residual threshold | 1 × 10⁻¹² (reached at iteration 81 800) |

## Results

| Quantity | Error | Tolerance |
|---|---|---|
| Field, max \|T − T<sub>exact</sub>\| | 5.2 × 10⁻¹⁰ K | 10⁻⁶ K |
| Faces 1, 2: max \|q<sub>w</sub> − q\| / \|q\| | ≤ 2.3 × 10⁻¹² | 10⁻⁶ |

Before October 2026 the face between two materials took the conductivity of the lower-index cell only. On this case the run then left the properties table at the first jacket cell ($i = 15$, $T = -128$ K) at VNN 0.25; with the face-conductivity time step it converged, but to a field 37.5 K off and a wall flux 11 % off.

```bash
cd test/steady-state/multilayer_Slab
python3 make_inputs.py && ./FUSS.sh solve && python3 verify.py
```
