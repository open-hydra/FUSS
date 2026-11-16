# 1-D Multilayer Slab

The test checks the face between two materials inside a block, with different cell widths on the two sides.

Setup and run: `test/steady-state/multilayer_Slab/README.md`.

---

## Problem

One block of 24 × 2 cells across a wall 15 mm thick:

| Layer | Thickness (mm) | Cells | Cell width (mm) | k (W/mK) | ρ (kg/m³) | c<sub>p</sub> (J/kgK) |
|---|---|---|---|---|---|---|
| coat | 2.5 | 8 | 0.3125 | 5 | 6000 | 500 |
| liner | 7.5 | 6 | 1.25 | 300 | 8900 | 390 |
| jacket | 5.0 | 10 | 0.5 | 15 | 7900 | 480 |

The walls are at 1000 K (x = 0) and 400 K (x = 15 mm). Each interface joins two materials and two
cell widths (ratios 4 and 2.5), with the jump in both directions: 5 → 300 and 300 → 15.

With k constant in each layer the field is linear in each layer and carries one heat flux:

$$
q = \frac{T_\text{left} - T_\text{right}}{\sum_i L_i / k_i} = 699\,029 \ \text{W/m}^2
$$

A face between two materials takes the series conductivity of its two half cells
([Face Conductivity](../theory/numerics.md#face-conductivity)). The scheme then reproduces this
field exactly, whatever the cell widths.

## Result

| Quantity | Error | Tolerance |
|---|---|---|
| Field, max \|T − T<sub>exact</sub>\| | 1.3 × 10⁻⁷ K | 10⁻⁶ K |
| Faces 1, 2: max \|q<sub>w</sub> − q\| / \|q\| | ≤ 5.5 × 10⁻¹⁰ | 10⁻⁶ |

The errors are those of the iteration at its stopping threshold.
