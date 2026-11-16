# Skewed Plate

The physics, the reference solution and the results are on the V&V page:
[Skewed Plate](../../../docs/vv/2D-steady-state-skewed-plate.md).

## Setup

| Parameter | Value |
|---|---|
| Mesh | one block, 20 × 20 cells, nodes at x = i/20 + y, y = j/20 (0.01 m thick) |
| Walls | 302, 303, 301, 305 on faces 1–4, data per cell from the exact field |
| Material | k = 52 W/mK, ρ = 7850 kg/m³, c<sub>p</sub> = 460 J/kgK |
| Initial condition | 1300 K |
| Time scheme | RK3, `vnn = 1.0`, steady |
| Integration variables | `prim` |
| IRS | `irs-beta = 0.5` |
| Stop | `res-threshold = 1e-12` (iteration 2 471) |

The wall data vary from cell to cell, so `make_inputs.py` writes `MESH/mesh.tec` and `INPUT/` instead of ATLAS; the files ship with the case.

## Check

`verify.py` passes when the field and the wall temperatures are within 10⁻⁶ K of the exact field and the wall fluxes within 10⁻⁶ relative.

## Run

```bash
cd test/steady-state/skewed_Plate
./FUSS.sh solve && python3 verify.py
```

CTest: `SkewedPlate` (labels fast, validation).
