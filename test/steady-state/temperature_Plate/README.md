# Temperature Plate

The physics, the reference solution and the results are on the V&V page:
[Temperature Plate](../../../docs/vv/2D-steady-state-temperature-plate.md).

## Setup

| Parameter | Value |
|---|---|
| Mesh | 60 × 40 cells, 0.6 m × 0.4 m, uniform (Δx = Δy = 0.01 m) |
| Face 1 (x = 0) | `type = wall`, `T = 1273.15` |
| Faces 2–4 | `type = wall`, `T = 273.15` |
| Faces 5–6 | `null` (2-D plane) |
| Material | k = 52 W/mK, ρ = 7850 kg/m³, c<sub>p</sub> = 460 J/kgK |
| Initial condition | 273.15 K |
| Time scheme | RK3, `vnn = 1.0`, steady |
| Integration variables | `cons` |
| IRS | `irs-beta = 0.5` |
| Stop | `res-threshold = 1e-8` (iteration 3 051) |

## Check

`verify.py` interpolates the cell temperatures bilinearly at point A (0.2 m, 0.2 m), passes below 1 % from 533.65 K, and writes `2D-steady-state-temperature-plate-field.svg`.

## Run

```bash
cd test/steady-state/temperature_Plate
./FUSS.sh solve && python3 verify.py
```

CTest: `TemperaturePlate` (labels fast, validation).
