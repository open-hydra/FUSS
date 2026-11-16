# Convective Plate

The physics, the reference solution and the results are on the V&V page:
[Convective Plate](../../../docs/vv/2D-steady-state-convective-plate.md).

## Setup

| Parameter | Value |
|---|---|
| Mesh | 60 × 100 cells, 0.6 m × 1.0 m, uniform (Δx = Δy = 0.01 m) |
| Face 1 (x = 0) | `type = wall`, `q = 0` |
| Faces 2, 4 | `type = wall`, `hconv = 750`, `Tref = 273.15`, `qrad = 0` (303) |
| Face 3 (y = 0) | `type = wall`, `T = 373.15` |
| Faces 5–6 | `null` (2-D plane) |
| Material | k = 52 W/mK, ρ = 7850 kg/m³, c<sub>p</sub> = 460 J/kgK |
| Initial condition | 373.15 K |
| Time scheme | RK3, `vnn = 1.0`, steady |
| Integration variables | `prim` |
| IRS | `irs-beta = 0.5` |
| Stop | `res-threshold = 1e-8` (iteration 21 838) |

## Check

`verify.py` averages the two wall cells around point A (0.6 m, 0.2 m) and passes below 1 % from 18.25 °C.

## Run

```bash
cd test/steady-state/convective_Plate
./FUSS.sh solve && python3 verify.py
```

CTest: `ConvectivePlate` (labels validation (about two minutes, not in `fast`)).
