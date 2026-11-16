# Volumetric Heat Plate

The physics, the reference solution and the results are on the V&V page:
[Volumetric Heat Plate](../../../docs/vv/2D-steady-state-qvol-plate.md).

## Setup

| Parameter | Value |
|---|---|
| Mesh | 60 × 40 cells, 0.6 m × 0.4 m, uniform |
| Faces 1–4 | `type = wall`, `T = 273.15` |
| Faces 5–6 | `null` (2-D plane) |
| Source | `[STB-Block1]` `qvol = 1e6` W/m³ (ATLAS STB writes `INPUT/st.tec`) |
| Material | k = 52 W/mK, ρ = 7850 kg/m³, c<sub>p</sub> = 460 J/kgK |
| Initial condition | 273.15 K |
| Time scheme | RK3, `vnn = 1.0`, steady |
| Integration variables | `cons` |
| IRS | `irs-beta = 0.5` |
| Stop | `res-threshold = 1e-8` (iteration 4 110) |

## Check

`verify.py` interpolates the cell temperatures at point B (0.3 m, 0.2 m), passes below 1 % from the NAFEMS value, and writes `2D-steady-state-qvol-plate-field.svg`.

## Run

```bash
cd test/steady-state/qvol_Plate
./FUSS.sh solve && python3 verify.py
```

CTest: `QvolPlate` (labels fast, validation).
