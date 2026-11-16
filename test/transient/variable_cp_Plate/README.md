# Variable-cp Plate

The physics, the reference solution and the results are on the V&V page:
[Variable-cp Plate](../../../docs/vv/2D-transient-variable-cp-plate.md).

## Setup

| Parameter | Value |
|---|---|
| Mesh | one block, 60 × 40 cells, 0.6 m × 0.4 m (as `qvol_Plate`) |
| Faces 1–4 | `q = 0` (301) |
| Faces 5–6 | `null` |
| Source | 5 × 10⁷ W/m³, uniform |
| Material | c<sub>p</sub> = 300 + 0.4 T, ρ = 7850 kg/m³, k = 52 W/mK |
| Initial condition | 300 K |
| Time scheme | RK3, `vnn = 0.4`, `time-accurate = true` (54 steps, Δt from 1.27 s to 2.49 s) |
| Integration variables | `prim` (set `cons` to check the other path) |
| End time | `time-threshold = 100` (last step ends at 101.48 s); field every 5 iterations |

`make_properties.py` regenerates `INPUT/properties.dat` with the Energy column h = ρ (A T + B T²/2); ATLAS GPB writes constant-property tables only.

## Check

`verify.py` reads every `OUTPUT/field*.tec`, compares the mean temperature with the closed form and the stored energy with q̇ t, checks that the field is uniform, and writes `2D-transient-variable-cp-plate.svg`. It passes below 10⁻³ of the temperature rise.

## Run

```bash
cd test/transient/variable_cp_Plate
./FUSS.sh solve && python3 verify.py
```

CTest: `VariableCpPlate` (labels fast, validation).
