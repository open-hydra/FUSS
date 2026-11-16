# Multi-Material Plate

The physics, the reference solution and the results are on the V&V page:
[Multi-Material Plate](../../../docs/vv/2D-transient-multimat-plate.md).

## Setup

| Parameter | Value |
|---|---|
| Mesh | three blocks, 15 × 40, 10 × 40 and 15 × 40 cells (Δx = Δy = 0.05 m), joined by `connection` (101) |
| Faces 3, 4 | `T = 101` (y = 0) and `T = 100` (y = 2 m) |
| Outer x faces | `q = 0` |
| Materials | insulator (k = 10⁻¹²) in blocks 1 and 3; conductor (k = 1) in block 2; ρ = c<sub>p</sub> = 1 |
| Initial condition | 100 K |
| Time scheme | RK3, `vnn = 0.4`, `time-accurate = true` (Δt = 5 × 10⁻⁴ s) |
| Integration variables | `prim` |
| End time | 0.5055 s; field every 10 iterations (`sol-overwrite = false`) |

## Check

`verify.py` reads column 6 of block 2 (x = 0.025 m) in `field1`, `field21`, …, `field101`, normalises it, and compares it with the 100-term Fourier series. It passes when the largest L∞ error is below 0.02.

## Run

```bash
cd test/transient/multimat_Plate
./FUSS.sh solve && python3 verify.py
```

CTest: `MultimatPlate` (labels fast, validation).
