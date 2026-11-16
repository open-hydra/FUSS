# Multilayer Slab

The physics, the reference solution and the results are on the V&V page:
[Multilayer Slab](../../../docs/vv/1D-steady-state-multilayer-slab.md).

## Setup

| Parameter | Value |
|---|---|
| Mesh | one block, 24 × 2 cells over x ∈ [0, 0.015] m (8 / 6 / 10 cells per layer) |
| Faces 1, 2 | 302 at 1000 K and 400 K |
| Faces 3–6 | adiabatic |
| Initial condition | 700 K |
| Time scheme | RK3, `vnn = 0.25`, steady |
| Integration variables | `cons` |
| Stop | `res-threshold = 1e-12` (iteration 64 200) |

`make_inputs.py` regenerates `MESH/mesh.tec` and `INPUT/`, which ship with the case (ATLAS does not write a material change inside a block).

## Check

`verify.py` passes when the field is within 10⁻⁶ K of the exact solution and the wall fluxes within 10⁻⁶ relative.

## Run

```bash
cd test/steady-state/multilayer_Slab
./FUSS.sh solve && python3 verify.py
```

CTest: `MultilayerSlab` (labels fast, validation).
