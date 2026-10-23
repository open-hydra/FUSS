# receding_front: the wall state a receding face sweeps out (plan 12 section D, P-D3)

Plan 08's P-D3: a receding ablating wall (BC 306) sweeps out the WALL state (`T_w`, the solid at the wall, `p_w` and
the stored gas at `p_w`), as PATO's patch values do. Stage C built it (`Remap_Pyrolysis_ALE`, `Wall_Q`); until this
case only the ATC 2.2/2.3 comparisons exercised it. `generate.py` writes four runs; `verify.py` checks them.

| run | what | checks |
|---|---|---|
| `front/` | a 0.1 m line of 200 cells (constant properties `rho_v` 280, `cp` 1000, `k` 0.4; TACOT kinetics frozen, `A = 0`) receding at the TEST ONLY `sdot-override` 1/700 m/s (law `nodeshift`) under an absorbed flux `q = 4e5 W/m2`; 280 steps | (G1) `T_w` against the steady closed form `T_w* = T_inf + q/(rho c sdot) = 1300 K`; (G2) the swept energy of every step at the face's `T_w` of the step before; (G3) the energy budget; (G4) mass and position; (G5) a numpy replica of the scheme; (G6) printed |
| `twin_pyro/` | `front` with a Robin wall, `hconv` 400 to 2300 K (`T_w*` again 1300 K), pyrolysis remap | within 2e-3 of `T_w*` |
| `twin_legacy/` | `twin_pyro` without `[FUSS-Pyrolysis]`: the legacy remap, which sweeps the wall state too (plan 12 U19) | within 2e-3 of `T_w*` |
| `t2/` | `front` under the type-2 closure, `K = 1e-12 m2`, `p_w = p-init = 101325 Pa`, 20 steps with a field file each | (T2-a) swept energy at (`T_w`, `p_w`), (T2-b) swept gas `eta_st p_w`, every step; (bal-m), (bal-e) |

```
python3 generate.py                                   # rewrites front/ twin_pyro/ twin_legacy/ t2/
python3 verify.py RUN/front RUN/twin_pyro RUN/twin_legacy
python3 verify.py --t2 RUN/t2
```

Why the closed form is exact enough to gate: in the surface frame the steady profile is
`T = T_inf + (q/(rho c sdot)) exp(-sdot xi/a)`, so the wall settles where the heat in raises the incoming solid from
`T_inf` to the temperature it leaves at, whatever `k`. Discretely, on a translating mesh, the remap's exact conservation
and the RK's exact flux integral give `q dt = delta rho c (T_swept - T_back)` per step: with the wall state swept,
`T_w = T_back + q/(rho c sdot)` to round-off at any `dx` and `dt`. FUSS's `nodeshift` compresses the line instead, so
the steady state is approached, not held: at `t* = sdot^2 t/a = 14.9` the transient left is about +-1e-3 of the 1000 K
rise, against +0.20 for a remap that sweeps the cell state.

Measured 2026-10-17 (FUSS `d6b6e91`, gnu, 1 thread, 3 s for `front`; to be re-measured on the branch tip it lands on):

| check | value | bound |
|---|---|---|
| (G1) | `T_w` 1300.013189 K, +1.32e-5 | 2e-3 |
| (G2) | 2.7e-13 | 1e-12 |
| (G3) | 1.4e-15 | 1e-12 |
| (G4) | 5.4e-14; `s_total` 2.3e-16; nodes 2.8e-17 m | 1e-13; 1e-13; 1e-15 m |
| (G5) | `T_w` 4.3e-14, `T` 1.3e-13 | 1e-9 |
| twin_pyro | -3.80e-4 | 2e-3 |
| twin_legacy | -3.80e-4 (2026-10-22, plan 12 U19; +9.04e-2 while the legacy remap swept the cell state); the twins' last fields agree to 1.0e-12 K | 2e-3 |
| (T2-a), (T2-b) | 1.8e-13, 1.8e-13 | 1e-12 |
| (bal-m), (bal-e) | 4.9e-15, 2.2e-15 | 1e-12 |
| `p_c - p_w` at step 20 | 20.5 Pa (2.0e-4 of `p_w`) | printed |

The per-step identities stop near 3e-13 because the swept volume is a difference of node positions near 0.1 m for a
step of 6e-5 m. With the flux raised by 1 % (the gate's selftest) exactly (G1) and (G5) fail.

Facts that shaped the case (plan 12 section D.7): step 1 of a new run sweeps the surface CELL's T (`Tw_prev` is 0 until
the first BC call), so (G2) uses `ic.tec` for it; under type 2 the swept stored gas is the cell's `eta_st` times `p_w`
(not `eta(T_w, p_w) p_w`); a receding face without `[FUSS-Pyrolysis]` swept the cell state (`twin_legacy`, +9.04e-2)
until plan 12 U19 gave the legacy remap the wall state (2026-10-22).
