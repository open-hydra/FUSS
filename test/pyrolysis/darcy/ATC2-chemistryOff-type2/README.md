# ATC2-chemistryOff under the type-2 gas closure (plan 07 sections 4.4 and 4.5)

Plan 06's `../../ATC2-chemistryOff` (TACOT with zero emissivity, heated face 303 `hconv 300`, `Tref 2800`, 60 s, 100 cells)
with `[FUSS-Pyrolysis] gas-closure = type2`, `[FUSS-Darcy] p-init = p-wall = 101325`, `solver = tridiag`, TACOT's permeability
(`K_v` 1.6e-11, `K_c` 2e-11 m², isotropic). The heated face is open, the back face impermeable. `.slow` (60 s simulated).
Run: `./FUSS.sh -p 1 solve`, then `python3 verify.py [--type1 <type-1 run of ../../ATC2-chemistryOff>]`;
`run_physics_gates.sh -c darcy_atc2_chemistry_off` does both. The oracle is the PATO run in `ref/` (see `ref/README.md`).

## Gate 4.4 (plan 17: tolerances from the references, no stored FUSS value)

One external code solved this variant, so the T tolerance is the PATO-Amaryllis T spread of Ablation Test Case 2.1 (same
material and heating), computed live from `../../ablation/ATC21/ref` (`tools/abl_tools.spread_atc21`); (acc) and (dir) compare
with the type-1 run the same gate makes (`--type1`; they print SKIP without it).

| check | binary b2b68f65 (2026-11-26) | bound |
|---|---|---|
| (T) max `|T − T_PATO|`, 7 depths, 0.1–60 s | 10.51 K at 1 mm, 0.9 s | the ATC 2.1 T spread, 25.35 K (1 mm, 0.3 s) |
| (acc) the same maximum | 10.51 K | the same gate's type-1 maximum, 14.82 K (1 mm, 1.1 s) |
| (dir) `|T − T_PATO|` at that type-1 maximum's point | 10.38 K | < 14.82 K |
| (p) max `|Δp_FUSS − Δp_PATO|` / PATO's peak `Δp` | 0.005 (11.2 Pa at 24 mm, 0.7 s; FUSS peak 2050.4 Pa, PATO 2049.0 Pa) | 0.1966 |
| (sign) `T_type2 − T_type1` where `T_PATO ≥ 700 K` | max −3.115 K | ≤ +0.05 K |

The (p) bound is plan 07 section 4.4's sensitivity of the pressure rise to T, μ and W, `δT/T + δμ/μ + δW/W` = 8.4 % at
δT = 20 K and 850 K (TACOT data; 0.42 %/K), scaled to δT = the ATC 2.1 T spread (10.65 %), plus ATC 2.1's mdot_g spread over
PATO's mdot_g peak (2.343e-3 / 2.599e-2 = 9.02 %) in place of plan 07's 5 % allowance for mdot_pg: 19.66 %, live.

History: 2026-10-13 the gas leaving through the open face started to carry `h_g(p_w, T_wall)`, the enthalpy at the wall
temperature its BC computes (it took the cell temperature before, an O(dx) flux error that made the exact transpiration
benchmark first order, `test/pyrolysis/exact/transpiration`): type 2 vs PATO at 1 mm went 15.79 → 10.52 K.

Budgets: the solid and gas mass balance closes to 7e-15 relative over the run (`M_s + M_g + M_gas_out`); the energy identity
`E(t) − E(0) + Q_out − Q_src` drifts linearly, −6e-7 J at 10 s to −9.7e-6 J at 60 s (5.4e-11 of the heat input), where type 1
plateaus at −1.4e-6 J after 20 s. The drift is consistent with, not proven to be, the stop band of the temperature inversion
(`|f| ≤ 1e-13 |E|` per cell and inversion, four inversions a step: a bound of order 1e-3 J over the run); a future 1e-10 identity gate on a long type-2 run would
sit within a factor 2 of it.

## Solver health, gate 4.5 (plan 01 §7 Q4's measurement; 2026-10-09, -p 1, machine load ~70 on 96 cores)

Binary `463c99a5`: the committed source differs from it only in I/O (the restart column `T_prev`, the `p_w` sentinel) and a
registry description, not in the pressure solver.

The case for exactly 2000 material steps (`iter-threshold = 2000`, `res-threshold = 1e-30`, t = 6.20 s), `tol = 1e-8`:

| solver | sweeps per step | max resid / tol | solver wall time |
|---|---|---|---|
| `tridiag` | 1 on every step | 1.6e-12 / 1e-8 | 2.15 s (1.08 ms/step) |
| `pseudo` (`pseudo-vnn 0.5`, `max-iter 1000000`) | max 45 851, median 35 101 (mean 35 676, min 15 239; 7.1e7 in all) | 1.0e-8 / 1e-8 (it stops at the tolerance) | 3538 s (1.77 s/step) |

Wall-clock ratio `pseudo/tridiag`: **1643**. The plan's estimate for this 100-cell sealed/open line was 5.6e4 (k = 3) to
1.1e5 (k = 6) sweeps per step (§3.3, f = 4); the measured median is 0.63 of the lower end. The sweep count rises from 1.5e4 at
the start to its maximum at step 208 (t = 0.77 s, 4.6e4) and falls slowly after (4.3e4 at step 500, 2.7e4 at step 2000). The two solvers
give the same trajectory to the pseudo tolerance: over the 2000 steps the probes differ by at most 1.4e-3 K and 6.4e-3 Pa.

Tolerance independence: the full 60 s case rerun with `tol = 1e-10` is byte-identical to the `tol = 1e-8` run (all 30 output
files, `darcy.txt` included): on a one-dimensional block `tridiag` solves the line exactly in its one sweep, so the plan's
1e-6 K / 1e-6 `Δp` bounds hold with zero.
