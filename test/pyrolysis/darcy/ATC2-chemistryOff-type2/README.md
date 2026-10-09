# ATC2-chemistryOff under the type-2 gas closure (plan 07 sections 4.4 and 4.5)

Plan 06's `../../ATC2-chemistryOff` (TACOT with zero emissivity, heated face 303 `hconv 300`, `Tref 2800`, 60 s, 100 cells)
with `[FUSS-Pyrolysis] gas-closure = type2`, `[FUSS-Darcy] p-init = p-wall = 101325`, `solver = tridiag`, TACOT's permeability
(`K_v` 1.6e-11, `K_c` 2e-11 m², isotropic). The heated face is open, the back face impermeable. `.slow` (60 s simulated).
Run: `./FUSS.sh -p 1 solve`, then `python3 verify.py [--type1 <type-1 run of ../../ATC2-chemistryOff>] [--record]`;
`run_physics_gates.sh -c darcy_atc2_chemistry_off` does both. The oracle is the PATO run in `ref/` (see `ref/README.md`).

## Gate 4.4 (2026-10-09, `expected.json`)

| check | measured | budget |
|---|---|---|
| (T) max `|T − T_PATO|`, 7 depths, 0.1–60 s | 15.79 K at 1 mm, 2.1 s | B_T 24.75 K (type 1: 22.75 K) |
| (dir) `|T − T_PATO|` at type 1's peak (1 mm, 2.6 s) | 15.65 K | < 22.75 K |
| (p) max `|Δp_FUSS − Δp_PATO|` / PATO's peak `Δp` | 24.2 Pa / 2048 Pa = 0.012 (FUSS peak 2052.6 Pa) | 0.15 |
| (sign) `T_type2 − T_type1` where `T_PATO ≥ 700 K` | −1.25 … −9.76 K (1 mm −7.3…−1.3, 8 mm −9.8…−8.9) | ≤ +0.05 K |

Per depth, type 2 vs PATO (type 1 vs PATO in parentheses): 1 mm 15.79 (22.75), 2 mm 9.72 (18.05), 4 mm 5.50 (14.56),
8 mm 2.77 (12.24), 12 mm 1.69 (11.13), 16 mm 1.27 (9.30), 24 mm 0.42 (2.80) K: removing the type-1 approximation brings every
depth closer to PATO, most of all where the char layer is thick.

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
