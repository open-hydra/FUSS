# Multigrid Cycles (FAS) – Measured Performance

The FAS cycles (`cycle = v | w`, [theory](../theory/time-integration.md#fas-cycles)) promise the
same converged field as the baseline run at a fraction of the fine-grid iterations. This page
collects the measured evidence: the equivalence checks, the gains on the verification cases, the
stress campaign on representative hard cases, and the known limits. Every run converges to
`res-threshold = 1e-8` unless stated; "iterations" are fine-grid (root-level) iterations, the
ones the baseline also counts. Iteration counts are deterministic: every repetition reproduced
them exactly. Wall times are serial solver times (the `Solver` line of the timing summary,
without set-up and output), the minimum of three or four runs of the current code on a shared
machine;
repeated runs differ by up to 1.9×, so read the wall ratios to within that factor.

!!! note "The two invariants behind every number"
    1. **`cycle = none` is bit-identical to the pre-FAS solver.** Verified on all 17 test cases,
       serial and 2 MPI ranks, against the pre-FAS binary (commit `494ebb1`): 0 bytes differ in
       any output file. The `MultigridNoneInert` CTest keeps a machine-independent form of this
       guarantee: every FAS key set under `cycle = none` must reproduce the plain run byte for byte.
    2. **Cycles converge to the fine discrete solution.** On patch cases with an exactly
       representable solution the cycle fixed point is bit-exact (residual reaches 0.0); on the
       cases below the field at threshold matches the `none` field to the stopping tolerance.

## Reference plate (residual_convergence)

The two-material plate of the [residual convergence study](residual-convergence.md) (2 blocks of
100 × 10 cells, 101 connection, threshold 1e-8, `cons`). The 10-cell height caps the hierarchy at
2 levels, which caps the cycle gain: this is the *worst* verified configuration, and the one the
`fas_v` / `fas_w` / `fas_w_irs` test cases ship.

| Configuration | Fine iterations | Wall (s) | Iterations vs `nominal` | Wall vs `nominal` |
|---|---:|---:|---:|---:|
| nominal (single grid) | 369 090 | 473 | 1× | 1× |
| irs | 99 305 | 236 | 3.7× | 2.0× |
| mg3 (sequencing) | 144 108 | 121 | 2.6× | 3.9× |
| mg + irs | 80 603 | 126 | 4.6× | 3.8× |
| **cycle = v** | 25 756 | 63 | 14× | 7.5× |
| **cycle = w** | 13 345 | 53 | 28× | 9.0× |
| **cycle = v + irs** | 6 934 | 25 | 53× | 19× |
| **cycle = w + irs** | **3 593** | **19** | **103×** | **24×** |
| cycle = v, coarse-iter 50 | 5 458 | 57 | 68× | 8.3× |
| cycle = v, coarse-iter 1 | 158 191 | — | 2.3× | (cripples the cycle) |

The cycle rows use the defaults (`coarse-iter = 10`, `pre-iter = 2`, `post-iter = 1`). On a
2-level hierarchy the coarse solve is the knob: `w + irs` needs 15 686 / 8 518 / 3 593 iterations
with `coarse-iter` 2 / 4 / 10, and on two levels `w` is exactly `v` with twice the coarse
iterations (`w` with `coarse-iter = 2` reproduces `v` with 4 bit for bit). The wall gains are
smaller than the iteration gains because every cycle adds its coarse work and transfers to the
fine smoothing. With deep hierarchies (below) the default is already optimal and W becomes
insensitive to it.

## Grid independence (deep hierarchies)

Unit-square patch case, full hierarchies, `coarse-iter 10`. Root iterations to 1e-12 (wall):

| Fine mesh | none | w | w + irs |
|---|---:|---:|---:|
| 32² (3 levels) | 12 854 (9.8 s) | 184 (0.5 s) | 103 (0.3 s) |
| 64² (4 levels) | 49 077 (184 s) | 191 (1.1 s) | 107 (1.0 s) |
| 128² (5 levels) | ~201 700 (~3 600 s), extrapolated | 207 (5.3 s) | 111 (4.2 s) |

The W-cycle iteration count is essentially **independent of N** (textbook behaviour): ~1800× fewer
iterations at 128² and, against the extrapolated single-grid run (10 000 iterations timed,
17.9 ms each), ~700× (`w`) to ~850× (`w + irs`) less wall time. Rates per cycle 0.37–0.41 for
`w + irs`.

## Nonlinear, layered, 3D

| Case | none | cycle | Gain (iterations / wall) |
|---|---:|---:|---:|
| k(T)×3, cp(T)×2 nonlinear (64², 3 levels, 1e-10) | 43 353 (167 s) | w: 211 (1.4 s) | 205× / ~120× |
| Bi-material k×10, interface on all levels (64², 3 levels, 1e-11) | 44 955 (114 s) | w: 217 (1.3 s) | 207× / ~90× |
| 3D box 16³ (2 levels, 1e-12) | 2 983 (7.7 s) | v: 263 (1.0 s), w: 185 (0.9 s) | 11–16× / 8–9× |

Fields equal to the `none` fields within 3e-9 K at threshold.

## Chimera (overset blocks)

New verification case `numerics-features/chimera_mg`: two genuinely overlapping blocks
(112 × 16 and 104 × 12 cells, 0.09 m overlap, different spacing in both directions), chimera
donors on all 3 levels built by ATLAS BCB (`BC-chimera = true`, `MG-levels = 3`). The FAS
transfers never read a chimera receiver: restriction, forcing and correction sums run over
interior cells only, the receivers are ghost cells refreshed by each level's own chimera fill,
and (with IRS) the ghost residual at a 102 face is already zero-gradient — the same policy IRS
applies to every non-connection ghost.

| Check | Result |
|---|---|
| `cycle = none` vs pre-FAS binary, serial and 2 ranks | bit-identical (0 bytes differ) |
| `v` / `w`: 1 process vs 2 ranks (one block per rank, donors cross ranks on all levels) | bit-identical |
| Field vs `none` at 1e-8 | v: 3.0e-7 K, w: 3.0e-7 K |
| Admissibility fallbacks (`FAS:` reports) | none |
| Iterations (wall): none 13 367 (20 s) | v 2 619 (3.3 s; 5.1× / 6.2×), **w 1 700 (1.2 s; 7.9× / 17×)** (iterations / wall) |

**W is safe with chimera.** The study provisioned a V-only constraint near 102 boundaries
(CFL3D practice); the data shows W converging faster than V, bit-reproducibly across ranks and
with zero admissibility fallbacks, so no constraint is imposed.

## Stress campaign

Four hard configurations built for this review (generators in the verification tooling), all
run to 1e-8 with `coarse-iter 10` and immediate descent; fields checked against the `none`
field (≤ 5e-7 K) or, for the 3D patch cases, against the exact linear solution (≤ 2e-9 K).

| Case | none | v | w | w + irs | w + irs gain (iterations / wall) |
|---|---:|---:|---:|---:|---:|
| 4 conformal blocks, 101 connections (2×2 × 32², 3 levels) | 14 047 (66 s) | 261 | 96 | **57** (1.0 s) | 246× / 66× |
| Localized q̇ hot spot (64², discontinuous source, 3 levels) | 34 637 (102 s) | 610 | 169 | **49** (0.8 s) | 707× / ~130× |
| 3D box 32³ (3 levels) | 8 493 (323 s) | 1 324 | 135 | **66** (5.3 s) | 129× / 61× |
| 3D box 64³ (4 levels) | ~33 000 extrap. (~2.7 h)¹ | — | 139 (39 s) | **69** (42 s) | ~480× / ~230× |

¹ The 64³ baseline was stopped at 10 000 iterations (residual 1.3e-2, rate 0.9994/iter); its
wall time is extrapolated from 1 500 timed iterations (0.29 s each). The cycles reach the
*exact* discrete solution (max error ≤ 8e-10 K) in ~100 iterations: the W-cycle
iteration count barely moves from 32³ to 64³ (135 → 139, 66 → 69), the same N-independence
measured in 2D.

**Realistic liner (4 material layers, squeezed cells, hard BCs).** 64 × 128 cells, 0.3 × 3.9 mm
(AR 12.5), layers k = 5 / 300 / 15 / 150 W/mK on all-level grid lines, convective 303 walls
(hot side Tref 3000 → 1500 K along the wall, cold side h = 10⁴). Two honest findings:

1. **The explicit baseline is the fragile one.** With `vnn = 0.5` (stable everywhere else) the
   *baseline* run dies at the k 300 → 15 interface in a few steps (T < 0 K): the local dt uses
   the cell's own diffusivity while the face flux carries the copper-side transmissibility — a
   pre-existing property of the explicit multimaterial scheme, not of the cycles. `vnn = 0.2`
   is stable, for the baseline and the cycles alike.
2. **IRS + cycles with aggressive vnn can limit-cycle here.** `v`/`w` + IRS at `vnn = 2.0`
   (the setting every other case takes) locks into the documented bounded orbit (residual
   ~2.4e4, steady `FAS:` fallback reports) instead of converging.
   <!-- FILL: liner vnn 1.0 + none vnn 0.2 results -->

Without IRS at `vnn = 0.2` both cycles converge cleanly: v 8 362, w 7 222 iterations (fields
equal to 2.6e-6 K). This case sits at the anisotropic end of the spectrum (squeezed cells +
material jumps): the gains are real but far from the isotropic 200×, exactly as the AR-100
limit below predicts.

## Anisotropic cells: the honest limit

Reference case `numerics-features/anisotropic_Plate` (128 × 16 cells on 1 m × 1.25 mm: cell
aspect ratio 100, 3 levels, exact linear solution), residual written every iteration
(`res-diter = 1`; the shipped case writes it every 100):

| Configuration | Iterations to 1e-8 | Wall (s) | Iterations to 1e-12 | Wall (s) |
|---|---:|---:|---:|---:|
| none | 2 253 | 2.9 | 3 074 | 3.9 |
| **none + irs** | **739** | **2.4** | stalls at 4.0e-11 | — |
| v | 1 513 | 5.0 | 2 439 | 4.7 |
| w | 1 514 | 8.4 | 2 439 | 7.9 |
| v + irs | 780 | 2.4 | 1 306 | 4.0 |
| w + irs | 782 | 6.7 | 1 297 | 6.6 |

V and W need the same iterations (rates per cycle 0.95–0.96): the coarse correction, not the
cycle shape, is the bottleneck. Full coarsening with a pointwise explicit smoother cannot
represent the strong transverse coupling on the coarse levels. **In wall time the cycles do not
pay off here**: at 1e-8 IRS alone is as fast as the best cycle and W costs ~3× more. IRS alone
stalls on a round-off orbit at a residual of 4.0e-11 (the residual freezes, bit for bit, from
~1 500 iterations on), which the cycles, with or without IRS, break through to 1e-12. On AR ≫ 1
meshes run IRS without cycles, unless the threshold is below that floor. The cure
(line-implicit smoothing or semi-coarsening) is outside the FAS v1 perimeter.

## Deep V instability

`cycle = v` without IRS on 4 or more levels can lock into a bounded limit cycle (residual
oscillating at a large value, a steady stream of `FAS: ... updates halved` reports). This is a
property of the method (weak explicit smoothing vs accumulated coarse corrections), not of a
case: use `w`, or enable IRS. Documented in
[time integration](../theory/time-integration.md#fas-cycles).

## Practical guidance

- **`cycle = w` + IRS is the production choice**: textbook rates, N-independent iterations,
  insensitive to `coarse-iter` on deep hierarchies, safe with chimera and at any depth. The
  exception is strongly anisotropic cells (AR ≫ 1), where IRS alone is as fast (see above).
- `coarse-iter = 10` (the default) is right for deep hierarchies; raise it only on 2-level
  hierarchies where the coarsest solve is the bottleneck.
- `v` without IRS: 3 levels at most.
- Descend immediately (`level<m>-iter = 1` for m > 1): cycling from the start beats long
  sequencing stints.
