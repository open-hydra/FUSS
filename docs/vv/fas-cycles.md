# Multigrid Cycles (FAS) – Measured Performance

The FAS cycles (`cycle = v | w`, [theory](../theory/time-integration.md#fas-cycles)) promise the
same converged field as the baseline run at a fraction of the fine-grid iterations. This page
collects the measured evidence: the equivalence checks, the gains on the verification cases, the
stress campaign on representative hard cases, and the known limits. Every run converges to
`res-threshold = 1e-8` unless stated; "iterations" are fine-grid (root-level) iterations, the
ones the baseline also counts.

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

| Configuration | Fine iterations | Wall (s) | vs `nominal` |
|---|---:|---:|---:|
| nominal (single grid) | 369 090 | 2 213 | 1× |
| irs | 99 305 | 1 051 | 3.7× |
| mg3 (sequencing) | 144 108 | 806 | 2.6× |
| mg + irs | 80 603 | 388 | 4.6× |
| **cycle = v** | 58 285 | 523 | 6.3× |
| **cycle = w** | 31 642 | 342 | 11.7× |
| **cycle = v + irs** | 15 686 | 225 | 23.5× |
| **cycle = w + irs** | **8 518** | **152** | **43×** |
| cycle = v, coarse-iter 50 | 5 458 | 292 | 68× |
| cycle = v, coarse-iter 1 | 158 191 | — | (cripples the cycle) |

On a 2-level hierarchy the coarse solve is the knob: `coarse-iter 50` buys iterations at extra
coarse work. With deep hierarchies (below) the default `coarse-iter = 10` is already optimal and
W becomes insensitive to it.

## Grid independence (deep hierarchies)

Unit-square patch case, full hierarchies, `coarse-iter 10`. Root iterations to 1e-8:

| Fine mesh | none | w | w + irs |
|---|---:|---:|---:|
| 32² (3 levels) | 12 854 | 184 | 103 |
| 64² (4 levels) | 49 077 | 191 | 107 |
| 128² (5 levels) | ~201 700 (extrapolated) | 207 | 111 |

The W-cycle iteration count is essentially **independent of N** (textbook behaviour): ~1800× fewer
iterations at 128², ≥200× in wall time. Rates per cycle 0.37–0.41 for `w + irs`.

## Nonlinear, layered, 3D

| Case | none | cycle | Gain |
|---|---:|---:|---:|
| k(T)×3, cp(T)×2 nonlinear (64²) | 43 353 (708 s) | w: 211 (11 s) | 205× |
| Bi-material k×10, interface on all levels (64²) | 44 955 | w: 217 | 207× |
| 3D box 16³ (3 levels) | 2 983 | v: 263, w: 185 | 11–16× |

Fields equal to the `none` fields within 3e-9 K at threshold.

## Chimera (overset blocks)

New verification case `numerics-features/chimera_mg`: two genuinely overlapping blocks
(112 × 16 and 104 × 12 cells, 0.09 m overlap, different spacing in both directions), chimera
donors on all 3 levels built by ATLAS BCB (`BC-chimera = true`, `MG-levels = 3`). The FAS
transfers never read a chimera receiver: restriction, forcing and correction sums run over
interior cells only, the receivers are ghost cells refreshed by each level's own chimera fill,
and (with IRS) the ghost residual at a 102 face is already zero-gradient — the same policy IRS
applies to every non-connection ghost.

<!-- FILL: chimera numbers -->

| Check | Result |
|---|---|
| `cycle = none` vs pre-FAS binary, serial and 2 ranks | bit-identical (0 bytes differ) |
| `v` / `w`: 1 process vs 2 ranks (one block per rank, donors cross ranks on all levels) | bit-identical |
| Field vs `none` at 1e-8 | v: 3.0e-7 K, w: 3.0e-7 K |
| Admissibility fallbacks (`FAS:` reports) | none |
| Iterations (wall): none 13 367 (70 s) | v 2 619 (13 s, 5.1×), **w 1 700 (5 s, 7.9×)** |

**W is safe with chimera.** The study provisioned a V-only constraint near 102 boundaries
(CFL3D practice); the data shows W converging faster than V, bit-reproducibly across ranks and
with zero admissibility fallbacks, so no constraint is imposed.

## Stress campaign

<!-- FILL: stress numbers -->

## Anisotropic cells: the honest limit

Reference case `numerics-features/anisotropic_Plate` (128 × 16 cells on 1 m × 1.25 mm: cell
aspect ratio 100, 3 levels, exact linear solution).

<!-- FILL: anisotropic numbers -->

Full coarsening with a pointwise explicit smoother cannot represent the strong transverse
coupling on the coarse levels: the cycles converge but gain little. The cure (line-implicit
smoothing or semi-coarsening) is outside the FAS v1 perimeter; use IRS, which roughly doubles
the gain, and expect single-digit speed-ups on AR ≫ 1 meshes.

## Deep V instability

`cycle = v` without IRS on 4 or more levels can lock into a bounded limit cycle (residual
oscillating at a large value, a steady stream of `FAS: ... updates halved` reports). This is a
property of the method (weak explicit smoothing vs accumulated coarse corrections), not of a
case: use `w`, or enable IRS. Documented in
[time integration](../theory/time-integration.md#fas-cycles).

## Practical guidance

- **`cycle = w` + IRS is the production choice**: textbook rates, N-independent iterations,
  insensitive to `coarse-iter`, safe with chimera and at any depth.
- `coarse-iter = 10` (the default) is right for deep hierarchies; raise it only on 2-level
  hierarchies where the coarsest solve is the bottleneck.
- `v` without IRS: 3 levels at most.
- Descend immediately (`level<m>-iter = 1` for m > 1): cycling from the start beats long
  sequencing stints.
