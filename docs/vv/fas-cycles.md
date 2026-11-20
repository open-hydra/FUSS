# Multigrid Cycles (FAS) – Measured Performance

The FAS cycles (`cycle = v | w`, [theory](../theory/time-integration.md#fas-cycles)) promise the
same converged field as the baseline run at a fraction of the fine-grid iterations. This page
collects the measured evidence: the six configurations none / IRS / V / W / V + IRS / W + IRS on
four cases of `test/` (two figures), the equivalence checks, the gains on the verification cases,
the stress campaign on representative hard cases, and the known limits. "Iterations" are
finest-level (root-level) iterations as the residual history counts them, the ones the baseline
also counts; with cycles they advance by `pre-iter + post-iter` per cycle and include the
one-iteration descent on each coarse level. Iteration counts are deterministic: every
repetition reproduced them exactly. Wall times come from a shared machine and are only as
reproducible as it is idle: repeated runs differ by up to 1.9×, so read the wall ratios to
within that factor. Each section states its threshold and which timing line it reports.

The iteration counts and wall times (‡) of this page were re-measured on 2026-10-29 with the binary of the post-merge
branch (`post-merge/sync-followups`, ifx 2024.1 with IPO, one thread per run, `vnn` in the current convention), after
`res-threshold` became relative to the first residual of the finest level: every count is that of a run that ends when its
residual has fallen by its threshold from that first residual, which ends most runs earlier than the absolute threshold of
the release before (the counts of that release are quoted where they help). The exceptions are the liner rows that describe
the code before the fix and the residual-floor column of the anisotropic plate. The binary carries the series face
conductivity at material jumps inside a block and the face-conductivity time step (both always on), the stencil face metric
and the exact cell volume. Wall times are valid only on an idle node; the node of the re-measurement was loaded (120 to 150
of 192 cores), so read the wall ratios as those of one loaded node.
Blocks of one material are bit-identical under the two
face-conductivity changes, and at a block connection between two materials (the reference plate) only
the time step of the cells on the less conductive side shortens; a material jump inside a block (the
bi-material patch of the nonlinear table, the liner) now carries the series face flux instead, and
the liner, the case the defect had broken, is the one described in detail below.

!!! note "The two invariants behind every number"
    1. **`cycle = none` is the solver without FAS.** When the cycles were introduced this was
       verified on all 17 test cases, serial and 2 MPI ranks, against the pre-FAS binary (commit
       `494ebb1`): 0 bytes differ in any output file. The face-conductivity time step of October
       2026 has since changed, by design, the 7 cases with a block connection between two
       materials (`nominal`, `irs`, `chimera`, `mg1`–`mg3`, `mg_irs`); the 9 one-material cases
       were byte-identical across that change, and `time_convergence/simdir` differs only in its
       `diagnostic.tec`, by round-off in cells still at the initial condition. The `MultigridNoneInert` CTest keeps the guarantee in a
       machine-independent form: every FAS key set under `cycle = none` must reproduce the plain
       run byte for byte.
    2. **Cycles converge to the fine discrete solution.** On patch cases with an exactly
       representable solution the cycle fixed point is bit-exact (residual reaches 0.0); on the
       cases below the field at threshold matches the `none` field to the stopping tolerance.

## Six configurations on four cases

The runner and the plot script live in `test/numerics-features/convergence_figure/`
(`run_configs.sh`, `plot_convergence.py`; the runs take about an hour, two at a time). Four
cases of `test/numerics-features/`, every run RK3, `cons`, `res-threshold = 1e-10`, the residual
written every iteration, one thread:

| Case | Mesh | Levels |
|---|---|---|
| plate (`residual_convergence/fas_v`) | two materials, 2 blocks of 100 × 10 cells joined by a 101 connection | 2 |
| liner (`multimat_liner`) | four-layer liner wall, 64 × 128 cells of 0.31 × 3.9 mm in 2 blocks (101), material jumps inside the blocks (k 5 / 300 / 15 / 150 W/mK), convective walls | 4 |
| chimera (`chimera_mg`) | two overlapping blocks (112 × 16 and 104 × 12 cells, chimera donors on every level) | 3 |
| anisotropic (`anisotropic_Plate`) | 128 × 16 cells of aspect ratio 100, one material | 3 |

`none` and `IRS` run on one grid; the cycles on the case's levels, descending with one
iteration per coarse level and the default `pre-iter = 2`, `post-iter = 1`, `coarse-iter = 10`.
Without IRS `vnn = 0.5`; with IRS `vnn = 1.0` and `irs-beta = 0.5`. All `vnn` values of this page are in the
current convention. The iteration counts and the wall times (‡) of this page were re-measured on 2026-10-29, after `res-threshold` became relative to the first residual of the finest level (see [Time Integration](../theory/time-integration.md)): every run ends when its residual has fallen by `res-threshold` from that first residual. The wall times are the `Elapsed` or `Solver` line of one run, one core per run, several runs at a time on a node loaded at 120 to 150 of its 192 cores: they compare the runs of one table, not machines.

![Residual against finest-level iterations, six configurations on four cases](images/fas-cycles-convergence-iterations.png)

*Residual (‖ΔT‖ per iteration) against finest-level iterations. Colour is the cycle (blue
`none`, orange `v`, green `w`), a dashed line is the same configuration with IRS. The top row shows
the full range of each case, the bottom row zooms on the range of the four FAS cycles (the
anisotropic plate, where the cycles end with the single-grid runs, keeps its full range). The
numbers at the top right are the iterations at which each run ends, when its residual has fallen
by 1e-10 from the first residual of its finest level.*

![Residual against wall time, six configurations on four cases](images/fas-cycles-convergence-time.png)

*The same histories against wall time: the `Elapsed` line of each run's timing summary
(set-up, iterations and output), spread uniformly over the iterations of the run. One core per
run, several runs at a time on a node loaded at 120 to 150 of its 192 cores (‡): the times carry
that noise, the iteration counts do not.*

Iterations to a fall of the residual by 1e-10 (wall time, s ‡):

| Configuration | plate | liner | chimera | anisotropic |
|---|---:|---:|---:|---:|
| none | 166 578 (184) | 244 325 (1 266) | 9 209 (12.8) | 1 022 (1.0) |
| IRS | 84 215 (142) | 123 712 (995) | 4 679 (11.9) | 565 (1.0) |
| V | 12 264 (27) | 2 020 (17.8) | 178 (0.6) | 749 (1.3) |
| W | 6 355 (20) | 2 061 (28.7) | 77 (0.6) | 749 (2.3) |
| V + IRS | 6 656 (21) | 2 704 (34.2) | 99 (0.6) | 680 (2.1) |
| **W + IRS** | **3 449 (17)** | **2 521 (49.5)** | **71 (0.6)** | **682 (3.1)** |
| W + IRS vs none, iterations / wall ‡ | 48× / 11× | 97× / 26× | 130× / 21× | 1.5× / 0.3× |

‡ Wall times and wall ratios: the `Elapsed` line, 2026-10-29, on a loaded node (see above).

The iteration counts and both figures are from the 24 runs of
`test/numerics-features/convergence_figure/run_configs.sh` of 2026-10-29 (`res-threshold = 1e-10`, relative; none of
the runs prints a `FAS:` or `IRS:` report; branch `post-merge/sync-followups`, ifx 2024.1 with IPO). With the absolute
threshold of the release before, the same runs took 238 496 / 347 456 / 12 841 / 1 412 iterations (`none`) and
4 450 / 3 615 / 94 / 1 030 (`w + irs`) on the four cases. The FAS acceptance set (`fuss-stress/results/fas-vv-post-i1b/fas_vv.md`)
reruns some of the same configurations at 1e-8 (for instance liner `none` 173 843, `w` 1 370, `w + irs` 3 211 iterations;
chimera `none` 9 153, `v` 1 914, `w` 1 241), at another threshold, so its counts do not replace these.

What the figures say:

- **Every cycle is one to two orders of magnitude below the single grid** in iterations (14× to
  130×) and 7× to 71× in time on the three isotropic cases; IRS alone halves the single-grid count.
- **W + IRS needs the fewest iterations on the plate and the chimera case**, where the order
  none > IRS > V > W holds and IRS helps each cycle. **On the liner the four cycles end within
  684 iterations of each other** (2 020 to 2 704) and
  the plain `v` is the one with the fewest: neither W nor IRS buys a lower count on this
  four-level, two-block hierarchy.
- **Wall time follows the iterations only where the cycle is cheap relative to the fine
  smoothing.** A W cycle visits its coarse sub-tree twice and IRS adds two Jacobi sweeps per
  stage, so the cost per fine iteration grows from `v` to `w + IRS`: on the plate (2 levels, a
  cheap coarse grid) and the chimera case the time still falls with the iterations, on the
  liner `v` is the fastest on the clock (17.8 s ‡ against 49.5 s ‡ for
  `w + IRS`, both against 1 266 s ‡ for the single grid).
- **On the anisotropic plate nothing pays**: all six runs end within 565–1 022 iterations, IRS
  alone the fewest, and the single-grid runs are the cheapest in time ‡ ([below](#anisotropic-cells-the-honest-limit)).

## Reference plate (residual_convergence)

The two-material plate of the [residual convergence study](residual-convergence.md) (2 blocks of
100 × 10 cells, 101 connection, `cons`, the case's `vnn = 0.25`), to `res-threshold` 1e-8; wall times
(‡) are the `Solver` line (without set-up and output) of one run. The
10-cell height caps the hierarchy at 2 levels, which caps the cycle gain: this is the *worst*
verified configuration, and the one the `fas_v` / `fas_w` / `fas_w_irs` test cases ship.

| Configuration | Fine iterations | Wall (s) | Iterations vs `nominal` | Wall vs `nominal` |
|---|---:|---:|---:|---:|
| nominal (single grid) | 236 927 | 226 ‡ | 1× | 1× ‡ |
| irs | 60 855 | 96 ‡ | 3.9× | 2.4× ‡ |
| mg3 (sequencing) | 163 642 | 110 ‡ | 1.4× | 2.1× ‡ |
| mg + irs | 55 523 | 68 ‡ | 4.3× | 3.3× ‡ |
| **cycle = v** | 17 155 | 36 ‡ | 14× | 6.3× ‡ |
| **cycle = w** | 8 890 | 27 ‡ | 27× | 8.3× ‡ |
| **cycle = v + irs** | 5 023 | 16 ‡ | 47× | 14× ‡ |
| **cycle = w + irs** | **2 605** | **13** ‡ | **91×** | **18×** ‡ |
| cycle = v, coarse-iter 50 | 3 637 | 21 ‡ | 65× | 11× ‡ |
| cycle = v, coarse-iter 1 | 105 296 | 124 ‡ | 2.3× | 1.8× ‡ (cripples the cycle) |

‡ Wall times and wall ratios: 2026-10-29, on a loaded node (see the top of the page). The iteration counts and times: the rows
`rc-*` of the FAS acceptance set (`fuss-stress/results/fas-vv-post-i1b/fas_vv.md`) and, for `v + irs` and the two `coarse-iter`
rows, runs of the same decks with those keys (`fuss-stress/results/p2vv-post-i1b/extras`). Every run ends at 1e-8 of the first
residual of its finest level: the cycles after a one-iteration descent (451 K, against 688 K of `nominal`), the sequencing rows
`mg3` and `mg + irs` after a long coarse stint, from a first level-1 residual about 60 times smaller, so they end about two
orders deeper ([Residual Convergence](residual-convergence.md)). With the absolute threshold of the release before: 369 498,
99 415, 212 414, 79 816, 25 817, 13 378, 6 950, 3 602, 5 471 and 158 481 iterations.

The cycle rows use the defaults (`coarse-iter = 10`, `pre-iter = 2`, `post-iter = 1`). On a
2-level hierarchy the coarse solve is the knob: `w + irs` needs 11 365 / 6 172 / 2 605 iterations
with `coarse-iter` 2 / 4 / 10, and on two levels `w` is exactly `v` with twice the coarse
iterations (`w` with `coarse-iter = 2` reproduces `v` with 4 bit for bit: residual history and field identical,
38 814 iterations). The wall gains are
smaller than the iteration gains because every cycle adds its coarse work and transfers to the
fine smoothing. With deep hierarchies (below) the default is already optimal and W becomes
insensitive to it.

## Grid independence (deep hierarchies)

Unit-square patch case, full hierarchies, `coarse-iter 10`. Root iterations to a fall of the residual by 1e-12 (wall ‡):

| Fine mesh | none | w | w + irs |
|---|---:|---:|---:|
| 32² (3 levels) | 10 627 (5.4 s ‡) | 145 (0.19 s ‡) | 81 (0.17 s ‡) |
| 64² (4 levels) | 41 044 (85 s ‡) | 148 (0.71 s ‡) | 77 (0.57 s ‡) |
| 128² (5 levels) | 158 048 (1 395 s ‡) | 152 (2.8 s ‡) | 77 (2.3 s ‡) |

The 128² single grid runs to the threshold (`fuss-stress/results/p2vv-post-i1b/extras`; the FAS acceptance set caps it at 40 000).

The W-cycle iteration count is essentially **independent of N** (textbook behaviour): ~1 000× (`w`) to
~2 000× (`w + irs`) fewer iterations than the single grid at 128² and ~500× (`w`) to ~600× (`w + irs`) less
wall time ‡. Rates per cycle 0.34–0.36 for `w + irs`.

## Nonlinear, layered, 3D

| Case | none | cycle | Gain (iterations / wall) |
|---|---:|---:|---:|
| k(T)×3, cp(T)×2 nonlinear (64², 3 levels, 1e-10) | 32 236 (78 s ‡) | w: 160 (0.90 s ‡) | 201× / 87× ‡ |
| Bi-material k×10, interface on all levels (64², 3 levels, 1e-11) | 35 118 (83 s ‡) | w: 172 (1.07 s ‡) | 204× / 78× ‡ |
| 3D box 16³ (2 levels, 1e-12) | 2 441 (5.5 s ‡) | v: 200 (0.74 s ‡), w: 143 (0.69 s ‡) | 12–17× / 7.5–8.0× ‡ |

The thresholds are falls of the residual from its first value. Fields equal to the `none` fields within 1.3e-6 K at threshold
(1.3e-6, 7.9e-8 and 3.1e-9 K; with the absolute thresholds of the release before, which ended these runs deeper, 2.8e-9,
2.7e-10 and 2.2e-11 K). The threshold of the 16³ box, 1e-12 of its first residual (515 K), lies above the round-off floor of
the field, so its counts no longer depend on the last bits (with the absolute 1e-12 they did: 267 and 195 iterations under `v`
and `w`, 261 and 186 before the sync, see the [re-validation](#re-validation-at-the-sync-of-october-2026)).

## Chimera (overset blocks)

Verification case `numerics-features/chimera_mg`: two genuinely overlapping blocks
(112 × 16 and 104 × 12 cells, 0.09 m overlap, different spacing in both directions), chimera
donors on all 3 levels built by ATLAS BCB (`BC-chimera = true`, `MG-levels = 3`; the files
of ATLAS main `7d35c3a`, whose sampled overlap weights close every receiver to 1). The FAS
transfers never read a chimera receiver: restriction, forcing and correction sums run over
interior cells only, the receivers are ghost cells refreshed by each level's own chimera fill,
and (with IRS) the ghost residual at a 102 face is already zero-gradient — the same policy IRS
applies to every non-connection ghost. Threshold 1e-8 at the case's `vnn = 0.25`
(the [campaign above](#six-configurations-on-four-cases) runs it at 1e-10 and `vnn = 0.5`):

| Check | Result |
|---|---|
| `cycle = none` vs pre-FAS binary, serial and 2 ranks | bit-identical (0 bytes differ) |
| `v` / `w`: 1 process vs 2 ranks (one block per rank, donors cross ranks on all levels) | bit-identical |
| Field vs `none` at 1e-8 | v: 4.0e-5 K, w: 4.1e-5 K |
| Admissibility fallbacks (`FAS:` reports) | none |
| Iterations (wall): none 9 153 (10.5 s ‡) | v 1 914 (2.0 s ‡; 4.8× / 5.2× ‡), **w 1 241 (0.74 s ‡; 7.4× / 14× ‡)** (iterations / wall, loaded node) |

**W is safe with chimera.** The study provisioned a V-only constraint near 102 boundaries
(CFL3D practice); the data shows W converging faster than V, bit-reproducibly across ranks and
with zero admissibility fallbacks, so no constraint is imposed.

## Stress campaign

Four hard configurations built for the FAS review (generators in the verification tooling), all
run to `res-threshold` 1e-8 with `coarse-iter 10` and immediate descent; fields checked against the `none`
field (≤ 1.5e-4 K, both runs ending at 1e-8 of their first residuals).

| Case | none | v | w | w + irs | w + irs gain (iterations / wall) |
|---|---:|---:|---:|---:|---:|
| 4 conformal blocks, 101 connections (2×2 × 32², 3 levels) | 9 681 (20 s ‡) | 183 | 69 | **41** (0.36 s ‡) | 236× / 57× ‡ |
| Localized q̇ hot spot (64², discontinuous source, 3 levels) | 28 736 (58 s ‡) | 507 | 142 | **40** (0.30 s ‡) | 718× / 189× ‡ |
| 3D box 32³ (3 levels) | 5 580 (108 s ‡) | 931 | 94 | **47** (2.0 s ‡) | 119× / 53× ‡ |
| 3D box 64³ (4 levels) | —¹ | — | 95 (21 s ‡) | **48** (15 s ‡) | — |

¹ The 64³ single grid was not run (the earlier campaign stopped it at 10 000 iterations and extrapolated ~33 000 to its
absolute threshold). The W-cycle iteration count barely moves from 32³ to 64³ (94 → 95, 47 → 48), the same N-independence
measured in 2D.

### The liner: a face defect, not a cycle defect

**Realistic liner (4 material layers, squeezed cells, hard BCs).** 64 × 128 cells, 0.3 × 3.9 mm
(aspect ratio 12.5), layers coat / liner / jacket / shell with k = 5 / 300 / 15 / 150 W/mK on
all-level grid lines, convective 303 walls (hot side Tref 3000 → 1500 K along the wall, cold
side h = 10⁴). The campaign version was one block; `test/numerics-features/multimat_liner` is
the same wall split along its length into two 64 × 64 blocks joined by a 101 connection, with
the material jumps inside both blocks (so 2 ranks each own a block).

On the code of the FAS review this case was the one that would not behave: the explicit
baseline died at the liner → jacket interface within a few steps at `vnn = 0.25` (T < 0 K),
plain IRS at `vnn = 1.0` or `0.5` locked into a bounded residual orbit with the IRS admissibility
fallback of that time (removed in October 2026, decision D2) firing on ~4·10⁵ cell updates per report, the cycles inherited the orbit whenever IRS
was on, and only `v` / `w` without IRS at `vnn = 0.1` converged (8 362 and 7 222 iterations,
against 1 259 270 for the baseline). **Both symptoms were one defect in the face flux, not in
the cycles.** An interior face took the conductivity of the lower-index cell's material: at the
300 → 15 jump the first jacket cell exchanged heat through 20 times its own k (the face flux
carried the copper-side transmissibility, 10.5 times the series value of the two half cells)
while its time step was built on k = 15, so every explicit update at that face violated the
stability bound — a crash without IRS, the admissibility orbit with it (the fallback was the
symptom; since its removal such an update stops the run). The steady field was wrong too: against a reference with each layer in its own block,
where the 101 connections always balanced the two half cells in series, the single-block
solution was off by 35.7 K.

The resolution is in the discretisation and needs no key: a face between two materials inside a
block carries the series conductivity of its two half cells, k = (h₁ + h₂) / (h₁/k₁ + h₂/k₂)
([Face Conductivity](../theory/numerics.md#face-conductivity)), and the local time step takes
the severest face conductivity of each direction, evaluated as the flux evaluates it
([VNN](../theory/time-integration.md#diffusive-stability-condition-vnn)); both are always on,
and blocks of one material keep their step bit for bit. Measured on the one-block liner at `res-threshold` 1e-8
(RK3, `cons`, 4 levels with a one-iteration descent, `cycle = none` for the rows without a
cycle; wall times ‡ the `Elapsed` line). The rows from `none, vnn = 0.6` to `w + IRS` were re-measured on 2026-10-29 with the
relative threshold; the first and the last row were measured when the face flux was fixed, against the per-layer reference
and with the old time step, and cannot be repeated with this code:

| Configuration | Result |
|---|---|
| Field vs the per-layer-block reference | 1.2e-8 K (was 35.7 K) |
| none, `vnn = 0.6` | converges: 153 677 iterations, 807 s ‡ |
| none, ceiling | `vnn = 0.625` stable (20 000 iterations), `vnn = 0.65` fails at iteration 91 (T = −43 K, outside the properties table): the RK3 limit is 0.628 ([VNN](../theory/time-integration.md#diffusive-stability-condition-vnn)) |
| IRS, `vnn = 1.0` | converges cleanly, no halved update, but slowly: 94 682 iterations (774 s ‡; with the absolute threshold it had not converged after 100 000) |
| w, `vnn = 0.5` / `0.6` | 808 iterations (12 s ‡) / 684 (9.8 s ‡) |
| w + IRS, `vnn = 0.5` | 1 267 iterations (25 s ‡), **zero** halved updates |
| series faces with the old cell time step, `vnn` 0.25 / 0.5 / 0.6 | stable over 3 000 iterations (residual 2.3 / 1.4 / 1.1, still descending): the face flux was the defect, the face step is what lets the run go at `vnn ≥ 0.5` |

The two-block `multimat_liner` converges in all six configurations on the current code
([campaign table](#six-configurations-on-four-cases): to a fall by 1e-10 none 244 325
iterations, IRS 123 712, the cycles 2 020 to 2 704, every run
without a single halved update). The liner is now a normal case, and the robustness ranking
the review had measured on it — cycles first, IRS last — was an artefact of the face flux.

## Anisotropic cells: the honest limit

Reference case `numerics-features/anisotropic_Plate` (128 × 16 cells on 1 m × 1.25 mm: cell
aspect ratio 100, 3 levels, exact linear solution), the campaign settings above (RK3, `cons`,
residual every iteration; the shipped case writes it every 100 and runs `w` + IRS):

| Configuration | Iterations (fall by 1e-10) | Wall (s) ‡ | Residual floor (iterations run)¹ |
|---|---:|---:|---:|
| none | 1 022 | 1.0 | 7.9e-13 (400 000) |
| IRS | 565 | 1.0 | 3.5e-12 (400 000) |
| v | 749 | 1.3 | 3.9e-12 (100 000) |
| w | 749 | 2.3 | 3.9e-12 (100 000) |
| v + IRS | 680 | 2.1 | 8.0e-13 (100 000) |
| w + IRS | 682 | 3.1 | 1.1e-12 (100 000) |

¹ The floor runs of the release before: an absolute threshold of 1e-13 that no run reached, so each ran its budget and the
floor is that of the arithmetic, not of the test (with the relative test the same runs stop at 1e-13 of their first residual,
1.8e-10 to 3.3e-10, above the floor).

V and W need the same iterations, IRS lowers each by about a tenth and IRS alone needs the fewest: the coarse
correction, not the cycle shape, is the bottleneck. Full coarsening with a pointwise explicit
smoother cannot represent the strong transverse coupling on the coarse levels. **In wall time ‡
the cycles cost more than the single grid here.** The threshold is reached early on this case
for a reason that has nothing to do with the solvers: the time step of an aspect-ratio-100 cell
is tiny, so the residual ‖ΔT‖ per iteration is small long before the field has converged, and at
1e-13 every configuration sits on a round-off floor of ~1e-12 (last column; all cross 1e-11
between 1 157 and 1 545 iterations and none reaches 1e-13 within its budget). On AR ≫ 1 meshes
run the single grid or IRS; the cure (line-implicit smoothing or semi-coarsening) is outside the
FAS v1 perimeter.

## Deep V instability

`cycle = v` without IRS on 4 or more levels can lock into a bounded limit cycle (residual
oscillating at a large value, a steady stream of `FAS: ... updates halved` reports). This is a
property of the method (weak explicit smoothing vs accumulated coarse corrections), not of a
case: use `w`, or enable IRS. Documented in
[time integration](../theory/time-integration.md#fas-cycles). The four-level liner above runs
`v` without IRS to a fall by 1e-10 in 2 020 iterations without a halved update: the instability
is possible, not certain.

## Re-validation at the sync of October 2026

The sync of October 2026 (MOSE's pins, parallel, input/stops/restart, metrics/dt and face-metric themes, then the removal
of the IRS admissibility fallback) was accepted against a 45-row FAS set run with `fuss-stress/fas_vv.sh` on the binary
before the sync (237c9c2, `fuss-stress/ref/237c9c2/fas_vv.md`) and on the integrated binary with the halved `vnn` values
(the same steps; `fuss-stress/results/fas-vv-s10/fas_vv.md`): the plates, chimera, liner, anisotropic and textbook cases
above, the 16³ and 32³ boxes, the k(T)/c_p(T) and bi-material plates, the stress cases and the guard probes, under `none`,
`v`, `w` and IRS. The set was run once more with the final binary of the sync (`verify/sync-2026-10` a616ef8, `USE_IPO`
on, ifx 2024.1; `fuss-stress/results/fas-vv-a616ef8/fas_vv.md`) and compared with `tools/fas_vv_summary.py compare`.

- **Against the binary before the sync (237c9c2):** 41 of the 45 rows reproduce their iteration counts and rates (38
  exactly, three within one iteration: `residual_convergence` `mg3` 212 413 → 212 414, `mg_irs` 79 817 → 79 816, the
  narrow-table `v` probe 942 → 941). The two probes built to fire the removed fallback (IRS on a table of 995–1705 K
  started at 1700 K, alone and under `v`) stop with exit code 1 instead of running, as the rule above says; every other
  row reports the same FAS guard events as before the sync (none, except the two probes under `v`: 4 corrections halved on
  the narrow table and 43 603 corrections and 181 143 coarse updates on the deep V without IRS). The 16³ box under `v` and
  `w` to the threshold 1e-12 takes 267 and 195 iterations instead of 261 and 186, the only rows outside the rule, and a
  round-off effect: the residual histories of the two binaries differ by at most 1.0e-12 K (absolute) while the residual
  is above 1e-9 (1.6e-7 and 3.2e-8 relative above 1e-6), they separate in the round-off tail (4.1e-12 against 4.5e-12 at
  iteration 250 on `v`) where the threshold 1e-12 lies, the single-grid run of the same box keeps its 2 983 iterations along
  an identical history, and the final fields differ by 5.7e-12 K (`v`) and 5.5e-12 K (`w`) (the last-bit changes of the cell
  volume and of the stencil metric of the sync).
- **Against the integrated binary without IPO (s10, 354bf28):** 45 rows, none failing, every iteration count
  and rate equal (the narrow-table `v` probe 939 → 941, with the same 4 halved corrections). 31 of the 45 residual histories
  are identical to the byte; in the other 14 the printed residuals differ in their last digits (by at most 1.0e-10 K, one unit
  of the tenth digit of a residual of order 1; 4e-5 relative in the tail near 1e-8), and the final fields differ by at most
  6.6e-12 K: the effect of IPO on the rounding (`CHANGELOG.md`), not of the code.
- **Ranks × threads and restart:** the bit identity under cycles (1×4, 2×1, 2×4, 4×1, 4×4 against 1×1 on `fas_w_irs`,
  chimera `w`, the liner `w` + IRS, the 128² plate `w` + IRS and the 16³ box `w`) holds on the final binary, 25 of 25 layouts,
  and the restart of a `v` run (302 + the rest against 535 in one go) reproduces the field to the last bit.
- **The relative convergence test (post-merge, 2026-10-29):** with `res-threshold` relative to the first residual of the
  finest level (as MOSE 3f3a4d9) the set was run again (`fuss-stress/results/fas-vv-post-i1b/fas_vv.md`, compared with the
  set above): every row that ends on its threshold ends 16 to 43 % earlier (the first residuals of the set are all above 1 K),
  the fields differ from the `none` ones by up to 4.4e-3 K accordingly, exit codes and guard reports are unchanged (the probe
  capped at `level1-iter` 5 000 ends its last cycle: 5 004 iterations), the bit identity holds on 25 of 25 layouts and the
  restart of the `v` run (302 + the rest against 448 in one go) is exact. The counts of this page are those of that run and of
  the 24 runs of the first table.

## Practical guidance

- **Recommended: `cycle = w` with IRS (`irs-beta = 0.5`, `vnn = 1.0`)**: the fewest iterations
  on the plate, the chimera case and every patch case of the review, textbook N-independent
  rates, insensitive to `coarse-iter` on deep hierarchies, safe with chimera and at any depth,
  and the configuration that stays stable where the others can fail (`v` on 4+ levels, the
  explicit single grid at the VNN ceiling). The default `cycle = none` is the bit-identical
  pre-FAS run, kept for that reason.
- **Wall time is not iterations.** `w` + IRS costs more per fine iteration than `v` or `w`
  alone; where the cycles converge within a few hundred iterations of each other the cheapest
  cycle per iteration wins the clock, and on the two-block liner plain `v` even needs the
  fewest iterations (`v` 17.8 s ‡, `w` + IRS 49.5 s ‡, both against
  1 266 s ‡ for the single grid). For a case that will be run many times, measure the
  candidates on it.
- **Strongly anisotropic cells (AR ≫ 1):** the cycles do not pay; the single grid or IRS alone
  is as fast ([above](#anisotropic-cells-the-honest-limit)).
- Multimaterial meshes with large conductivity jumps inside a block are ordinary cases: the
  series face conductivity and the face-conductivity time step are always on, no key to set;
  the liner runs every configuration at `vnn = 0.5` (`1.0` with IRS).
- `coarse-iter = 10` (the default) is right for deep hierarchies; raise it only on 2-level
  hierarchies where the coarsest solve is the bottleneck.
- `v` without IRS: 3 levels at most, unless measured on the case.
- Descend immediately (`level<m>-iter = 1` for m > 1): cycling from the start beats long
  sequencing stints. With `levels > 1` the stop is `level1-iter` or `res-threshold`:
  `iter-threshold` is ignored, with a warning at set-up.
