# Regression baseline

`manifest.txt` is the md5 fingerprint of every artefact the default test suite produces
(`test/run_regression.sh`, no `-a`, one thread), recorded from a **fresh preset build**.

## Check against it

```bash
./test/run_regression.sh -p 1 -k          # runs the suite, diffs, exits 1 on any difference
```

## Re-record it — only when a result change is INTENDED

```bash
./test/run_regression.sh -B -p 1 -o /tmp/rec   # -B wipes build/ and rebuilds from the preset
cp /tmp/rec/manifest.txt test/baseline/manifest.txt
git add test/baseline/manifest.txt             # and say in the commit WHY results changed
```

Never record from an incremental build. On 2026-09-14 an incremental relink produced a binary
that differed from a fresh build of the same source on 90 of 290 artefacts; `-B` exists so this
cannot happen by forgetting.

## Why this file is committed

Every baseline recorded before 2026-09-14 lived in a session scratch directory and was lost when
that directory rotated. The cost was a 20-minute regeneration from the previous commit before the
next change could be judged inert. A baseline that is not in the repository is not a baseline.

## What "identical" means here

`SOLUTIONTIME` in the `.tec` headers is deliberately **not** stripped: for a steady run it carries
the converged iteration count, so a change in it is a real result change. Restart fidelity is
capped by the 16-digit ASCII solution format (~3e-15 relative), but that affects only
`run_restart_gate.sh`, not this manifest — every case here starts from its initial condition.

## Re-recordings (what moved, and why)

- **2026-10-07, ORION v1.7.0 (`5c32b73`), binary `7fc284fe`, 294 rows.** 271 rows changed, 23 identical. Every changed
  artefact is a Tecplot file written by ORION and differs in **line 1 only**: ORION `9ce4e12` (in v1.7.0) stopped quoting the
  `VARIABLES` names (`VARIABLES ="x" "y" "z"  "T" …` → `VARIABLES = x y z   T …`); every other byte is written by unchanged
  code. The 23 identical rows are the files ORION does not write (21 `residual-history.dat`, 2 probe files). Attribution
  evidence: a merge-only build (upstream `a418a26` at ORION `a770748`) reproduced the previous manifest on 294/294 rows, and
  the post-bump manifest equals the md5s predicted from the old artefacts with their first line rewritten
  (`plan-bucket/records/2026-10-07-submodule-bump/orion-delta.predicted-manifest.txt`). Row order is now the `C`-locale
  `sort` order (`run_regression.sh` exports `LC_ALL=C` since the same day; the old file was sorted under another collation,
  which moved `multimat_Plate`'s `field1/field10/field2` rows).
- **2026-10-07 (evening), plan 06 — new cases added; no pre-existing result changed. 480 rows.** The 294 rows of the previous
  record are byte-identical (fresh preset build of the final plan-06 source, binary `1c5d2934`, `run_regression.sh -B -o` on a
  scratch snapshot of the tree, `diff | grep '^<'` empty); the 186 added rows are `test/pyrolysis/kinetics_isothermal/T0700` (60),
  `T0900` (60), `T1200` (58) and `test/pyrolysis/energy_audit` (8). `ATC1` and `ATC2-chemistryOff` carry `.slow` and are not in the
  default sweep. The feature-off identity was also recorded after each of the three legacy-path edits separately (gates A/B/C in
  `plan-bucket/records/2026-10-07-plan06/`).
- **2026-10-08, verification gates — three cases added; no pre-existing result changed. 691 rows.** The 480 rows of the 2026-10-07
  record are byte-identical (`run_regression.sh -B -o` on a scratch snapshot of the tree, fresh preset build, `diff | grep '^<'` empty);
  the 211 added rows are `test/pyrolysis/adiabatic0d/T1000` (66), `test/pyrolysis/heated0d` (120) and
  `test/pyrolysis/conduction/graded100` (25), the exact-solution gates of `plan-bucket/records/2026-10-08-verification/` (no Fortran
  change; the binary is the 2026-10-07 source).
- **2026-10-09, plan 07 (the type-2 gas closure) — Darcy cases added; no pre-existing result changed. 914 rows.** The 691 rows
  of the 2026-10-08 record are byte-identical (fresh preset build of the final plan-07 source, binary `800a1206`,
  `run_regression.sh -B -k -o` on a scratch snapshot of the tree, `diff | grep '^<'` empty); the same identity held on an earlier
  snapshot of the work (binary `463c99a5`). The 223 added rows are `test/pyrolysis/darcy/heated0d_type2` (121),
  `adiabatic0d_type2` (67), `mass_balance_1d` (5), `sealed_nosource` (5), `sealed_profile/N020`, `N040`, `N080`, `N160` (5 each) and
  `tensor_linear_2d` (5). `darcy/ATC2-chemistryOff-type2` carries `.slow`. Between the two snapshots only the type-2 field files
  moved (the restart column `T_prev` became `TO`, record `plan-bucket/records/2026-10-09-plan07/`); probes, budgets and
  `darcy.txt` were identical.
- **2026-10-13, the benchmark fixes and gates: two type-2 cases moved at the 1e-10 level, six benchmark cases added. 1633 rows.**
  - **How measured:** each step ran `run_regression.sh -B -k -o` on a scratch snapshot (fresh preset build). The record is
    `plan-bucket/records/2026-10-13-benchmark-fixes/`.
  - **`061160c`** (BC 305 keeps h_conv, the η* limiter, the deficit supply; binary `7c6aefba`): 914/914 identical. No suite
    case uses BC 305, and neither the limiter nor the supply fires in any suite case.
  - **`54acf19`** (open faces carry h_g at the wall temperature, the quasi-steady carry) changed 68 rows in two cases:
    - `darcy/adiabatic0d_type2`, all 67 rows: T by at most 1.2e-10 relative, ṁ_pg and u_pg_1 by 1.2e-7 of their peak.
    - `darcy/heated0d_type2`: `energy-budget.txt` only, one entry, 4.4e-16 relative.
    - **Why they moved:** the open faces of both cases are adiabatic 301 faces, where T_w = T_cell, so the face enthalpy is
      the same. What changed is the form of the flux: h_w·½(Γ_c + Γ_w)·∇p replaces ½(h_c Γ_c + h_w Γ_w)·∇p. The Newton
      stop band of the temperature inversion carries last-bit changes in these cases to about 1e-10.
    - The suite's quasi-steady cases (`darcy/sealed_profile/N*`, `darcy/tensor_linear_2d`) did not move. They freeze T
      (solid cp 1e12), so η is constant and the carry is zero.
  - **`a10a2a3`** (cases and gates, no Fortran change) reproduced those 914 rows and added 719, each case running under 8 s:
    - `test/pyrolysis/exact/pressure_transient/N020` (53), `exact/sharp_front/N050` (113), `exact/transpiration/N020` (11);
    - `lit/baum_atreya2007/full_N100` (68), `qs_N100` (68);
    - `lit/zhang2010/N050` (406).
  - **Cross-checks:** the determinism gate's pass A (the worktree build of `a10a2a3`, binary `ac89d28c`) gives the same 1633 rows. The
    final source reproduces all 1633 (`1cc8fea` adds the restart reader, the later commits are test scripts).
  - **Binaries:** the final source's scratch build is `5a32f9bb` and its worktree build is `456af9c5`. They differ because the
    binary embeds its source paths, so an md5 identifies a build only within one tree.
- **2026-10-13 (night), plan 08 stage A (the ablating wall BC 306, static mesh): four type-1 cases moved by round-off, four
  ablation cases added. 1689 rows.**
  - **How measured:** `run_regression.sh -B -p 1 -k -o` on a scratch snapshot of the tree (fresh preset build, binary `7bf8ff24`;
    the worktree build of the same source is `15cece56`). 491 `<` lines, 547 `>` lines.
  - **Added, 56 rows:** `test/pyrolysis/ablation/demo_base` (5), `demo_base_pyro` (7), `seb_closure` (22), `zero_blowing` (22).
    `ATC21` and `atc2off_306_static(_t2)` carry `.slow`.
  - **Moved, 491 rows, all type-1 pyrolysis:** `adiabatic0d/T1000` (42), `energy_audit` (6), `heated0d` (80) and
    `lit/zhang2010/N050` (363). Maximum change, old binary (`d9d2c21`) against new, each field over all its files and
    scaled by the field's own maximum:
    - T 6.2e-15;
    - h and χ 1.9e-12 (`heated0d`: the decomposition carries the last bits through its 3 s);
    - ṁ_pg 7.4e-13 (`zhang2010`), other fields ≤ 1.6e-14;
    - energy budgets ≤ 1.9e-12 pointwise, probes ≤ 2.4e-15;
    - `residual-history.dat` 5.9e-10 (a difference of nearly equal numbers).
    `SOLUTIONTIME` moves in its last printed digit (dt depends on T). No type-2 or non-pyrolysis case moved.
  - **Why they moved:** LTO inlining under FMA contraction (`-march=native`). Bisected on `energy_audit` with scratch builds:
    - `d9d2c21` plus the appended, unused `Surface_Optics`/`Char_Enthalpy`: the old bits;
    - the full new source without the `Gas_State` call, or without the BC-kernel calls: still the new bits;
    - the full new source without the `Surface_Optics` and `Char_Enthalpy` calls: the old bits.
    The BC 306 code's new calls into `Tab_Lookup`/`Tau_Of` change how LTO inlines those kernels into the type-1 energy
    path, and the fused multiply-adds then round differently. No line of the type-1 path changed.
  - **Not done:** copying the two kernels to keep the bits would hold for one stage only (plan 08 stage C's remap calls
    `Solid_State`/`Energy_Of`). `-ffp-contract=off` would make results independent of inlining, at the cost of one global
    re-record; that is a project decision, left open.
  - **Gates:** all 23 physics gates pass on the new binary (`run_physics_gates.sh`, 8 min 55 s), the guard demos and the
    restart gate too.
- **2026-10-13 (night), plan 08 stage B (the receding mesh, law = nodeshift): one demo base changes on purpose, two
  cases added; no other row moved. 1700 rows.**
  - **How measured:** `run_regression.sh -B -p 1 -k -o` on a scratch snapshot of the stage-B tree (fresh preset build,
    binary `353f96f8`). 5 `<` lines, 16 `>` lines.
  - **Changed on purpose, 5 rows:** `test/pyrolysis/ablation/demo_base` now recedes (`law = nodeshift`,
    `sdot-override = 1e-5`), as plan 08 §2.10 specifies.
  - **Added, 11 rows:** `test/moving-mesh/conjugate_receding` (6) and `test/pyrolysis/ablation/demo_base_2blk` (5).
  - **The other 1684 rows are byte-identical.** Unlike stage A, these stage-B changes moved no type-1 bits.
- **2026-10-13 (night), plan 08 stage C (pyrolysis on the receding mesh): one case added; no other row moved. 1707 rows.**
  - **How measured:** `run_regression.sh -B -p 1 -k -o` on the stage-C tree (fresh preset build, binary `42dd54fa`).
    0 `<` lines, 7 `>` lines.
  - **Added, 7 rows:** `test/pyrolysis/ablation/remap_uniform` (gate V2). ATC22 and ATC23 are `.slow` and add no row.
  - **The other 1700 rows are byte-identical.** The T_prev shift moved out of `Copy_State` (phase 0) and the remap,
    budget and guard changes touch no pre-existing case.
- **2026-10-13 (night), plan 08 stage D (307, smoothing, multi-block pyrolysis): 15 type-2 cases move, three cases
  added. 1771 rows.**
  - **How measured:** `run_regression.sh -B -p 1 -k -o` on the stage-D tree (fresh preset build, binary `fed0baec`):
    506 `<` lines, 570 `>` lines. `test/pyrolysis/ablation/smoothing/{w1,w5,w5_2blk}` are `.slow` (their generator
    writes the markers).
  - **Moved, 506 rows, every one attributed to the Darcy line-solver change of stage D.**
    - The change: a face across the tridiag lines, and every boundary row, has its conductance on the line diagonal.
      An additive correction then removes the error that is uniform across the lines.
    - Why: with the whole cross flux on the right-hand side, a 2-D slab diverged.
    - The same fixed point by another arithmetic route. The moves are round-off on the line cases (T 1e-14 on
      `mass_balance_1d`, 9e-15 on `sharp_front`, 5e-13 on `seb_closure`) and up to 6e-9 (η_st) in the stiff 0-D cases.
    - The 15 cases: `darcy/{adiabatic0d_type2, heated0d_type2, mass_balance_1d, sealed_profile/N020..N160}`,
      `exact/{pressure_transient/N020, sharp_front/N050, transpiration/N020}`, `lit/baum_atreya2007/{full,qs}_N100`,
      `ablation/{seb_closure, zero_blowing, remap_uniform}`.
    - Proof, case by case: the stage-D sources with only the old splitting restored (without it the correction is
      inert) reproduce the stage-C binary byte for byte on all 15, every output file.
  - **Added, 64 rows:** `test/pyrolysis/ablation/coupled307/{run307,run306}` (56) and
    `test/moving-mesh/conjugate_receding_t2` (8).
  - **The other 1201 pre-existing rows are byte-identical.**

