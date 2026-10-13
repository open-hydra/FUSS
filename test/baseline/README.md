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
