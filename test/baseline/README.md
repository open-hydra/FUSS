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
