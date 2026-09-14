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
