# PATO oracle runs for ATC 2.1 (plan 08 §3.3, 2026-10-13, sprop2)

Each run is a copy of `$PATO_DIR/tutorials/1D/AblationTestCase_2.x` (PATO 3.1, OpenFOAM 7), edited by
`../../edit_pato_atc2x.py`, then `./Allrun`. Each directory below holds the run's `edits.diff` against the tutorial, as
`pato-edits.diff`.

| directory | edits | wall time | role |
|---|---|---|---|
| `.` (this one) | `--case 2.1 --max-dt 0.0005` | 331 s | the oracle of `verify.py` |
| `pato-dt0.005/` | `--case 2.1` (the tutorial's 5 ms step) | 81 s | the step study |
| `pato-dt0.00025-tight/` | `--case 2.1 --max-dt 0.00025 --ta-tol 1e-12` | 645 s | the step study |

Edit meanings:
- `--case 2.1`: `movingMesh no` and the preconditioners of Ta, p and Xsii set to DIC (the static mesh's symmetric matrices
  make the shipped DILU abort);
- `--max-dt`: `maxDeltaT` in `system/controlDict`;
- `--ta-tol`: the absolute tolerance of the Ta solve.

Files, every 0.1 s, as `output/porousMat/`:
- `Ta_plot`: T at 1, 2, 4, 8, 12, 16, 24 mm below the heated face;
- `Ta_surfacePatch`: the top face (T_w) and the bottom face;
- `mass`: `t ṁ_g ṁ_c virgin-98% char-2% recession`.

The oracle run's `controlDict` line reads `maxDeltaT   0.0005;`: its edit was a `sed` that dropped the line's trailing
comment, while the script keeps the comment. The two are the same setting.

Run health: no `FOAM FATAL`, and `mass` reaches 120 s in every run. Every Ta and p solve iterated at least once in the
0.5 ms and 0.25 ms-tight runs; the Xsi solves are diagonal, so they report 0 iterations by construction. The same
0.25 ms run at the shipped Ta tolerance (1e-6) did zero iterations in 30035 Ta solves between 40 and 60 s and drifted
38 K, so it is not kept.

`Amaryllis_{Energy,Mass}_TestCase_2.1.txt`: the organisers' reference, copied from the tutorial's `data/ref/Amaryllis/`
(header plus 1201 rows at 0.1 s).
