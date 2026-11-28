# PATO oracle runs for ATC 2.2 (plan 08 §3.3, 2026-10-13, sprop2)

Each run is a copy of `$PATO_DIR/tutorials/1D/AblationTestCase_2.x` (PATO 3.1, OpenFOAM 7), edited by
`../../edit_pato_atc2x.py`, then `./Allrun`. Each directory holds the run's `edits.diff` against the tutorial, as
`pato-edits.diff`.

| directory | edits | wall time | role |
|---|---|---|---|
| `.` (this one) | `--case 2.2 --max-dt 0.0005` | 600 s | the oracle of `verify.py` |
| `pato-dt0.005/` | `--case 2.2` (the tutorial's 5 ms step) | 113 s | the step study (first-order Richardson with the oracle) |
| `mesh-N400/` | `--case 2.2 --max-dt 0.0005 --ncells 400 --end-time 8` | 103 s | the mesh study of the first 8 s (the same run at 100 cells and 8 s reproduces the oracle's first 8 s bit for bit) |

Files, every 0.1 s, as `output/porousMat/`:
- `Ta_plot`: T at 1, 2, 4, 8, 12, 16, 24 mm below the initial surface (0 once the surface has passed the point);
- `Ta_surfacePatch`: the top face (T_w) and the bottom face;
- `mass`: `t ṁ_g ṁ_c virgin-98% char-2% recession`.

Run health: no `FOAM FATAL`, and every Ta and p solve iterated at least once (`solver.health`).

`Amaryllis_{Energy,Mass}_TestCase_2.2.txt`: the organisers' reference, copied from the tutorial's `data/ref/Amaryllis/`
(header plus 1201 rows at 0.1 s). A probe the surface has passed reads 0, or a stray value below 1 K (8 mm at 40.7 s
of 2.2 reads 0.01191).
