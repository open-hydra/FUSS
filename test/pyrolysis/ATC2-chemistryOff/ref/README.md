# Reference data for the ATC2-chemistryOff gate (plan 06 section 4.5)

`Ta_plot` and `Ta_surfacePatch`: PATO 3.1 `tutorials/1D/AblationTestCase_2.x_chemistryOff` run on 2026-10-07 with the
edits of `../edit_pato_case.py` (every `BoundaryConditions` row at `hconv = 300`, `Tedge = 2800`, `chemistryOn = 0`;
`MaterialPropertiesDirectory` pointing at a copy of the TACOT tables whose emissivity column is 0; `endTime 60`), so
that PATO's surface energy balance is the closed form `hconv (Tedge - T_w)` that FUSS's BC 303 evaluates. `Ta_plot`
columns: `t` then 1, 2, 4, 8, 12, 16, 24 mm below the heated face every 0.1 s; `Ta_surfacePatch`: `t`, surface, back
face. `verify.py` reads `$FUSS_PATO_ATC2_OUT` (a PATO `output/porousMat/scalar` directory) when set, else this copy.
