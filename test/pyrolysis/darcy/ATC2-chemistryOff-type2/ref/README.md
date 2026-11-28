# Reference data for the ATC2-chemistryOff type-2 gate (plan 07 section 4.4)

The same PATO 3.1 run as `../../../ATC2-chemistryOff/ref/` (plan 06 section 4.5: `tutorials/1D/AblationTestCase_2.x_chemistryOff`
with edits (1)-(3) of `../edit_pato_case.py`, zero-emissivity TACOT, `endTime 60`), repeated on 2026-10-08 with edit (4): `p`
added to the probed `fields` of `system/porousMat/plotDict`. PATO solves the type-2 model (`DarcyLawMassModel`,
`PyrolysisEnergyModel`), the closure FUSS's `gas-closure = type2` implements.

- `Ta_plot`, `Ta_surfacePatch`: byte-identical to the plan 06 run (md5 `7fdf47a7…`, `d253ab8e…`): adding a probed field changes nothing else.
- `p_plot`: the gas pressure [Pa], columns `t` then 1, 2, 4, 8, 12, 16, 24 mm below the heated face, every 0.1 s to 60 s. The
  heated face is open at 101325 Pa, the back face impermeable; the excursion `p - p_w` grows with the char layer and reaches
  2049 Pa at 24 mm at 60 s (371 Pa at 1 mm).

`verify.py` reads `$FUSS_PATO_ATC2_OUT` (a PATO `output/porousMat/scalar` directory with `Ta_plot` and `p_plot`) when set,
else this copy.
