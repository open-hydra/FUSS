# Equilibrium B′ table for the Ablation Test Case 2.x cases (plan 08 §2.2, §3.4)

`bprime_tacot_earth.dat` holds TACOT pyrolysis gas and carbon char in Earth air at 101325 Pa, in PATO's 5-column layout
`p B'g B'c T h_w` (SI). It has 1 pressure × 30 B′g × 200 T = 6000 rows, with T uniform from 300 to 4000 K (18.59 K steps).
The file starts with the two header lines `FUSS_Lib_Bprime` reads: `// kind: equilibrium` and `// axes: p Bprime-g T`.
The `atc2x306` cases copy it into their `INPUT/` through `test/pyrolysis/tools/make_case.py --bprime`; this file is the default.
md5 `f38e9a36e0cfe94643d1304490e3a9f2`.

## How it was made (2026-10-13, sprop2, PATO 3.1 / OpenFOAM 7, PATO's mutation++)

```bash
set +u; source /opt/anaconda3/etc/profile.d/conda.sh && conda activate pato
. $PATO_DIR/src/applications/utilities/runFunctions/RunFunctions && pato_init
mkdir <s> && cp createBprimeTable.inputs <s>/ && cd <s>
createBprimeTable createBprimeTable.inputs > log.createBprimeTable 2>&1     # 3-45 s; writes p_101325/
python3 <this dir>/assemble_uniform_T.py <s> bprime_tacot_earth.dat         # 6000 rows
```

`createBprimeTable.inputs` follows plan 08 §3.4. The elemental MASS fractions are PATO's own `XE_TO_YE` conversions of the
ATC 2.x mole fractions, as the shipped ATC 2.2 run prints them:
- Earth air `O 0.21, N 0.79` gives `z_e[O] 0.232918`, `z_e[N] 0.767082`;
- TACOT gas `C 0.206, H 0.679, O 0.115` gives `z_pg[C] 0.494992`, `z_pg[H] 0.136918`, `z_pg[O] 0.36809`.

The 30 B′g values are those of PATO's Mars example. `assemble_uniform_T.py` writes the run's uniform-T scans
`p_101325/{Bc,hw}_Bg<v>` in `BgList` order, rows sorted by T. A rerun of the script on the run above reproduces this file byte
for byte.

## Why the uniform-T assembly and not createBprimeTable's own table

createBprimeTable's own table places its rows at log-uniform B′c (`numberBc` per sub-table). Across the oxidation plateau
one sub-block therefore jumps from 1093 to 2187 K between two rows. PATO interpolates linearly in T, so it reads a wrong
`h_w` there.

PATO ATC 2.2 was run three times with `BprimeFile` in place of `mixtureMutationBprime tacot26;` (PATO refuses both keys
together) and compared with the shipped on-the-fly run. Probes were compared only while the receding surface is above
them in both runs.

| table | rows | max \|ΔT_w\| | max \|ΔT probe\| | recession, relative (t ≥ 5 s) |
|---|---|---|---|---|
| createBprimeTable, `numberBc 100`, `nTemp 200` | 2943 | 15.69 K | 14.66 K | 1.9e-2 |
| createBprimeTable, `numberBc 2000`, `nTemp 1000` | 59586 | 0.26 K | 0.28 K | 2.5e-4 |
| **this file** (uniform-T assembly of the first run) | 6000 | **0.16 K** | **0.14 K** | **3.0e-5** |

The `numberBc 2000` table passes as well, but it is ten times larger. It also prints T with 6 significant digits, which
leaves 3633 pairs of equal consecutive T inside ascending sub-blocks. `FUSS_Lib_Bprime` accepts such ties with PATO's
first-match rule.

## Cross-check against mutation++ `bprime` standalone (plan 08 §3.4)

`bprime` reads the edge and pyrolysis compositions by name from the mixture file. `tacot26.xml` names none, so the check
uses `tacot26_bp.xml` (this directory), which carries the same 26 species and the ATC 2.x mole fractions:

```bash
bprime -T <T> -P 101325 -b <B'g> -m tacot26_bp -bl BLedge -py Gas     # in a directory holding tacot26_bp.xml
```

| T [K] | B′g | B′c table | B′c bprime | rel. Δ | h_w table [J/kg] | h_w bprime [J/kg] | Δ [J/kg] |
|---|---|---|---|---|---|---|---|
| 671.859 | 0.001 | 0.0875435 | 0.0875433 | 2.3e-6 | −2.23826e6 | −2.23825e6 | 10 |
| 1415.58 | 0.1 | 0.152956 | 0.152956 | 0 | −48368.8 | −48363.7 | 5 |
| 2345.23 | 0.8 | 0.0320505 | 0.03205 | 1.6e-5 | 2.5815e6 | 2.58151e6 | 10 |
| 3088.94 | 2.15 | 0.131062 | 0.131056 | 4.6e-5 | 9.19336e6 | 9.1933e6 | 60 |
| 3832.66 | 20 | 69.3437 | 69.341 | 3.9e-5 | 2.91373e7 | 2.91372e7 | 100 |

The bounds are `|ΔB′c|/B′c ≤ 1e-3` and `|Δh_w| ≤ 1e3 J/kg`; all five rows pass. The differences come from the printed digits
(6 in the table, as `bprime` prints).
