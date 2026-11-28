# qvol_Manufactured_Plate

Manufactured-solution gate for the volumetric-source reader (plan 04 item 2; fix `19db31a`).
The physics and tolerances are documented in `verify.py`; this file records the numbers the
plan asked for and the one-line check that the case is capable of seeing the defect.

## Observed error before and after the fix

| binary | L-inf [K] | L-2 [K] | source read-back L-inf [W/m3] |
|---|---|---|---|
| `2237fdb` (reads `vars(3)` of a `(1:1)` array — out-of-bounds) | **4.8501** | — | **48310** |
| `19db31a` and later (reads `vars(1)`) | 0.0426 | 0.0213 | 0.000 |

Budget: `TOL_LINF = 1.0 K`, `TOL_L2 = 0.5 K`, `TOL_QVOL = 1.0 W/m3` (see `verify.py` for the derivation).

## Why this case can see the defect and `qvol_Plate` cannot

A uniform source field cannot expose an indexing error in the reader. This case's `INPUT/st.tec`
carries 2398 distinct values over 2400 cells (`qvol_Plate`: 1 distinct value). Check it with

```bash
python3 -c "v=open('INPUT/st.tec').read().split('\n',2)[2].split(); print(len(set(v[-2400:])), 'distinct source values of 2400')"
```

(BLOCK packing puts the cell-centred `qvol` column last, after the 3 x 61 x 41 x 2 nodal coordinates).

## Run

```bash
./FUSS.sh -p 1 solve && python3 verify.py
```
